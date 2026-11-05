#!/usr/bin/env python3
"""Final report: turn anchor + quant results into the deliverable tables.

Outputs (all under `<out>/pri/`):

* `pri_mirna.tsv`            one row per pri-miRNA locus - the annotation itself
* `pri_pairs.tsv`            one row per hairpin x transcript pair - the evidence
* `pri_expression_tpm.tsv`   copied from the quant stage
* `pri_expression_counts.tsv`
* `manifest.json`            inputs, parameters, tool versions, stage stats

`genomic_env` and `class_code` are optional: without a reference GTF / SQANTI3
classification they are written as `unknown` / `NA` instead of failing.
"""

from __future__ import annotations

import bisect
import json
import shutil
import sys
import time
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path

from .. import __schema_version__, __version__
from ..config import RunConfig
from ..utils import env_with_env_bin
from .anchor import TX_ID_RE, parse_gtf

DETECTED_TPM = 1.0     # a sample counts as "detected" above this pri_sum TPM


def _status(msg: str) -> None:
    sys.stderr.write(f"[report] {msg}\n")


# --------------------------------------------------------------------------
# readers
# --------------------------------------------------------------------------
def read_tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]


def read_expression(path: Path) -> tuple[list[str], dict[str, dict[str, float]]]:
    """Returns (sample columns, locus -> {column: value})."""
    if not path.exists():
        return [], {}
    with path.open(encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cols = header[1:]
        rows: dict[str, dict[str, float]] = {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            row = {}
            for name, value in zip(cols, f[1:]):
                try:
                    row[name] = float(value)
                except ValueError:
                    row[name] = float("nan")
            rows[f[0]] = row
    return cols, rows


def sample_columns(cols: list[str]) -> list[str]:
    return [c for c in cols if not c.startswith(("mean_", "pri_max"))]


# --------------------------------------------------------------------------
# genomic environment (relative to the reference annotation)
# --------------------------------------------------------------------------
def build_exon_index(gtf: Path) -> dict[tuple[str, str], tuple[list[int], list[int]]]:
    """(chrom, strand) -> (sorted starts, sorted ends) of annotated exons."""
    exons: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    spans: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    tx_span: dict[str, tuple[str, str, int, int]] = {}
    with gtf.open(encoding="utf-8-sig", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9:
                continue
            m = TX_ID_RE.search(f[8])
            if not m:
                continue
            start, end = int(f[3]) - 1, int(f[4])
            if f[2] == "exon":
                exons[(f[0], f[6])].append((start, end))
            tid = m.group(1)
            chrom, strand = f[0], f[6]
            if tid in tx_span:
                c, s, st, en = tx_span[tid]
                tx_span[tid] = (c, s, min(st, start), max(en, end))
            else:
                tx_span[tid] = (chrom, strand, start, end)
    for chrom, strand, start, end in tx_span.values():
        spans[(chrom, strand)].append((start, end))

    def pack(d):
        return {k: ([s for s, _ in sorted(v)], sorted(e for _, e in v)) for k, v in d.items()}

    return {"exons": pack(exons), "genes": pack(spans)}   # type: ignore[return-value]


def classify_env(hp_chrom: str, hp_start: int, hp_end: int, hp_strand: str,
                 idx: dict) -> str:
    for kind, label in (("exons", "genic-exonic"), ("genes", "genic-intronic")):
        entry = idx[kind].get((hp_chrom, hp_strand))
        if not entry:
            continue
        starts, ends = entry
        i = bisect.bisect_right(starts, hp_end)
        for j in range(max(0, i - 50), i):
            if ends[j] >= hp_start:
                return label
    return "intergenic"


def read_classification(catalog_dir: Path) -> dict[str, str]:
    """SQANTI3 `*classification.txt`: isoform -> structural_category."""
    hits = sorted(catalog_dir.glob("*classification.txt"))
    if not hits:
        return {}
    out: dict[str, str] = {}
    with hits[0].open(encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if "isoform" not in header or "structural_category" not in header:
            return {}
        i_id, i_cat = header.index("isoform"), header.index("structural_category")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > i_cat:
                out[f[i_id]] = f[i_cat]
    return out


def read_fl_counts(catalog_dir: Path) -> dict[str, int]:
    p = catalog_dir / "fl_counts.tsv"
    if not p.exists():
        return {}
    out: dict[str, int] = {}
    for row in read_tsv(p):
        try:
            out[row["transcript_id"]] = int(float(row["n_fl"]))
        except (KeyError, ValueError):
            continue
    return out


# --------------------------------------------------------------------------
# assembly of the final tables
# --------------------------------------------------------------------------
def build_locus_rows(cfg: RunConfig, pairs: list[dict], intermed: list[dict],
                     expr_cols: list[str], expr: dict[str, dict[str, float]],
                     counts: dict[str, dict[str, float]],
                     env_idx: dict | None, classes: dict[str, str],
                     fl: dict[str, int], backend: str,
                     chains: dict[str, tuple]) -> list[dict]:
    samples = sample_columns(expr_cols)
    groups = sorted({s["group"] for s in _samples_of(cfg) if s["group"]}) if cfg.rna else []

    by_locus: dict[str, list[dict]] = defaultdict(list)
    for p in pairs:
        by_locus[p["hairpin_id"]].append(p)

    rows: list[dict] = []
    for locus, items in by_locus.items():
        first = items[0]
        n_spliced = sum(1 for i in items if i["anchor_type"] == "spliced")
        n_intronic = sum(1 for i in items if i["anchor_type"] == "intronic")
        n_partial = sum(1 for i in items if i["anchor_type"] == "partial")
        sources = {s for i in items for s in i["transcript_source"].split("+")}
        if "Iso-Seq" in sources:
            evidence = "L1_iso"
        elif "RNA-Seq" in sources:
            evidence = "L2_rnaseq"
        else:
            evidence = "L3_annotation"

        row = {
            "schema_version": __schema_version__,
            "locus_id": locus,
            "chrom": first["chrom"],
            "hp_start": int(first["hp_start"]),
            "hp_end": int(first["hp_end"]),
            "strand": first["strand"],
            "hp_length": int(first["hp_length"]),
            "n_transcripts": len({i["transcript_id"] for i in items}),
            "n_structures": len({chains[t] for t in {i["transcript_id"] for i in items}
                                 if t in chains}) if chains else "NA",
            "n_spliced": n_spliced,
            "n_intronic": n_intronic,
            "n_partial": n_partial,
            "anchor_type": _dominant_type(n_spliced, n_intronic, n_partial),
            "transcripts": ",".join(sorted({i["transcript_id"] for i in items})),
            "sources": "+".join(sorted(sources)),
            "evidence_level": evidence,
            "flank_5p_max": max(int(i["flank_5p"]) for i in items),
            "flank_3p_max": max(int(i["flank_3p"]) for i in items),
        }

        if env_idx is not None:
            row["genomic_env"] = classify_env(row["chrom"], row["hp_start"],
                                              row["hp_end"], row["strand"], env_idx)
        else:
            row["genomic_env"] = "unknown"

        cats = {classes.get(i["transcript_id"], "NA") for i in items}
        row["class_code"] = "/".join(sorted(cats))

        n_int = sum(1 for i in intermed if i["hairpin_id"] == locus)
        row["n_intermediates"] = n_int
        row["n_fl"] = sum(fl.get(i["transcript_id"], 0) for i in items) if fl else "NA"

        vals = expr.get(locus, {})
        if vals:
            per_sample = [vals.get(s, 0.0) for s in samples]
            row["expr_backend"] = backend
            row["pri_sum_tpm"] = round(sum(per_sample) / len(per_sample), 4) if per_sample else "NA"
            row["pri_sum_tpm_max"] = round(max(per_sample), 4) if per_sample else "NA"
            row["pri_sum_tpm_max_sample"] = samples[per_sample.index(max(per_sample))] if per_sample else "NA"
            row["n_samples_detected"] = sum(1 for v in per_sample if v >= DETECTED_TPM)
            for g in groups:
                gv = [vals.get(s["sample_id"], 0.0) for s in _samples_of(cfg) if s["group"] == g]
                row[f"tpm_{g}"] = round(sum(gv) / len(gv), 4) if gv else "NA"
            cvals = counts.get(locus, {})
            if cvals and backend == "kallisto":
                cs = [cvals.get(s, 0.0) for s in samples]
                row["pri_sum_counts"] = round(sum(cs) / len(cs), 2)
            else:
                row["pri_sum_counts"] = "NA"
        else:
            row["expr_backend"] = "none"
            row["pri_sum_tpm"] = "NA"
            row["pri_sum_tpm_max"] = "NA"
            row["pri_sum_tpm_max_sample"] = "NA"
            row["n_samples_detected"] = "NA"
            for g in groups:
                row[f"tpm_{g}"] = "NA"
            row["pri_sum_counts"] = "NA"

        rows.append(row)
    rows.sort(key=lambda r: (r["chrom"], r["hp_start"]))
    return rows


def _dominant_type(n_spliced: int, n_intronic: int, n_partial: int) -> str:
    pairs = [("intronic", n_intronic), ("spliced", n_spliced), ("partial", n_partial)]
    pairs.sort(key=lambda x: -x[1])
    return pairs[0][0] if pairs[0][1] else "NA"


def _samples_of(cfg: RunConfig) -> list[dict]:
    from .quant import load_samples

    return load_samples(cfg)


PAIR_COLS = ["locus_id", "transcript_id", "transcript_source", "tx_start", "tx_end",
             "n_exons", "anchor_type", "exon_index", "flank_5p", "flank_3p"]


def write_table(path: Path, rows: list[dict], cols: list[str] | None = None) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    cols = cols or list(rows[0].keys())
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "NA")) for c in cols) + "\n")


# --------------------------------------------------------------------------
# manifest
# --------------------------------------------------------------------------
def build_manifest(cfg: RunConfig, rep, anchor_stats: dict, quant_stats: dict,
                   n_loci: int, outputs: list[Path]) -> dict:
    params = {k: (str(v) if isinstance(v, Path) else v)
              for k, v in asdict(cfg).items() if k != "iso"}
    return {
        "tool": {"name": "pripipe", "version": __version__,
                 "schema_version": __schema_version__},
        "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "mode": cfg.mode,
        "params": params,
        "inputs": {name: {"path": str(p), "size": p.stat().st_size}
                   for name, p in (("genome", cfg.genome), ("mirna_bed", cfg.mirna_bed),
                                   ("gtf", cfg.gtf), ("rna", cfg.rna))
                   if p is not None and Path(p).exists()},
        "iso_inputs": [{"path": str(p), "size": p.stat().st_size}
                       for p in cfg.iso if p.exists()],
        "validation": {"warnings": rep.warnings, "notes": rep.notes},
        "stages": {"anchor": anchor_stats, "quant": quant_stats},
        "n_loci_reported": n_loci,
        "outputs": [str(p) for p in outputs],
    }


# --------------------------------------------------------------------------
# stage entry
# --------------------------------------------------------------------------
def stage_report(cfg: RunConfig, rep) -> None:
    anchor_dir = cfg.stage_dir("anchor")
    quant_dir = cfg.stage_dir("quant")
    pri_dir = cfg.out / "pri"
    pri_dir.mkdir(parents=True, exist_ok=True)

    pairs = read_tsv(anchor_dir / "anchor_pairs.tsv")
    intermed = read_tsv(anchor_dir / "intermediates.tsv")
    if not pairs:
        _status("no anchor pairs found - nothing to report")
        raise SystemExit(5)

    expr_cols, expr = read_expression(quant_dir / "pri_expression_tpm.tsv")
    _, counts = read_expression(quant_dir / "pri_expression_counts.tsv")
    quant_stats = {}
    qs = quant_dir / "quant_stats.json"
    if qs.exists():
        quant_stats = json.loads(qs.read_text(encoding="utf-8"))
    if quant_dir.exists():
        for name in ("pri_expression_tpm.tsv", "pri_expression_counts.tsv"):
            src = quant_dir / name
            if src.exists():
                shutil.copy2(src, pri_dir / name)

    env_idx = None
    if cfg.gtf and cfg.gtf.exists():
        env_idx = build_exon_index(cfg.gtf)
    else:
        _status("no GTF: genomic_env = unknown")

    classes = read_classification(cfg.stage_dir("catalog"))
    fl = read_fl_counts(cfg.stage_dir("catalog"))
    backend = quant_stats.get("backend", "none") if quant_stats else "none"

    chains: dict[str, tuple] = {}
    catalog_gtf = cfg.stage_dir("catalog") / "transcripts.gtf"
    if catalog_gtf.exists():
        for tid, tx in parse_gtf(catalog_gtf).items():
            chains[tid] = (tx.chrom, tx.strand, tuple(tx.exons))

    rows = build_locus_rows(cfg, pairs, intermed, expr_cols, expr, counts,
                            env_idx, classes, fl, backend, chains)
    write_table(pri_dir / "pri_mirna.tsv", rows)

    pair_rows = [{"locus_id": p["hairpin_id"], **{k: p[k] for k in PAIR_COLS if k != "locus_id"}}
                 for p in pairs]
    write_table(pri_dir / "pri_pairs.tsv", pair_rows, PAIR_COLS)

    anchor_stats = {}
    ajs = anchor_dir / "anchor_stats.json"
    if ajs.exists():
        anchor_stats = json.loads(ajs.read_text(encoding="utf-8"))

    outputs = sorted(pri_dir.glob("*"))

    qc_dir = cfg.out / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    bam_qc = rep.facts.get("bam_qc", [])
    if bam_qc:
        write_table(qc_dir / "flagstat_summary.tsv", bam_qc,
                    ["file", "total", "mapped", "unmapped", "duplicates",
                     "sampled", "spliced_sampled", "mapq_max"])
    (qc_dir / "validation.json").write_text(rep.to_json(), encoding="utf-8")
    outputs = sorted(pri_dir.glob("*")) + sorted(qc_dir.glob("*"))

    manifest = build_manifest(cfg, rep, anchor_stats, quant_stats, len(rows), outputs)
    (cfg.out / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False),
                                           encoding="utf-8")

    rep.facts["report"] = {"n_loci": len(rows), "n_pairs": len(pair_rows)}
    _status(f"{len(rows)} loci / {len(pair_rows)} pairs -> {pri_dir}")
