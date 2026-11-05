#!/usr/bin/env python3
"""Quantification of pri-miRNA loci.

Two backends, chosen by `--quant-backend auto`:

* ``kallisto``   - pseudo-alignment, needs FASTQ, gives EM-split counts + TPM.
* ``stringtie``  - ``-e -B`` on the existing BAM, gives coverage-based TPM.

The catalog is quantified as a whole, then aggregated per pri-miRNA locus
(sum / max over its supporting transcripts). Two matrices are produced when
an Iso-Seq-only catalog is available: the union catalog and the Iso-Seq-only
catalog, so RNA-Seq-derived pri loci can be spotted for the circularity bias.
"""

from __future__ import annotations

import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from ..config import RunConfig, _read_sheet
from ..utils import env_with_env_bin, missing_tool_message

SE_FRAG_LEN = 200
SE_FRAG_SD = 20
KALLISTO_K = 31


def _status(msg: str) -> None:
    sys.stderr.write(f"[quant] {msg}\n")


def _run(cmd: list[str], log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("a", encoding="utf-8") as fh:
        fh.write("\n$ " + " ".join(cmd) + "\n")
        fh.flush()
        try:
            proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                                  env=env_with_env_bin())
        except FileNotFoundError:
            raise RuntimeError(missing_tool_message(cmd[0])) from None
    if proc.returncode != 0:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(cmd)} (see {log})")


# --------------------------------------------------------------------------
# samples
# --------------------------------------------------------------------------
def load_samples(cfg: RunConfig) -> list[dict]:
    if not cfg.rna:
        return []
    header, rows = _read_sheet(cfg.rna)
    idx = {k: i for i, k in enumerate(header)}
    out: list[dict] = []
    for r in rows:
        def get(key: str) -> str | None:
            i = idx.get(key)
            if i is None or i >= len(r):
                return None
            v = r[i].strip()
            return None if v in ("", "-", ".") else v

        out.append({"sample_id": get("sample_id"), "r1": get("r1"), "r2": get("r2"),
                    "bam": get("bam"), "group": get("group")})
    return out


# --------------------------------------------------------------------------
# catalog fasta (kallisto needs it)
# --------------------------------------------------------------------------
def ensure_catalog_fasta(cfg: RunConfig, gtf: Path) -> Path:
    fa = gtf.with_suffix(".fa")
    if fa.exists() and fa.stat().st_size > 0:
        return fa
    _status(f"building transcript FASTA with gffread: {fa.name}")
    _run(["gffread", "-w", str(fa), "-g", str(cfg.genome), str(gtf)],
         cfg.stage_dir("quant") / "logs" / "gffread.log")
    return fa


# --------------------------------------------------------------------------
# kallisto
# --------------------------------------------------------------------------
def kallisto_index(fa: Path, idx: Path, log: Path) -> None:
    if idx.exists():
        return
    _run(["kallisto", "index", "-i", str(idx), "-k", str(KALLISTO_K), str(fa)], log)


def kallisto_sample(sample: dict, idx: Path, outdir: Path, log: Path, threads: int) -> dict[str, tuple[float, float]]:
    ab = outdir / "abundance.tsv"
    if not ab.exists():
        outdir.parent.mkdir(parents=True, exist_ok=True)   # kallisto only makes one level
        cmd = ["kallisto", "quant", "-i", str(idx), "-o", str(outdir),
               "-t", str(threads), "--plaintext"]
        if sample["r2"]:
            cmd += [sample["r1"], sample["r2"]]
        else:
            cmd += ["--single", "-l", str(SE_FRAG_LEN), "-s", str(SE_FRAG_SD), sample["r1"]]
        _run(cmd, log)
    vals: dict[str, tuple[float, float]] = {}
    with ab.open() as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            vals[f[0]] = (float(f[3]), float(f[4]))   # est_counts, tpm
    return vals


# --------------------------------------------------------------------------
# stringtie -e -B (BAM only)
# --------------------------------------------------------------------------
def stringtie_sample(sample: dict, gtf: Path, outdir: Path, log: Path, threads: int) -> dict[str, tuple[float, float, float]]:
    """Returns tid -> (count=nan, tpm, cov).

    StringTie 3.x Ballgown tables carry `cov` and `FPKM` but no `TPM`, so TPM
    is derived as FPKM / sum(FPKM) * 1e6 over the quantified catalog.
    """
    tdata = outdir / "t_data.ctab"
    if not tdata.exists():
        outdir.mkdir(parents=True, exist_ok=True)
        _run(["stringtie", "-e", "-B", "-G", str(gtf), "-p", str(threads),
              "-o", str(outdir / "quant.gtf"), sample["bam"]], log)
    rows: list[tuple[str, float, float]] = []
    with tdata.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        id_col = "t_name" if "t_name" in col else "transcript_id"
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append((f[col[id_col]], float(f[col["FPKM"]]), float(f[col["cov"]])))
    total = sum(fpkm for _, fpkm, _ in rows) or 1.0
    return {tid: (float("nan"), fpkm / total * 1e6, cov) for tid, fpkm, cov in rows}


# --------------------------------------------------------------------------
# aggregation
# --------------------------------------------------------------------------
def load_locus_map(pairs_tsv: Path) -> dict[str, set[str]]:
    """locus -> set(transcript ids) from the anchor table."""
    out: dict[str, set[str]] = defaultdict(set)
    with pairs_tsv.open(encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        c_hp, c_tx = header.index("hairpin_id"), header.index("transcript_id")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[c_hp]].add(f[c_tx])
    return dict(out)


def aggregate(locus_map: dict[str, set[str]],
              per_sample: dict[str, dict[str, tuple[float, float, float]]],
              idx: int = 1) -> dict[str, dict[str, float]]:
    """locus -> {sample: sum over supporting transcripts}; idx 1=TPM, 0=counts."""
    rows: dict[str, dict[str, float]] = {}
    for locus, tids in locus_map.items():
        row: dict[str, float] = {"n_isoforms": len(tids)}
        for sample, vals in per_sample.items():
            v = [vals[t][idx] for t in tids if t in vals]
            row[sample] = sum(v) if v else 0.0
            row[f"_max_{sample}"] = max(v) if v else 0.0
        rows[locus] = row
    return rows


def _sample_columns(samples: list[dict]) -> list[str]:
    return [s["sample_id"] for s in samples]


def write_matrix(path: Path, rows: dict[str, dict[str, float]], samples: list[dict],
                 value: str, na_counts: bool = False) -> None:
    def fmt(v: float) -> str:
        if na_counts or v != v:          # nan check
            return "NA"
        return f"{v:.4f}"

    ids = _sample_columns(samples)
    groups = sorted({s["group"] for s in samples if s["group"]})
    cols = ["hairpin_id", *ids, *[f"mean_{g}" for g in groups], "mean_all"]
    if value == "sum":
        cols.append("pri_max_mean")
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for locus, row in sorted(rows.items()):
            vals = [row.get(s, 0.0) for s in ids]
            line = [locus] + [fmt(v) for v in vals]
            for g in groups:
                gv = [row[s["sample_id"]] for s in samples
                      if s["group"] == g and s["sample_id"] in row]
                line.append(fmt(sum(gv) / len(gv)) if gv else "NA")
            line.append(fmt(sum(vals) / len(vals)) if vals else "NA")
            if value == "sum":
                line.append(fmt(max(row.get(f'_max_{s}', 0.0) for s in ids)))
            fh.write("\t".join(line) + "\n")


# --------------------------------------------------------------------------
# stage entry
# --------------------------------------------------------------------------
def stage_quant(cfg: RunConfig, rep) -> None:
    if cfg.mode == "iso":
        _status("no RNA-Seq input: expression matrix skipped")
        rep.facts["quant"] = {"skipped": "iso-only run"}
        return

    samp_dir = cfg.stage_dir("quant")
    logs = samp_dir / "logs"
    samples = [s for s in load_samples(cfg) if s["sample_id"]]
    if not samples:
        raise SystemExit("quant: no samples in --rna sheet")

    backend = cfg.resolve_quant_backend()
    if backend == "kallisto" and not any(s["r1"] for s in samples):
        _status("no FASTQ in sample sheet, falling back to stringtie")
        backend = "stringtie"
    if backend == "stringtie" and not any(s["bam"] for s in samples):
        raise SystemExit("quant: stringtie backend needs BAM, but no bam column in the sheet")

    gtf = cfg.stage_dir("catalog") / "transcripts.gtf"
    if not gtf.exists():
        raise SystemExit(f"quant: catalog not found: {gtf}")

    per_sample: dict[str, dict[str, tuple[float, float, float]]] = {}
    if backend == "kallisto":
        fa = ensure_catalog_fasta(cfg, gtf)
        idx = samp_dir / "index" / "catalog.idx"
        idx.parent.mkdir(parents=True, exist_ok=True)
        kallisto_index(fa, idx, logs / "kallisto_index.log")
        for s in samples:
            _status(f"kallisto: {s['sample_id']}")
            vals = kallisto_sample(s, idx, samp_dir / "kallisto" / s["sample_id"],
                                   logs / f"{s['sample_id']}.log", cfg.threads)
            per_sample[s["sample_id"]] = {t: (c, tpm, float("nan")) for t, (c, tpm) in vals.items()}
    else:
        for s in samples:
            if not s["bam"]:
                _status(f"{s['sample_id']}: no BAM, skipped")
                continue
            _status(f"stringtie -e -B: {s['sample_id']}")
            per_sample[s["sample_id"]] = stringtie_sample(
                s, gtf, samp_dir / "stringtie" / s["sample_id"],
                logs / f"{s['sample_id']}.log", cfg.threads)

    pairs_tsv = cfg.stage_dir("anchor") / "anchor_pairs.tsv"
    locus_map = load_locus_map(pairs_tsv)
    rows_tpm = aggregate(locus_map, per_sample, idx=1)
    write_matrix(samp_dir / "pri_expression_tpm.tsv", rows_tpm, samples, "sum")

    rows_cnt = aggregate(locus_map, per_sample, idx=0)
    write_matrix(samp_dir / "pri_expression_counts.tsv", rows_cnt, samples, "counts",
                 na_counts=(backend != "kallisto"))
    if backend != "kallisto":
        _status("stringtie backend: counts written as NA (no EM split); TPM only")

    st = {
        "backend": backend,
        "n_samples": len(per_sample),
        "n_loci": len(rows_tpm),
        "counts_available": backend == "kallisto",
        "note": ("counts unavailable: BAM-only quantification has no EM split; "
                 "supply FASTQ to get kallisto counts") if backend != "kallisto" else "",
    }
    (samp_dir / "quant_stats.json").write_text(json.dumps(st, indent=2), encoding="utf-8")
    rep.facts["quant"] = st
    _status(f"{st['n_samples']} samples x {st['n_loci']} loci ({backend})")
