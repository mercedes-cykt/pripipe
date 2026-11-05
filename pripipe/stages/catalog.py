#!/usr/bin/env python3
"""Catalog stage: turn assembled transcripts into the pripipe transcript catalog.

Three ways in:

* ``--existing-catalog DIR`` - reuse a finished SQANTI3 output directory
  (``*_corrected.gtf`` + ``*_classification.txt``). Nothing is recomputed.
* Iso-Seq and/or RNA-Seq transcripts - SQANTI3 is run on each source when a
  reference annotation is available, then the sources are merged.
* nothing - error.

After merging, the catalog is filtered (length, span, SQANTI3 class) and made
non-redundant (structural clustering by default), so that downstream evidence
counts are not inflated by near-identical models.

Outputs (``<out>/catalog/``):
    transcripts.gtf        filtered, non-redundant catalog
    transcripts.raw.gtf    catalog before filtering (kept for auditing)
    transcripts.iso.gtf    Iso-Seq-only models   (when both sources exist)
    transcripts.rna.gtf    RNA-Seq-only models   (when both sources exist)
    meta.json              provenance + filter/dedup statistics
    dropped.tsv            every removed transcript with its reason
    flags.tsv              transcripts kept but flagged (e.g. long mono-exon)
    <name>_classification.txt (copied from SQANTI3)
    fl_counts.tsv          Iso-Seq FL counts (when available)
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from ..config import RunConfig
from ..models import Transcript
from ..utils import env_with_env_bin, missing_tool_message
from .anchor import merge_catalogs, parse_gtf

ATTR_TX = re.compile(r'transcript_id\s+"?([^";]+)"?')


def _status(msg: str) -> None:
    sys.stderr.write(f"[catalog] {msg}\n")


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
# GTF normalisation
# --------------------------------------------------------------------------
def normalize_gtf(src: Path, dst: Path) -> dict:
    """Drop features without transcript_id and repair broken gene_id values."""
    kept = dropped = fixed = 0
    dst.parent.mkdir(parents=True, exist_ok=True)
    with src.open(encoding="utf-8-sig", errors="replace") as fin, \
            dst.open("w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("#") or not line.strip():
                fout.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                dropped += 1
                continue
            attrs = f[8]
            m = ATTR_TX.search(attrs)
            if not m:
                dropped += 1
                continue
            tid = m.group(1)
            gene = re.search(r'gene_id\s+"?([^";]+)"?', attrs)
            if gene is None or gene.group(1) in ("PB", "", "NA"):
                attrs = re.sub(r'gene_id\s+"?[^";]*"?;?\s*', "", attrs).strip().rstrip(";")
                attrs = f'gene_id "{tid}"; ' + attrs
                fixed += 1
            f[8] = attrs
            fout.write("\t".join(f) + "\n")
            kept += 1
    return {"kept": kept, "dropped": dropped, "gene_id_fixed": fixed}


def pick_sqanti3_gtf(directory: Path) -> Path | None:
    for pattern in ("*_corrected.gtf", "*.gtf"):
        hits = sorted(directory.glob(pattern))
        if hits:
            return hits[0]
    return None


def copy_classification(src_dir: Path, dst_dir: Path) -> str | None:
    hits = sorted(src_dir.glob("*classification.txt"))
    if not hits:
        return None
    src = hits[0]
    dst = dst_dir / src.name
    if src.resolve() != dst.resolve():
        shutil.copy2(src, dst)
    return src.name


def read_abundance_fl(src_dir: Path, dst: Path) -> int:
    """Iso-Seq copy number per transcript (`*.abundance.txt`: tid <tab> count)."""
    hits = sorted(list(src_dir.glob("*.abundance.txt")) + list(src_dir.glob("*_fl_count*")))
    if not hits:
        return 0
    rows: dict[str, float] = {}
    for line in hits[0].read_text(encoding="utf-8", errors="replace").splitlines():
        f = line.split("\t")
        if len(f) < 2:
            continue
        try:
            rows[f[0]] = float(f[-1])
        except ValueError:
            continue
    if not rows:
        return 0
    with dst.open("w", encoding="utf-8") as fh:
        fh.write("transcript_id\tn_fl\n")
        for tid, n in sorted(rows.items()):
            fh.write(f"{tid}\t{int(n)}\n")
    return len(rows)


# --------------------------------------------------------------------------
# SQANTI3
# --------------------------------------------------------------------------
def run_sqanti3(cfg: RunConfig, isoforms: Path, outdir: Path, tag: str) -> Path:
    qc = cfg.sqanti3_dir / "sqanti3_qc.py"
    if not qc.exists():
        raise SystemExit(f"catalog: SQANTI3 not found at {cfg.sqanti3_dir}")
    outdir.mkdir(parents=True, exist_ok=True)
    cached = pick_sqanti3_gtf(outdir)
    if cached is not None:
        _status(f"reusing SQANTI3 output: {cached.name}")
        return cached
    cmd = [sys.executable, str(qc),
           "--isoforms", str(isoforms),
           "--refGTF", str(cfg.gtf),
           "--refFasta", str(cfg.genome),
           "--skipORF", "-t", str(cfg.threads),
           "-o", tag, "-d", str(outdir)]
    _status(f"SQANTI3 on {isoforms.name} -> {tag}")
    _run(cmd, cfg.stage_dir("catalog") / "logs" / f"sqanti3_{tag}.log")
    gtf = pick_sqanti3_gtf(outdir)
    if gtf is None:
        raise SystemExit(f"catalog: SQANTI3 produced no GTF in {outdir}")
    return gtf


def read_classification(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open(encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if "isoform" not in header or "structural_category" not in header:
            return out
        i_id, i_cat = header.index("isoform"), header.index("structural_category")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > i_cat:
                out[f[i_id]] = f[i_cat]
    return out


def collect_classes(cfg: RunConfig) -> dict[str, str]:
    classes: dict[str, str] = {}
    for path in sorted(cfg.stage_dir("catalog").rglob("*classification.txt")):
        classes.update(read_classification(path))
    return classes


# --------------------------------------------------------------------------
# filtering
# --------------------------------------------------------------------------
def apply_filters(txs: dict[str, Transcript], classes: dict[str, str],
                  cfg: RunConfig) -> tuple[dict[str, Transcript], list[tuple], list[tuple], dict]:
    """Remove implausible models; flag (but keep) suspicious ones.

    Returns (kept, dropped, flagged, stats). Nothing is removed silently: every
    removal lands in dropped.tsv with a reason.
    """
    drop_classes = {c.strip() for c in cfg.drop_classes.split(",") if c.strip()}
    kept: dict[str, Transcript] = {}
    dropped: list[tuple[str, str, int]] = []
    flagged: list[tuple[str, str, int]] = []

    for tid, tx in txs.items():
        length = tx.end - tx.start
        if length < cfg.min_transcript_length:
            dropped.append((tid, "shorter_than_min", length))
            continue
        if length > cfg.max_transcript_length:
            dropped.append((tid, "longer_than_max", length))
            continue
        category = classes.get(tid)
        if category and category in drop_classes:
            dropped.append((tid, f"class:{category}", length))
            continue
        if tx.n_exons == 1 and length > cfg.max_mono_exon_span:
            flagged.append((tid, "mono_exon_span", length))
        kept[tid] = tx

    stats = {
        "n_in": len(txs),
        "n_out": len(kept),
        "n_dropped": len(dropped),
        "n_flagged": len(flagged),
        "min_transcript_length": cfg.min_transcript_length,
        "max_transcript_length": cfg.max_transcript_length,
        "max_mono_exon_span": cfg.max_mono_exon_span,
        "drop_classes": sorted(drop_classes),
        "dropped_by_reason": dict(sorted(
            ((reason, sum(1 for _, r, _ in dropped if r == reason))
             for reason in {r for _, r, _ in dropped}), key=lambda kv: -kv[1])),
    }
    return kept, dropped, flagged, stats


# --------------------------------------------------------------------------
# redundancy removal
# --------------------------------------------------------------------------
def dedup_structural(txs: dict[str, Transcript], tss_tol: int,
                     pas_tol: int) -> tuple[dict[str, Transcript], list[tuple]]:
    """Collapse transcripts sharing a junction chain and a close TSS/PAS.

    Keeps the longest model of each cluster and records what it absorbed.
    """
    groups: dict[tuple, list[Transcript]] = defaultdict(list)
    for tx in txs.values():
        inner = tuple(tx.exons[1:-1]) if tx.n_exons > 2 else ()
        groups[(tx.chrom, tx.strand, tx.n_exons, inner)].append(tx)

    kept: dict[str, Transcript] = {}
    collapsed: list[tuple[str, str]] = []

    def flush(cluster: list[Transcript]) -> None:
        rep = max(cluster, key=lambda t: t.end - t.start)
        kept[rep.tid] = rep
        for member in cluster:
            if member.tid != rep.tid:
                collapsed.append((member.tid, rep.tid))

    for members in groups.values():
        members.sort(key=lambda t: t.start)
        open_clusters: list[list[Transcript]] = []
        for tx in members:
            # clusters that fell out of the TSS window are finished: flush them
            still_open: list[list[Transcript]] = []
            for cluster in open_clusters:
                if tx.start - cluster[0].start <= tss_tol:
                    still_open.append(cluster)
                else:
                    flush(cluster)
            open_clusters = still_open
            for cluster in open_clusters:
                if abs(tx.end - cluster[0].end) <= pas_tol:
                    cluster.append(tx)
                    break
            else:
                open_clusters.append([tx])
        for cluster in open_clusters:
            flush(cluster)
    return kept, collapsed


def dedup_cdhit(cfg: RunConfig, txs: dict[str, Transcript], work: Path,
                logs: Path) -> tuple[dict[str, Transcript], list[tuple]] | None:
    """Optional CD-HIT-EST pass; returns None when cd-hit-est is unavailable."""
    if shutil.which("cd-hit-est", path=str(Path(sys.executable).parent)) is None \
            and shutil.which("cd-hit-est") is None:
        _status("cd-hit-est not found, falling back to structural dedup")
        return None

    work.mkdir(parents=True, exist_ok=True)
    gtf = work / "dedup_in.gtf"
    _write_gtf(txs, gtf)
    fasta = work / "dedup_in.fa"
    _run(["gffread", "-w", str(fasta), "-g", str(cfg.genome), str(gtf)],
         logs / "dedup_gffread.log")

    out_fa = work / "dedup_out.fa"
    _run(["cd-hit-est", "-i", str(fasta), "-o", str(out_fa),
          "-c", str(cfg.cdhit_identity), "-n", "9", "-d", "0", "-M", "0",
          "-T", str(cfg.threads)], logs / "cdhit.log")

    members: dict[str, list[str]] = defaultdict(list)
    current: list[str] = []
    for line in Path(str(out_fa) + ".clstr").read_text().splitlines():
        if line.startswith(">Cluster"):
            if current:
                members[current[0]] = current
            current = []
        else:
            m = re.search(r">(\S+)\.\.\.", line)
            if m:
                current.append(m.group(1))
    if current:
        members[current[0]] = current

    kept = {tid: txs[tid] for tid in members if tid in txs}
    collapsed = [(member, rep) for rep, group in members.items()
                 for member in group if member != rep and member in txs and rep in txs]
    return kept, collapsed


# --------------------------------------------------------------------------
# stage entry
# --------------------------------------------------------------------------
def stage_catalog(cfg: RunConfig, rep) -> None:
    cat_dir = cfg.stage_dir("catalog")
    cat_dir.mkdir(parents=True, exist_ok=True)
    logs = cat_dir / "logs"

    if cfg.existing_catalog:
        src_dir = cfg.existing_catalog
        gtf = pick_sqanti3_gtf(src_dir)
        if gtf is None:
            raise SystemExit(f"catalog: no *_corrected.gtf / *.gtf in {src_dir}")
        _status(f"reusing existing SQANTI3 output: {gtf.name}")
        if cfg.catalog_sources:
            sources = [{"iso": "Iso-Seq", "rna": "RNA-Seq"}[s]
                       for s in cfg.catalog_sources.split(",")]
        else:
            sources = {"iso": ["Iso-Seq"], "rna": ["RNA-Seq"],
                       "both": ["Iso-Seq", "RNA-Seq"]}[cfg.mode]
        normalize_gtf(gtf, cat_dir / "transcripts.raw.gtf")
        cls = copy_classification(src_dir, cat_dir)
        n_fl = (read_abundance_fl(src_dir, cat_dir / "fl_counts.tsv")
                or read_abundance_fl(src_dir.parent, cat_dir / "fl_counts.tsv"))
    else:
        iso_gtf = _find(cfg.stage_dir("iso"), ["transcripts.gtf", "*.corrected.gtf"])
        rna_gtf = _find(cfg.stage_dir("assemble"), ["merged.gtf", "transcripts.gtf"])
        if iso_gtf is None and rna_gtf is None:
            raise SystemExit("catalog: nothing to work with - run `iso` and/or `assemble` first")

        sources = [s for s, p in (("Iso-Seq", iso_gtf), ("RNA-Seq", rna_gtf)) if p]
        parts: list[tuple[str, Path]] = []
        n_fl = 0
        for tag, source, path in (("iso", "Iso-Seq", iso_gtf), ("rna", "RNA-Seq", rna_gtf)):
            if path is None:
                continue
            dst = cat_dir / f"transcripts.{tag}.gtf"
            if cfg.gtf:
                corrected = run_sqanti3(cfg, path, cat_dir / f"sqanti3_{tag}", tag)
                normalize_gtf(corrected, dst)
                cls_src = corrected.parent / f"{tag}_classification.txt"
                if cls_src.exists():
                    shutil.copy2(cls_src, cat_dir / cls_src.name)
            else:
                _status(f"no --gtf: SQANTI3 skipped for {source} (no class_code)")
                normalize_gtf(path, dst)
            parts.append((source, dst))
            if tag == "iso":
                n_fl = read_abundance_fl(path.parent, cat_dir / "fl_counts.tsv")

        if len(parts) == 1:
            shutil.copy2(parts[0][1], cat_dir / "transcripts.raw.gtf")
        else:
            merged = merge_catalogs([(src, parse_gtf(p, source=src)) for src, p in parts])
            _write_gtf(merged, cat_dir / "transcripts.raw.gtf")
        cls = copy_classification(cat_dir, cat_dir)

    # ---- filter + deduplicate -------------------------------------------
    raw = parse_gtf(cat_dir / "transcripts.raw.gtf", source="+".join(sources))
    classes = collect_classes(cfg)
    kept, dropped, flagged, fstats = apply_filters(raw, classes, cfg)

    collapsed: list[tuple[str, str]] = []
    if cfg.dedup == "structural":
        kept, collapsed = dedup_structural(kept, cfg.dedup_tss_tolerance, cfg.dedup_pas_tolerance)
    elif cfg.dedup == "cdhit":
        result = dedup_cdhit(cfg, kept, cat_dir / "dedup_work", logs)
        if result is None:
            kept, collapsed = dedup_structural(kept, cfg.dedup_tss_tolerance,
                                               cfg.dedup_pas_tolerance)
        else:
            kept, collapsed = result

    expected = fstats["n_out"] - len(kept)
    if expected != len(collapsed):
        _status(f"WARNING: accounting mismatch, {expected} removed but "
                f"{len(collapsed)} collapses recorded")

    _write_gtf(kept, cat_dir / "transcripts.gtf")
    _write_rows(cat_dir / "dropped.tsv", ["transcript_id", "reason", "length"], dropped)
    _write_rows(cat_dir / "flags.tsv", ["transcript_id", "flag", "length"], flagged)
    _write_rows(cat_dir / "collapsed.tsv", ["transcript_id", "representative"], collapsed)

    dstats = {
        "mode": cfg.dedup,
        "tss_tolerance": cfg.dedup_tss_tolerance,
        "pas_tolerance": cfg.dedup_pas_tolerance,
        "n_before": fstats["n_out"],
        "n_after": len(kept),
        "n_collapsed": len(collapsed),
        "redundancy_before": round(fstats["n_out"] / len(kept), 3) if kept else 1.0,
    }
    (cat_dir / "filter_stats.json").write_text(
        json.dumps({"filter": fstats, "dedup": dstats}, indent=2), encoding="utf-8")

    txs = kept
    meta = {
        "sources": sources,
        "n_transcripts": len(txs),
        "sqanti3": {"dir": str(cfg.sqanti3_dir), "version": _sqanti3_version(cfg)},
        "inputs": {"existing_catalog": str(cfg.existing_catalog) if cfg.existing_catalog else None},
        "classification": cls,
        "fl_counts": int(n_fl),
        "filter": fstats,
        "dedup": dstats,
    }
    (cat_dir / "meta.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False),
                                       encoding="utf-8")
    rep.facts["catalog"] = meta
    _status(f"{len(raw)} -> {len(txs)} transcripts "
            f"(dropped {fstats['n_dropped']}, collapsed {len(collapsed)}, "
            f"flagged {fstats['n_flagged']}); sources={sources}")


def _find(directory: Path, patterns: list[str]) -> Path | None:
    if not directory.exists():
        return None
    for pattern in patterns:
        hits = sorted(directory.glob(pattern))
        if hits:
            return hits[0]
    return None


def _sqanti3_version(cfg: RunConfig) -> str:
    try:
        proc = subprocess.run(["git", "-C", str(cfg.sqanti3_dir), "describe", "--tags"],
                              capture_output=True, text=True, env=env_with_env_bin())
        return proc.stdout.strip() or "unknown"
    except (OSError, FileNotFoundError):
        return "unknown"


def _write_gtf(txs: dict[str, Transcript], dst: Path) -> None:
    with dst.open("w", encoding="utf-8") as fh:
        for tx in txs.values():
            # no gene feature: gffread would emit a second FASTA record under the
            # same name (gene_id == transcript_id), which breaks kallisto indexing
            for i, (s, e) in enumerate(tx.exons, 1):
                fh.write(f'{tx.chrom}\tpripipe\texon\t{s + 1}\t{e}\t.\t{tx.strand}\t.\t'
                         f'gene_id "{tx.tid}"; transcript_id "{tx.tid}"; exon_number "{i}";\n')
            fh.write(f'{tx.chrom}\tpripipe\ttranscript\t{tx.start + 1}\t{tx.end}\t.\t'
                     f'{tx.strand}\t.\tgene_id "{tx.tid}"; transcript_id "{tx.tid}";\n')


def _write_rows(path: Path, header: list[str], rows: list[tuple]) -> None:
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(v) for v in row) + "\n")
