#!/usr/bin/env python3
"""Catalog stage: turn assembled transcripts into the pripipe transcript catalog.

Two ways in:

* ``--existing-catalog DIR`` - reuse a finished SQANTI3 output directory
  (``*_corrected.gtf`` + ``*_classification.txt``). Nothing is recomputed.
* otherwise - run SQANTI3 (mandatory dependency) on the transcripts produced by
  the ``iso`` and/or ``assemble`` stages, then build the union catalog.

Outputs (``<out>/catalog/``):
    transcripts.gtf        union catalog (all sources)
    transcripts.iso.gtf    Iso-Seq-only catalog   (when both sources exist)
    transcripts.rna.gtf    RNA-Seq-only catalog   (when both sources exist)
    meta.json              provenance: sources, inputs, versions, counts
    <name>_classification.txt (copied from SQANTI3)
    fl_counts.tsv          Iso-Seq FL counts (when available)
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

from ..config import RunConfig
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
def run_sqanti3(cfg: RunConfig, isoforms: Path, outdir: Path) -> Path:
    qc = cfg.sqanti3_dir / "sqanti3_qc.py"
    if not qc.exists():
        raise SystemExit(f"catalog: SQANTI3 not found at {cfg.sqanti3_dir}")
    if cfg.gtf is None:
        raise SystemExit("catalog: SQANTI3 needs --gtf (reference annotation)")
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
           "-o", "pripipe", "-d", str(outdir)]
    _status(f"SQANTI3 on {isoforms.name}")
    _run(cmd, cfg.stage_dir("catalog") / "logs" / "sqanti3.log")
    gtf = pick_sqanti3_gtf(outdir)
    if gtf is None:
        raise SystemExit(f"catalog: SQANTI3 produced no GTF in {outdir}")
    return gtf


# --------------------------------------------------------------------------
# stage entry
# --------------------------------------------------------------------------
def stage_catalog(cfg: RunConfig, rep) -> None:
    cat_dir = cfg.stage_dir("catalog")
    cat_dir.mkdir(parents=True, exist_ok=True)

    if cfg.existing_catalog:
        src_dir = cfg.existing_catalog
        gtf = pick_sqanti3_gtf(src_dir)
        if gtf is None:
            raise SystemExit(f"catalog: no *_corrected.gtf / *.gtf in {src_dir}")
        _status(f"reusing existing SQANTI3 output: {gtf.name}")
        if cfg.catalog_sources:
            sources = [{"iso": "Iso-Seq", "rna": "RNA-Seq"}[s] for s in cfg.catalog_sources.split(",")]
        else:
            sources = {"iso": ["Iso-Seq"], "rna": ["RNA-Seq"],
                       "both": ["Iso-Seq", "RNA-Seq"]}[cfg.mode]
        normalize_gtf(gtf, cat_dir / "transcripts.gtf")
        cls = copy_classification(src_dir, cat_dir)
        n_fl = read_abundance_fl(src_dir, cat_dir / "fl_counts.tsv") or read_abundance_fl(src_dir.parent, cat_dir / "fl_counts.tsv")
        sq_version = _sqanti3_version(cfg)
    else:
        iso_gtf = _find(cfg.stage_dir("iso"), ["transcripts.gtf", "*.corrected.gtf"])
        rna_gtf = _find(cfg.stage_dir("assemble"), ["merged.gtf", "transcripts.gtf"])
        if iso_gtf is None and rna_gtf is None:
            raise SystemExit("catalog: nothing to work with - run `iso` and/or `assemble` first")

        sources = [s for s, p in (("Iso-Seq", iso_gtf), ("RNA-Seq", rna_gtf)) if p]
        parts: list[tuple[str, Path]] = []
        if iso_gtf:
            corrected = run_sqanti3(cfg, iso_gtf, cat_dir / "sqanti3_iso")
            normalize_gtf(corrected, cat_dir / "transcripts.iso.gtf")
            cls_src = corrected.parent / "pripipe_classification.txt"
            if cls_src.exists():
                shutil.copy2(cls_src, cat_dir / cls_src.name)
            parts.append(("Iso-Seq", cat_dir / "transcripts.iso.gtf"))
            n_fl = read_abundance_fl(iso_gtf.parent, cat_dir / "fl_counts.tsv")
        else:
            n_fl = 0
        if rna_gtf:
            normalize_gtf(rna_gtf, cat_dir / "transcripts.rna.gtf")
            parts.append(("RNA-Seq", cat_dir / "transcripts.rna.gtf"))

        if len(parts) == 1:
            shutil.copy2(parts[0][1], cat_dir / "transcripts.gtf")
        else:
            merged = merge_catalogs([(src, parse_gtf(p, source=src)) for src, p in parts])
            _write_gtf(merged, cat_dir / "transcripts.gtf")
        cls = copy_classification(cat_dir, cat_dir)
        sq_version = _sqanti3_version(cfg)

    txs = parse_gtf(cat_dir / "transcripts.gtf")
    meta = {
        "sources": sources,
        "n_transcripts": len(txs),
        "sqanti3": {"dir": str(cfg.sqanti3_dir), "version": sq_version},
        "inputs": {"existing_catalog": str(cfg.existing_catalog) if cfg.existing_catalog else None},
        "classification": cls,
        "fl_counts": int(n_fl),
    }
    (cat_dir / "meta.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False),
                                       encoding="utf-8")
    rep.facts["catalog"] = meta
    _status(f"{len(txs)} transcripts, sources={sources}, fl={n_fl}")


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


def _write_gtf(txs: dict, dst: Path) -> None:
    with dst.open("w", encoding="utf-8") as fh:
        for tx in txs.values():
            # no gene feature: gffread would emit a second FASTA record under the
            # same name (gene_id == transcript_id), which breaks kallisto indexing
            for i, (s, e) in enumerate(tx.exons, 1):
                fh.write(f'{tx.chrom}\tpripipe\texon\t{s + 1}\t{e}\t.\t{tx.strand}\t.\t'
                         f'gene_id "{tx.tid}"; transcript_id "{tx.tid}"; exon_number "{i}";\n')
            fh.write(f'{tx.chrom}\tpripipe\ttranscript\t{tx.start + 1}\t{tx.end}\t.\t'
                     f'{tx.strand}\t.\tgene_id "{tx.tid}"; transcript_id "{tx.tid}";\n')
