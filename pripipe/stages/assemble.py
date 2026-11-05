#!/usr/bin/env python3
"""Assemble RNA-Seq BAMs into transcript models with StringTie.

Per-sample assembly (recall matters: a pri that is tissue-specific must not be
diluted by merging BAMs first), then a single `--merge` across samples.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

from ..config import RunConfig
from ..utils import env_with_env_bin, missing_tool_message


def _status(msg: str) -> None:
    sys.stderr.write(f"[assemble] {msg}\n")


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


def collect_bams(cfg: RunConfig) -> list[tuple[str, Path]]:
    """Sample id -> BAM, preferring the rna stage output."""
    bam_dir = cfg.stage_dir("rna") / "bams"
    out: list[tuple[str, Path]] = []
    if bam_dir.exists():
        for bam in sorted(bam_dir.glob("*.bam")):
            out.append((bam.stem, bam))
    if out:
        return out
    from .quant import load_samples

    for s in load_samples(cfg):
        if s["bam"]:
            out.append((s["sample_id"], Path(s["bam"])))
    return out


def expressed_ids(gtf: Path, min_fpkm: float = 0.0) -> set[str]:
    """Transcript ids with FPKM > min in one StringTie per-sample GTF."""
    keep: set[str] = set()
    with gtf.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "\ttranscript\t" not in line:
                continue
            attrs = line.rstrip("\n").split("\t")[8]
            m = re.search(r'transcript_id\s+"([^"]+)"', attrs)
            fpkm = re.search(r'FPKM\s+"([0-9.eE+-]+)"', attrs)
            if not m or fpkm is None:
                continue
            if float(fpkm.group(1)) > min_fpkm:
                keep.add(m.group(1))
    return keep


def filter_expressed(gtfs: list[Path], merged: Path, dst: Path) -> dict:
    """Keep only merged transcripts that were actually assembled somewhere.

    `stringtie --merge -G ref.gtf` echoes every reference transcript back, which
    would make the RNA-Seq catalog indistinguishable from the annotation. The
    per-sample GTFs carry FPKM, so they decide what was really transcribed.
    """
    keep: set[str] = set()
    for gtf in gtfs:
        keep |= expressed_ids(gtf)

    total = kept = 0
    with merged.open(encoding="utf-8", errors="replace") as fin, \
            dst.open("w", encoding="utf-8") as fout:
        for line in fin:
            is_tx = "\ttranscript\t" in line
            if is_tx:
                total += 1
            m = re.search(r'transcript_id\s+"([^"]+)"', line)
            if m and m.group(1) in keep:
                fout.write(line)
                if is_tx:
                    kept += 1
    return {"transcripts_total": total, "kept": kept}


def stage_assemble(cfg: RunConfig, rep) -> None:
    out_dir = cfg.stage_dir("assemble")
    logs = out_dir / "logs"
    bams = collect_bams(cfg)
    if not bams:
        raise SystemExit("assemble: no BAM found - run the `rna` stage first")

    guide = ["-G", str(cfg.gtf)] if cfg.gtf else []
    per_sample: list[Path] = []
    for sid, bam in bams:
        gtf = out_dir / f"{sid}.gtf"
        if not gtf.exists():
            _status(f"stringtie: {sid}")
            _run(["stringtie", *guide, "-p", str(cfg.threads), "-o", str(gtf), str(bam)],
                 logs / f"{sid}.log")
        per_sample.append(gtf)

    merged = out_dir / "merged.gtf"
    if not merged.exists():
        lst = out_dir / "merge_list.txt"
        lst.write_text("\n".join(str(p) for p in per_sample) + "\n", encoding="utf-8")
        _status(f"stringtie --merge ({len(per_sample)} samples)")
        _run(["stringtie", "--merge", "-o", str(merged), str(lst)],
             logs / "merge.log")

    shutil.copy2(merged, out_dir / "transcripts.gtf")
    n = sum(1 for line in merged.open(encoding="utf-8", errors="replace")
            if "\ttranscript\t" in line)
    rep.facts["assemble"] = {"samples": len(per_sample), "transcripts": n}
    _status(f"{len(per_sample)} samples -> {n} merged transcripts")
