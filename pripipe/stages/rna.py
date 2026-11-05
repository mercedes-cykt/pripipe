#!/usr/bin/env python3
"""RNA-Seq stage: FASTQ -> filtered -> spliced BAM (or accept existing BAMs).

Outputs ``<out>/rna/bams/<sample>.bam`` which both `assemble` and `quant` read,
plus a library strandedness call in ``strandness.txt``.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

from ..config import RunConfig
from ..utils import env_with_env_bin, missing_tool_message


def _status(msg: str) -> None:
    sys.stderr.write(f"[rna] {msg}\n")


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


def ensure_index(cfg: RunConfig, logs: Path) -> Path:
    """Return the hisat2 index prefix, building it only when unavoidable."""
    if cfg.hisat2_index:
        prefix = cfg.hisat2_index
        if not Path(str(prefix) + ".1.ht2").exists():
            raise SystemExit(f"rna: hisat2 index not found at {prefix}")
        return Path(prefix)
    prefix = cfg.out / "index" / "hisat2" / "genome"
    if Path(str(prefix) + ".1.ht2").exists():
        return prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    _status(f"building hisat2 index (this can take 20-40 min): {prefix}")
    _run(["hisat2-build", "-p", str(cfg.threads), str(cfg.genome), str(prefix)],
         logs / "hisat2_build.log")
    return prefix


def detect_strandness(bams: list[tuple[str, Path]], out: Path) -> str:
    """Read1 vs XS-tag concordance on the first available BAM."""
    try:
        import pysam
    except ImportError:
        return "unknown"
    for _, bam in bams:
        concordant = total = 0
        try:
            with pysam.AlignmentFile(str(bam), "rb") as af:
                for rec in af.fetch(until_eof=True):
                    if rec.is_unmapped or rec.is_read2:
                        continue
                    xs = rec.get_tag("XS") if rec.has_tag("XS") else None
                    if xs is None:
                        continue
                    total += 1
                    if (rec.is_reverse and xs == "-") or (not rec.is_reverse and xs == "+"):
                        concordant += 1
                    if total >= 200_000:
                        break
        except (OSError, ValueError):
            continue
        if total >= 1000:
            ratio = concordant / total
            call = "rf" if ratio < 0.35 else ("fr" if ratio > 0.65 else "unstranded")
            out.write_text(f"{call}\t{ratio:.3f}\t{total}\n", encoding="utf-8")
            return call
    return "unknown"


def stage_rna(cfg: RunConfig, rep) -> None:
    out_dir = cfg.stage_dir("rna")
    bam_dir = out_dir / "bams"
    logs = out_dir / "logs"
    bam_dir.mkdir(parents=True, exist_ok=True)

    from .quant import load_samples

    samples = load_samples(cfg)
    if not samples:
        raise SystemExit("rna: empty --rna sample sheet")

    # samples that already come with a BAM are linked, not aligned
    aligned: list[tuple[str, Path]] = []
    for s in samples:
        if s["bam"]:
            dst = bam_dir / f"{s['sample_id']}.bam"
            if not dst.exists():
                shutil.copy2(s["bam"], dst)
                if Path(s["bam"] + ".bai").exists():
                    shutil.copy2(s["bam"] + ".bai", Path(str(dst) + ".bai"))
            aligned.append((s["sample_id"], dst))

    to_align = [s for s in samples if s["r1"] and not s["bam"]]
    if to_align:
        index = ensure_index(cfg, logs)
        for s in to_align:
            sid = s["sample_id"]
            dst = bam_dir / f"{sid}.bam"
            if dst.exists():
                aligned.append((sid, dst))
                continue
            _status(f"fastp + hisat2: {sid}")
            fq1 = out_dir / "fastq" / f"{sid}_1.fq.gz"
            fq2 = out_dir / "fastq" / f"{sid}_2.fq.gz"
            fq1.parent.mkdir(parents=True, exist_ok=True)
            fastp = ["fastp", "-i", s["r1"], "-o", str(fq1),
                     "-w", str(min(8, cfg.threads)), "-j", str(out_dir / "fastq" / f"{sid}.json")]
            if s["r2"]:
                fastp += ["-I", s["r2"], "-O", str(fq2)]
            _run(fastp, logs / f"{sid}_fastp.log")

            hisat = ["hisat2", "-p", str(cfg.threads), "-x", str(index),
                     "--dta", "--novel-splicesite-outfile", str(out_dir / f"{sid}.splicesites")]
            hisat += ["-1", str(fq1)] + (["-2", str(fq2)] if s["r2"] else [])
            log = logs / f"{sid}_hisat2.log"
            log.parent.mkdir(parents=True, exist_ok=True)
            with log.open("a", encoding="utf-8") as fh:
                fh.write("\n$ " + " ".join(hisat) + " | samtools sort\n")
                fh.flush()
                h = subprocess.Popen(hisat, stdout=subprocess.PIPE, stderr=fh,
                                     env=env_with_env_bin())
                srt = subprocess.run(["samtools", "sort", "-@", str(cfg.threads),
                                      "-o", str(dst)], stdin=h.stdout, env=env_with_env_bin())
                h.wait()
            if h.returncode or srt.returncode:
                raise RuntimeError(f"alignment failed for {sid} (see {log})")
            _run(["samtools", "index", str(dst)], logs / f"{sid}_index.log")
            aligned.append((sid, dst))

    strand_file = out_dir / "strandness.txt"
    call = cfg.strandedness
    if call == "auto":
        call = detect_strandness(aligned, strand_file)
        _status(f"strandedness detected: {call}")
    else:
        strand_file.write_text(f"{call}\tna\tna\n", encoding="utf-8")

    rep.facts["rna"] = {"samples": len(aligned), "strandedness": call}
    _status(f"{len(aligned)} BAM ready in {bam_dir}")
