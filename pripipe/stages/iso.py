#!/usr/bin/env python3
"""Iso-Seq stage: raw reads -> collapsed transcript models.

The entry point is chosen by ``--iso-type``:

    subreads : ccs -> cat -> lima -> refine -> cluster2
    ccs      :         lima -> refine -> cluster2
    flnc     :                       cluster2
    hq       : (already clustered)   -------

then, for every mode: minimap2 -> isoseq3 collapse -> gffread -T.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

from ..config import RunConfig
from ..utils import env_with_env_bin, missing_tool_message


def _status(msg: str) -> None:
    sys.stderr.write(f"[iso] {msg}\n")


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


def _ccs(cfg: RunConfig, work: Path, logs: Path) -> list[Path]:
    outs: list[Path] = []
    for i, src in enumerate(cfg.iso, 1):
        dst = work / f"ccs_{i}.bam"
        if not dst.exists():
            _status(f"ccs: {src.name}")
            _run(["ccs", "--min-rq", str(cfg.iso_min_rq), "-j", str(cfg.threads),
                  str(src), str(dst)], logs / f"ccs_{i}.log")
        outs.append(dst)
    return outs


def _cat(bams: list[Path], dst: Path, logs: Path) -> Path:
    if dst.exists():
        return dst
    if len(bams) == 1:
        shutil.copy2(bams[0], dst)
        return dst
    _status(f"samtools cat ({len(bams)} cells)")
    _run(["samtools", "cat", "-o", str(dst), *map(str, bams)], logs / "samtools_cat.log")
    return dst


def _lima(cfg: RunConfig, bam: Path, dst: Path, logs: Path) -> Path:
    if dst.exists():
        return dst
    if not cfg.iso_primers:
        raise SystemExit("iso: --iso-primers is required for subreads/ccs input")

    # lima writes <prefix>.<primer5p>--<primer3p>.bam, never the bare prefix
    parts = sorted(p for p in dst.parent.glob(dst.stem + ".*.bam") if p.is_file())
    if not parts:
        _status("lima --isoseq: removing primers")
        _run(["lima", "--isoseq", str(bam), str(cfg.iso_primers), str(dst),
              "-j", str(cfg.threads)], logs / "lima.log")
        parts = sorted(p for p in dst.parent.glob(dst.stem + ".*.bam") if p.is_file())
    if not parts:
        raise RuntimeError(f"lima produced no BAM next to {dst}")
    if len(parts) == 1:
        shutil.move(str(parts[0]), str(dst))
    else:
        _status(f"samtools merge ({len(parts)} primer pairs)")
        _run(["samtools", "merge", "-f", str(dst), *map(str, parts)],
             logs / "lima_merge.log")
    return dst


def _refine(cfg: RunConfig, bam: Path, dst: Path, logs: Path) -> Path:
    """isoseq3 refine (>=4.x) takes: <demux.bam> <primers.fasta> <flnc.bam>."""
    if dst.exists():
        return dst
    if not cfg.iso_primers:
        raise SystemExit("iso: --iso-primers is required for isoseq3 refine")
    _status("isoseq3 refine -> FLNC")
    _run(["isoseq3", "refine", "--require-polya", "-j", str(cfg.threads),
          str(bam), str(cfg.iso_primers), str(dst)], logs / "refine.log")
    return dst


def _cluster(cfg: RunConfig, flnc: Path, dst: Path, logs: Path) -> Path:
    if dst.exists():
        return dst
    _status(f"isoseq3 cluster2 (singletons={'on' if cfg.iso_singletons else 'off'})")
    # --singletons keeps FLNCs that cannot be clustered: they carry most
    # low-abundance pri-miRNA isoforms, but they are also the lowest-evidence
    # ones, so the choice is exposed as --iso-singletons / --no-iso-singletons.
    extra = ["--singletons"] if cfg.iso_singletons else []
    _run(["isoseq3", "cluster2", *extra, "-j", str(cfg.threads),
          str(flnc), str(dst)], logs / "cluster.log")
    return dst


def bam_to_gtf(bam: Path, gtf: Path) -> int:
    """Fallback: build transcript models straight from alignment blocks.

    Used when `isoseq3 collapse` returns nothing (it expects cluster2-style
    reads with isoform tags; already-polished hq transcripts often collapse
    to zero).
    """
    try:
        import pysam
    except ImportError:
        raise SystemExit("iso: pysam is required for the BAM->GTF fallback")
    n = 0
    with pysam.AlignmentFile(str(bam), "rb") as af, gtf.open("w", encoding="utf-8") as out:
        for rec in af.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                continue
            blocks = rec.get_blocks()
            if not blocks:
                continue
            blocks.sort()
            tid = rec.query_name
            start, end = blocks[0][0], blocks[-1][1]
            chrom, strand = af.get_reference_name(rec.reference_id), "-" if rec.is_reverse else "+"
            out.write(f'{chrom}\tpripipe\tgene\t{start + 1}\t{end}\t.\t{strand}\t.\t'
                      f'gene_id "{tid}";\n')
            for i, (s, e) in enumerate(blocks, 1):
                out.write(f'{chrom}\tpripipe\texon\t{s + 1}\t{e}\t.\t{strand}\t.\t'
                          f'gene_id "{tid}"; transcript_id "{tid}"; exon_number "{i}";\n')
            out.write(f'{chrom}\tpripipe\ttranscript\t{start + 1}\t{end}\t.\t{strand}\t.\t'
                      f'gene_id "{tid}"; transcript_id "{tid}";\n')
            n += 1
    return n


def _cupcake_collapse(cfg: RunConfig, fasta: Path, aligned: Path, work: Path, logs: Path) -> Path:
    """Cupcake collapse - what the 2026-09 analysis used (-c 0.99 -i 0.95).

    Cupcake lives in the legacy python-2.7 `isoseq` env, so it is invoked
    through that interpreter and falls back to `isoseq3 collapse` when absent.
    """
    script = cfg.cupcake_dir / "collapse_isoforms_by_sam.py"
    py = cfg.cupcake_dir / "python"
    if not script.exists() or not py.exists():
        _status("Cupcake not found, falling back to isoseq3 collapse")
        out = work / "collapsed_isoseq3.gff"
        _run(["isoseq3", "collapse", str(aligned), str(out)], logs / "collapse.log")
        return out

    prefix = work / "collapsed"
    out = work / "collapsed.collapsed.gff"
    if not out.exists():
        sam = work / "aligned.sam"
        if not sam.exists():
            _status("samtools view: BAM -> SAM (Cupcake wants -s SAM)")
            _run(["samtools", "view", "-h", "-o", str(sam), str(aligned)],
                 logs / "samtools_sam.log")
        _status(f"Cupcake collapse (-c {cfg.iso_collapse_c} -i {cfg.iso_collapse_i})")
        _run([str(py), str(script), "--input", str(fasta), "-s", str(sam),
              "-o", str(prefix), "-c", str(cfg.iso_collapse_c),
              "-i", str(cfg.iso_collapse_i), "--dun-merge-5-shorter"],
             logs / "collapse.log")
    if not out.exists():
        raise RuntimeError(f"Cupcake collapse produced no {out.name}")
    return out


def _collapse(cfg: RunConfig, transcripts: Path, work: Path, logs: Path) -> Path:
    """Align transcripts to the genome, collapse redundancy, return GTF."""
    gtf = work / "transcripts.gtf"
    if gtf.exists():
        return gtf

    # minimap2 wants FASTA/FASTQ, cluster2 hands us BAM
    if transcripts.suffix.lower() == ".bam":
        fa = work / "transcripts.fasta"
        if not fa.exists():
            _status("samtools fasta: BAM -> FASTA")
            _run(["samtools", "fasta", "-@", str(cfg.threads), "-0", str(fa), str(transcripts)],
                 logs / "samtools_fasta.log")
        transcripts = fa

    aligned = work / "aligned.bam"
    if not aligned.exists():
        _status("minimap2 -ax splice:hq")
        log = logs / "minimap2.log"
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open("a", encoding="utf-8") as fh:
            fh.write("\n$ minimap2 -ax splice:hq | samtools sort\n")
            fh.flush()
            mm = subprocess.Popen(["minimap2", "-ax", "splice:hq", "-t", str(cfg.threads),
                                   str(cfg.genome), str(transcripts)],
                                  stdout=subprocess.PIPE, stderr=fh, env=env_with_env_bin())
            srt = subprocess.run(["samtools", "sort", "-@", str(cfg.threads), "-o", str(aligned)],
                                 stdin=mm.stdout, env=env_with_env_bin())
            mm.wait()
        if mm.returncode or srt.returncode:
            raise RuntimeError(f"minimap2 failed (see {log})")
        _run(["samtools", "index", str(aligned)], logs / "samtools_index.log")

    collapsed_gff = work / "collapsed.gff"
    if not collapsed_gff.exists():
        if cfg.collapse_strategy == "cupcake":
            collapsed_gff = _cupcake_collapse(cfg, transcripts, aligned, work, logs)
        else:
            _status("isoseq3 collapse")
            _run(["isoseq3", "collapse", str(aligned), str(collapsed_gff)],
                 logs / "collapse.log")

    # isoseq collapse emits GFF3; everything downstream wants GTF
    _status("gffread -T -> GTF")
    _run(["gffread", "-T", "-o", str(gtf), str(collapsed_gff)], logs / "gffread.log")

    n = sum(1 for line in gtf.open(encoding="utf-8", errors="replace") if "\ttranscript\t" in line)
    if n == 0:
        _status("collapse produced no transcript; falling back to BAM -> GTF")
        gtf.unlink(missing_ok=True)
        n = bam_to_gtf(aligned, gtf)
        _status(f"fallback built {n} transcript(s) from alignment blocks")
    return gtf


def stage_iso(cfg: RunConfig, rep) -> None:
    out_dir = cfg.stage_dir("iso")
    work = out_dir / "work"
    logs = out_dir / "logs"
    work.mkdir(parents=True, exist_ok=True)
    if not cfg.iso:
        _status("no Iso-Seq input, nothing to do")
        return

    if cfg.iso_type == "subreads":
        ccs_bams = _ccs(cfg, work, logs)
        merged = _cat(ccs_bams, work / "ccs_merged.bam", logs)
        trimmed = _lima(cfg, merged, work / "trimmed.bam", logs)
        flnc = _refine(cfg, trimmed, work / "flnc.bam", logs)
        transcripts = _cluster(cfg, flnc, work / "clustered.bam", logs)
    elif cfg.iso_type == "ccs":
        merged = _cat([Path(p) for p in cfg.iso], work / "ccs_merged.bam", logs)
        trimmed = _lima(cfg, merged, work / "trimmed.bam", logs)
        flnc = _refine(cfg, trimmed, work / "flnc.bam", logs)
        transcripts = _cluster(cfg, flnc, work / "clustered.bam", logs)
    elif cfg.iso_type == "flnc":
        flnc = _cat([Path(p) for p in cfg.iso], work / "flnc.bam", logs)
        transcripts = _cluster(cfg, flnc, work / "clustered.bam", logs)
    else:  # hq: FASTA of polished transcripts
        transcripts = Path(cfg.iso[0])

    gtf = _collapse(cfg, transcripts, work, logs)
    dst = out_dir / "transcripts.gtf"
    if gtf != dst:
        shutil.copy2(gtf, dst)

    n = sum(1 for line in dst.open(encoding="utf-8", errors="replace") if "\ttranscript\t" in line)
    rep.facts["iso"] = {"iso_type": cfg.iso_type, "transcripts": n}
    _status(f"{cfg.iso_type}: {n} transcripts")
