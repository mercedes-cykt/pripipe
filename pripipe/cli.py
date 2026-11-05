#!/usr/bin/env python3
"""pripipe command line entry point.

Two commands:

    pripipe validate ...     check every input and stop
    pripipe run ...          validate, then run the pipeline
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .config import STAGES, RunConfig
from .validate import validate

EPILOG = """\
examples:
  # 1) Iso-Seq only, starting from raw subreads (primers required)
  pripipe run --genome Sta_v4.fasta --mirna-bed mirna_pre.bed --gtf Sta.gtf \\
      --iso cell1.subreads.bam cell2.subreads.bam --iso-type subreads \\
      --iso-primers primer.fasta --out results/ --threads 32

  # 2) Iso-Seq FLNC + RNA-Seq FASTQ, default quantification backend
  pripipe run --genome Sta_v4.fasta --mirna-bed mirna_pre.bed --gtf Sta.gtf \\
      --iso flnc.bam --iso-type flnc --rna samples.tsv --out results/

  # 3) reuse a finished SQANTI3 run, only recompute anchoring + report
  pripipe run --genome Sta_v4.fasta --mirna-bed mirna_pre.bed --gtf Sta.gtf \\
      --iso flnc.bam --iso-type flnc --existing-catalog sqanti3_out \\
      --catalog-sources iso --from catalog --to report --out results/

  # 4) check inputs only, do not touch anything
  pripipe validate --genome Sta_v4.fasta --mirna-bed mirna_pre.bed --json qc.json

stages (run in this order):
  iso       Iso-Seq reads -> collapsed transcript models
  rna       FASTQ -> filtered, spliced BAM (or accept existing BAM)
  assemble  StringTie per sample + merge
  catalog   SQANTI3 (mandatory) -> normalised transcript catalog + provenance
  anchor    pre-miRNA x transcript anchoring + processing intermediates
  quant     kallisto (FASTQ) or stringtie -e -B (BAM only)
  report    pri_mirna.tsv, pri_pairs.tsv, expression matrices, manifest.json

docs: README.md, SPEC.md
"""


def _common(p: argparse.ArgumentParser) -> None:
    req = p.add_argument_group("required inputs")
    ref = p.add_argument_group("reference (optional but recommended)")
    iso = p.add_argument_group("Iso-Seq input")
    rna = p.add_argument_group("RNA-Seq input")
    anc = p.add_argument_group("anchoring criteria")
    qnt = p.add_argument_group("quantification")
    cat = p.add_argument_group("catalog filtering")
    cat.add_argument("--min-transcript-length", type=int, default=200, metavar="N",
                     help="drop transcripts shorter than N bp (default: %(default)s)")
    cat.add_argument("--max-transcript-length", type=int, default=100_000, metavar="N",
                     help="drop transcripts spanning more than N bp (default: %(default)s)")
    cat.add_argument("--max-mono-exon-span", type=int, default=10_000, metavar="N",
                     help="flag single-exon transcripts spanning more than N bp "
                          "(flagged, never dropped; default: %(default)s)")
    cat.add_argument("--drop-classes", default="fusion", metavar="LIST",
                     help="comma-separated SQANTI3 categories to drop (default: %(default)s)")
    cat.add_argument("--dedup", choices=["structural", "cdhit", "none"], default="structural",
                     help="redundancy removal: structural clustering on exon chain + TSS/PAS "
                          "tolerance, cd-hit-est, or none (default: %(default)s)")
    cat.add_argument("--dedup-tss-tolerance", type=int, default=50, metavar="N",
                     help="TSS tolerance for structural dedup (default: %(default)s)")
    cat.add_argument("--dedup-pas-tolerance", type=int, default=50, metavar="N",
                     help="PAS tolerance for structural dedup (default: %(default)s)")
    cat.add_argument("--cdhit-identity", type=float, default=0.95, metavar="F",
                     help="cd-hit-est identity threshold (default: %(default)s)")

    dep = p.add_argument_group("external dependencies")
    run = p.add_argument_group("execution")

    req.add_argument("--genome", required=True, type=Path, metavar="FASTA",
                     help="reference genome; .fai is built automatically when missing")
    req.add_argument("--mirna-bed", required=True, type=Path, metavar="BED",
                     help="pre-miRNA loci, BED6 (name must be unique, strand must be +/-)")

    ref.add_argument("--gtf", type=Path, metavar="GTF",
                     help="reference annotation (GTF/GFF3). Without it: no class_code, "
                          "no genic/intergenic label, weaker strand inheritance")

    iso.add_argument("--iso", nargs="+", type=Path, default=[], metavar="FILE",
                     help="Iso-Seq input(s); several files are merged (multi-cell libraries)")
    iso.add_argument("--iso-type", choices=["subreads", "ccs", "flnc", "hq"],
                     default="subreads",
                     help="where to start: subreads (ccs+lima+refine+cluster2), ccs "
                          "(lima+refine+cluster2), flnc (cluster2) or hq (already clustered "
                          "FASTA/BAM). default: %(default)s")
    iso.add_argument("--iso-primers", type=Path, metavar="FASTA",
                     help="primer sequences (lima/refine); required for --iso-type subreads|ccs")
    iso.add_argument("--iso-barcode-set", type=Path, metavar="XML",
                     help="barcode set XML for multiplexed libraries (refine)")
    iso.add_argument("--iso-min-rq", type=float, default=0.9, metavar="F",
                     help="ccs --min-rq, i.e. minimum predicted accuracy (default: %(default)s)")
    iso.add_argument("--collapse-strategy", choices=["cupcake", "isoseq3"], default="cupcake",
                     help="redundancy removal: cupcake (-c/-i, the 2026-09 analysis) or "
                          "isoseq3 collapse. default: %(default)s")
    iso.add_argument("--iso-singletons", dest="iso_singletons", action="store_true",
                     default=True,
                     help="keep clustered-read singletons: higher recall, more "
                          "single-transcript loci (default: on)")
    iso.add_argument("--no-iso-singletons", dest="iso_singletons", action="store_false",
                     help="drop singletons: matches the stricter 2026-09 baseline")

    rna.add_argument("--rna", type=Path, metavar="TSV",
                     help="sample sheet (header required): sample_id, r1, r2, [bam], [group]; "
                          "see README for the exact format")
    rna.add_argument("--strandedness", choices=["auto", "rf", "fr", "unstranded"],
                     default="auto", help="library strandedness (default: %(default)s, detected)")
    rna.add_argument("--hisat2-index", type=Path, metavar="PREFIX",
                     help="reuse an existing hisat2 index instead of building one")
    rna.add_argument("--fastp-clean", dest="fastp_clean", action="store_true", default=True,
                     help="clean raw FASTQ with fastp before alignment (default: on)")
    rna.add_argument("--no-fastp-clean", dest="fastp_clean", action="store_false",
                     help="skip fastp and align the input FASTQ as given (already-cleaned data)")

    anc.add_argument("--extend-hairpin", type=int, default=0, metavar="N",
                     help="grow each hairpin by N bp before anchoring, as a boundary "
                          "tolerance (default: %(default)s)")
    anc.add_argument("--min-flank", type=int, default=50, metavar="N",
                     help="flag transcripts with less than N bp flanking the hairpin "
                          "(flagged, never dropped; default: %(default)s)")
    anc.add_argument("--allow-partial", action="store_true",
                     help="also report transcripts that cover the hairpin only partially")
    anc.add_argument("--multi-hairpin", choices=["per-hairpin", "merge"], default="per-hairpin",
                     help="how to report transcripts carrying several hairpins (default: %(default)s)")

    qnt.add_argument("--quant", dest="quant", action="store_true", default=True,
                     help="quantify pri-miRNA loci (default: on; skipped when no RNA-Seq)")
    qnt.add_argument("--no-quant", dest="quant", action="store_false",
                     help="annotation only, skip quantification")
    qnt.add_argument("--quant-backend", choices=["auto", "kallisto", "stringtie"],
                     default="auto",
                     help="auto: FASTQ -> kallisto, BAM only -> stringtie -e -B "
                          "(default: %(default)s)")

    dep.add_argument("--sqanti3-dir", type=Path,
                     default=Path("/nfs/users/zhaojw/software/SQANTI3"),
                     help="SQANTI3 source directory (mandatory dependency)")
    dep.add_argument("--cupcake-dir", type=Path,
                     default=Path("/nfs/users/zhaojw/miniconda3/envs/isoseq/bin"),
                     help="directory holding collapse_isoforms_by_sam.py (+ its python)")
    dep.add_argument("--existing-catalog", type=Path, metavar="DIR",
                     help="reuse a finished SQANTI3 output directory instead of running it")
    dep.add_argument("--catalog-sources", choices=["iso", "rna", "both"],
                     help="provenance to record when using --existing-catalog")

    run.add_argument("--out", type=Path, default=Path("pripipe_out"),
                     help="output directory (default: %(default)s)")
    run.add_argument("--threads", type=int, default=8, help="threads (default: %(default)s)")
    run.add_argument("--tmp-dir", type=Path, default=Path("./tmp"),
                     help="scratch directory (default: %(default)s; never use /tmp)")


def _config(args: argparse.Namespace) -> RunConfig:
    return RunConfig(
        genome=args.genome,
        mirna_bed=args.mirna_bed,
        gtf=args.gtf,
        iso=args.iso,
        iso_type=args.iso_type,
        iso_primers=args.iso_primers,
        iso_barcode_set=args.iso_barcode_set,
        iso_min_rq=args.iso_min_rq,
        rna=args.rna,
        strandedness=args.strandedness,
        fastp_clean=args.fastp_clean,
        extend_hairpin=args.extend_hairpin,
        min_flank=args.min_flank,
        allow_partial=args.allow_partial,
        multi_hairpin=args.multi_hairpin,
        quant=args.quant,
        quant_backend=args.quant_backend,
        out=getattr(args, "out", Path("pripipe_out")),
        tmp_dir=args.tmp_dir,
        threads=args.threads,
        sqanti3_dir=args.sqanti3_dir,
        existing_catalog=args.existing_catalog,
        hisat2_index=args.hisat2_index,
        catalog_sources=("iso,rna" if args.catalog_sources == "both" else args.catalog_sources),
        collapse_strategy=args.collapse_strategy,
        iso_singletons=args.iso_singletons,
        cupcake_dir=args.cupcake_dir,
        min_transcript_length=args.min_transcript_length,
        max_transcript_length=args.max_transcript_length,
        max_mono_exon_span=args.max_mono_exon_span,
        drop_classes=args.drop_classes,
        dedup=args.dedup,
        dedup_tss_tolerance=args.dedup_tss_tolerance,
        dedup_pas_tolerance=args.dedup_pas_tolerance,
        cdhit_identity=args.cdhit_identity,
        from_stage=getattr(args, "from_stage", None),
        to_stage=getattr(args, "to_stage", None),
        keep_intermediate=getattr(args, "keep_intermediate", False),
    )


def _print_report(rep) -> None:
    for line in rep.notes:
        print(f"[note]    {line}")
    for line in rep.warnings:
        print(f"[warning] {line}")
    for line in rep.errors:
        print(f"[ERROR]   {line}", file=sys.stderr)
    print(f"\nvalidation: {'PASS' if rep.ok else 'FAIL'}"
          f"  ({len(rep.errors)} errors, {len(rep.warnings)} warnings)")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="pripipe",
        description="pripipe - plant pri-miRNA annotation from Iso-Seq and/or RNA-Seq",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--version", action="version", version=f"pripipe {__version__}")
    sub = ap.add_subparsers(dest="command", required=True, metavar="{validate,run}")

    v = sub.add_parser("validate",
                       help="check every input and stop (no files are written except --json)",
                       description="Run all input gates and exit. Errors stop a run; warnings "
                                   "only reduce power.",
                       formatter_class=argparse.RawDescriptionHelpFormatter)
    _common(v)
    v.add_argument("--json", type=Path, metavar="FILE",
                   help="write the full validation report (errors, warnings, facts) as JSON")

    r = sub.add_parser("run",
                       help="validate then run the pipeline",
                       description="Run the pipeline stages, each writing <stage>/<stage>.done so "
                                   "a rerun skips what already succeeded.",
                       epilog=EPILOG,
                       formatter_class=argparse.RawDescriptionHelpFormatter)
    _common(r)
    r.add_argument("--from", dest="from_stage", choices=STAGES, metavar="STAGE",
                   help="start at this stage (default: iso)")
    r.add_argument("--to", dest="to_stage", choices=STAGES, metavar="STAGE",
                   help="stop after this stage (default: report)")
    r.add_argument("--keep-intermediate", action="store_true",
                   help="keep the scratch directory after a successful run")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    cfg = _config(args)
    rep = validate(cfg)
    _print_report(rep)
    if getattr(args, "json", None):
        args.json.write_text(rep.to_json(), encoding="utf-8")
    if not rep.ok:
        return 2
    if args.command == "validate":
        return 0

    cfg.prepare_dirs()
    from .stages import run_stages  # imported late so `validate` needs no stage deps

    return run_stages(cfg, rep, cfg.stages_to_run())


if __name__ == "__main__":
    raise SystemExit(main())
