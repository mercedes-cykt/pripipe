"""Run configuration and stage layout for pripipe."""

from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Optional

IsoType = Literal["subreads", "ccs", "flnc", "hq"]
Strandedness = Literal["auto", "rf", "fr", "unstranded"]
QuantBackend = Literal["auto", "kallisto", "stringtie"]

STAGES = [
    "iso",       # ccs/lima/refine/cluster2/align/collapse by --iso-type
    "rna",       # fastq -> hisat2, or accept aligned BAM
    "assemble",  # stringtie per-sample + merge
    "catalog",   # sqanti3 (mandatory) -> normalized transcript catalog
    "anchor",    # hairpin x transcript anchoring (+ intermediates, union)
    "quant",     # kallisto (FASTQ) or stringtie -e -B (BAM only)
    "report",    # pri_mirna.tsv, expression matrices, manifest.json
]


@dataclass
class RunConfig:
    # --- required inputs -------------------------------------------------
    genome: Path
    mirna_bed: Path

    # --- optional reference ---------------------------------------------
    gtf: Optional[Path] = None

    # --- Iso-Seq ---------------------------------------------------------
    iso: list[Path] = field(default_factory=list)
    iso_type: IsoType = "subreads"
    iso_primers: Optional[Path] = None
    iso_barcode_set: Optional[Path] = None
    iso_min_rq: float = 0.9
    iso_collapse_c: float = 0.99
    iso_collapse_i: float = 0.95
    collapse_strategy: Literal["cupcake", "isoseq3"] = "cupcake"
    iso_singletons: bool = True        # keep unclustered FLNC (high recall)
    cupcake_dir: Path = Path("/nfs/users/zhaojw/miniconda3/envs/isoseq/bin")

    # --- RNA-Seq ---------------------------------------------------------
    rna: Optional[Path] = None          # sample sheet (see SPEC.md)

    # --- shared options --------------------------------------------------
    strandedness: Strandedness = "auto"
    extend_hairpin: int = 0
    min_flank: int = 50
    allow_partial: bool = False
    multi_hairpin: Literal["per-hairpin", "merge"] = "per-hairpin"

    # --- quantification --------------------------------------------------
    quant: bool = True
    quant_backend: QuantBackend = "auto"

    # --- runtime ---------------------------------------------------------
    out: Path = Path("pripipe_out")
    tmp_dir: Path = Path("./tmp")
    threads: int = 8
    keep_intermediate: bool = False

    # --- reuse an existing run -------------------------------------------
    existing_catalog: Optional[Path] = None   # SQANTI3 output dir to reuse
    hisat2_index: Optional[Path] = None       # existing hisat2 index prefix
    catalog_sources: Optional[str] = None     # e.g. 'iso' / 'rna' / 'iso,rna'

    # --- external dependency (pinned contract) ---------------------------
    sqanti3_dir: Path = Path("/nfs/users/zhaojw/software/SQANTI3")
    sqanti3_pin: str = "v5.5.4"

    # --- execution control ----------------------------------------------
    from_stage: Optional[str] = None
    to_stage: Optional[str] = None

    # ------------------------------------------------------------------
    @property
    def mode(self) -> str:
        if self.iso and self.rna:
            return "both"
        if self.iso:
            return "iso"
        return "rna"

    @property
    def is_iso_from_reads(self) -> bool:
        """True when primers are required (raw subreads or CCS)."""
        return self.iso_type in ("subreads", "ccs")

    def stage_dir(self, stage: str) -> Path:
        return self.out / stage

    def resolve_quant_backend(self) -> str:
        """auto: FASTQ available -> kallisto, BAM only -> stringtie."""
        if self.quant_backend != "auto":
            return self.quant_backend
        return "kallisto" if _has_fastq(self.rna) else "stringtie"

    def stages_to_run(self) -> list[str]:
        selected = list(STAGES)
        if self.mode == "rna":
            selected = [s for s in selected if s != "iso"]
        if self.mode == "iso":
            # no RNA-Seq: nothing to align/assemble/quantify
            selected = [s for s in selected if s not in ("rna", "assemble", "quant")]
        if not self.quant:
            selected = [s for s in selected if s != "quant"]
        if self.from_stage:
            selected = selected[selected.index(self.from_stage):]
        if self.to_stage:
            selected = selected[: selected.index(self.to_stage) + 1]
        return selected

    def prepare_dirs(self) -> None:
        self.out.mkdir(parents=True, exist_ok=True)
        self.tmp_dir.mkdir(parents=True, exist_ok=True)
        for stage in STAGES:
            self.stage_dir(stage).mkdir(parents=True, exist_ok=True)

    def clean_tmp(self) -> None:
        if self.tmp_dir.exists() and not self.keep_intermediate:
            shutil.rmtree(self.tmp_dir, ignore_errors=True)


def _has_fastq(sheet: Optional[Path]) -> bool:
    """A sample sheet counts as FASTQ-bearing if any r1 column is non-empty."""
    if sheet is None or not sheet.exists():
        return False
    header, rows = _read_sheet(sheet)
    if "r1" not in header:
        return False
    idx = header.index("r1")
    return any(len(r) > idx and r[idx] not in ("", "-", ".") for r in rows)


def _read_sheet(sheet: Path) -> tuple[list[str], list[list[str]]]:
    with sheet.open(encoding="utf-8-sig") as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip() and not ln.startswith("#")]
    header = lines[0].split("\t")
    rows = [ln.split("\t") for ln in lines[1:]]
    return header, rows
