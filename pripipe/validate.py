#!/usr/bin/env python3
"""Input validation gates for pripipe.

Design rule: anything that would silently corrupt coordinates, strand, or
quantification is an *error* (hard stop). Anything that only reduces power
is a *warning*. Every finding ends up in the validation report and, later,
in manifest.json so a run can always be audited.

Gates implemented here (see SPEC.md 5):
  G1 chromosome naming consistency (BED / GTF / BAM vs genome .fai)
  G2 BED6 legality (0-based, bounds, strand, unique ids, length sanity)
  G3 BAM splice-awareness (a genome BAM must contain N in CIGAR)
  G4 BAM completeness (unmapped reads retained; MAPQ filtering hint)
  G5 overlapping / duplicated hairpins
  G6 primers required for iso_type in {subreads, ccs}
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

from .config import RunConfig, _read_sheet
from .models import Hairpin
from .utils import env_with_env_bin, missing_tool_message

MIN_HAIRPIN_LEN = 50
MAX_HAIRPIN_LEN = 500
BAM_SCAN_READS = 20000
# Below this many mapped reads the splice-awareness test cannot distinguish
# "aligner was not splice-aware" from "this data simply has no junction read"
# (synthetic test data, exon-only simulations). Demote to warning there.
MIN_SPLICE_SCAN = 5000


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------
@dataclass
class Report:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    facts: dict = field(default_factory=dict)

    def err(self, msg: str) -> None:
        self.errors.append(msg)

    def warn(self, msg: str) -> None:
        self.warnings.append(msg)

    def note(self, msg: str) -> None:
        self.notes.append(msg)

    @property
    def ok(self) -> bool:
        return not self.errors

    def to_json(self) -> str:
        return json.dumps(
            {
                "ok": self.ok,
                "errors": self.errors,
                "warnings": self.warnings,
                "notes": self.notes,
                "facts": self.facts,
            },
            ensure_ascii=False,
            indent=2,
        )


def _status(msg: str) -> None:
    sys.stderr.write(f"[validate] {msg}\n")


def _run(cmd: list[str]) -> subprocess.CompletedProcess:
    try:
        return subprocess.run(cmd, capture_output=True, text=True, env=env_with_env_bin())
    except FileNotFoundError:
        raise RuntimeError(missing_tool_message(cmd[0])) from None


# --------------------------------------------------------------------------
# G0 genome + .fai (auto-build)
# --------------------------------------------------------------------------
def ensure_fai(cfg: RunConfig, rep: Report) -> dict[str, int]:
    genome = cfg.genome
    if not genome.exists():
        rep.err(f"genome not found: {genome}")
        return {}
    fai = Path(str(genome) + ".fai")
    if not fai.exists():
        _status(f".fai missing, building: samtools faidx {genome}")
        proc = _run(["samtools", "faidx", str(genome)])
        if proc.returncode != 0:
            rep.err(f"failed to build .fai: {proc.stderr.strip()[:300]}")
            return {}
    lengths: dict[str, int] = {}
    with fai.open(encoding="utf-8-sig") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                lengths[parts[0]] = int(parts[1])
    rep.facts["genome"] = {"path": str(genome), "n_chrom": len(lengths)}
    rep.note(f"genome: {len(lengths)} sequences")
    return lengths


# --------------------------------------------------------------------------
# G1 + G2 + G5 hairpin BED
# --------------------------------------------------------------------------
def read_hairpins(cfg: RunConfig, lengths: dict[str, int], rep: Report) -> list[Hairpin]:
    path = cfg.mirna_bed
    if not path.exists():
        rep.err(f"miRNA BED not found: {path}")
        return []

    pins: list[Hairpin] = []
    seen: dict[str, int] = {}
    bad_cols = 0
    for lineno, raw in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
            continue
        f = line.split("\t")
        if len(f) < 6:
            bad_cols += 1
            if bad_cols <= 3:
                rep.err(f"BED line {lineno}: need 6 columns (chrom start end name score strand), got {len(f)}")
            continue
        chrom, start, end, name, score, strand = f[0], f[1], f[2], f[3], f[4], f[5]
        try:
            st, en = int(start), int(end)
        except ValueError:
            rep.err(f"BED line {lineno}: non-integer coordinates ({start!r}, {end!r})")
            continue
        # G2: coordinate legality
        if st < 0 or en <= st:
            rep.err(f"BED line {lineno} ({name}): illegal interval {st}-{en}")
            continue
        if strand not in ("+", "-"):
            rep.err(f"BED line {lineno} ({name}): strand must be '+' or '-', got {strand!r}")
            continue
        # G1: naming
        if chrom not in lengths:
            rep.err(f"BED line {lineno} ({name}): chromosome {chrom!r} not in genome")
            continue
        if en > lengths[chrom]:
            rep.err(f"BED line {lineno} ({name}): end {en} exceeds {chrom} length {lengths[chrom]}")
            continue
        if name in seen:
            rep.err(f"BED line {lineno}: duplicate hairpin id {name!r} (first at line {seen[name]})")
            continue
        seen[name] = lineno
        hp = Hairpin(chrom, st, en, name, score, strand)
        if not (MIN_HAIRPIN_LEN <= hp.length <= MAX_HAIRPIN_LEN):
            rep.warn(f"hairpin {name}: length {hp.length} outside [{MIN_HAIRPIN_LEN}, {MAX_HAIRPIN_LEN}]")
        pins.append(hp)

    if bad_cols > 3:
        rep.err(f"BED: {bad_cols} lines with fewer than 6 columns in total")

    # G1 refinement: give an actionable naming hint before drowning in errors
    bed_chroms = {p.chrom for p in pins}
    if not pins and lengths and bed_chroms == set():
        first = path.read_text(encoding="utf-8-sig").splitlines()[0].split("\t")[0] if path.stat().st_size else "?"
        genome_first = next(iter(lengths))
        rep.err(
            "BED produced no usable hairpins. If the first column looks like "
            f"{first!r} while the genome uses {genome_first!r}, this is a chromosome-naming "
            "mismatch (e.g. 'Chr01' vs 'chr1'); rename the BED or the genome headers."
        )

    # G5: overlapping hairpins are legal (tandem families) but must be declared
    pairs = 0
    by_chrom: dict[str, list[Hairpin]] = {}
    for p in pins:
        by_chrom.setdefault(p.chrom, []).append(p)
    for chrom, group in by_chrom.items():
        group.sort(key=lambda x: x.start)
        for a, b in zip(group, group[1:]):
            if b.start < a.end:
                pairs += 1
    if pairs:
        rep.warn(f"{pairs} overlapping hairpin pair(s); multi-hairpin rule = {cfg.multi_hairpin}")

    rep.facts["hairpins"] = {
        "n": len(pins),
        "strand_split": {
            "+": sum(1 for p in pins if p.strand == "+"),
            "-": sum(1 for p in pins if p.strand == "-"),
        },
        "median_length": sorted(p.length for p in pins)[len(pins) // 2] if pins else 0,
    }
    rep.note(f"hairpins: {len(pins)} loci ({rep.facts['hairpins']['strand_split']})")
    return pins


# --------------------------------------------------------------------------
# G6 + Iso-Seq inputs
# --------------------------------------------------------------------------
def check_iso_inputs(cfg: RunConfig, rep: Report) -> None:
    if not cfg.iso:
        return
    for p in cfg.iso:
        if not p.exists():
            rep.err(f"Iso-Seq input not found: {p}")
            continue
        if cfg.iso_type in ("subreads", "ccs", "flnc") and p.suffix.lower() != ".bam":
            rep.warn(f"iso_type={cfg.iso_type}: expect BAM, got {p.name}")

    if cfg.is_iso_from_reads and not cfg.iso_primers:
        rep.err(
            f"iso_type={cfg.iso_type} needs primer info: pass --iso-primers primers.fasta, "
            "or use --iso-type flnc / hq if your data is already refined/clustered."
        )
    if cfg.iso_primers and not cfg.iso_primers.exists():
        rep.err(f"primer file not found: {cfg.iso_primers}")
    if cfg.iso_barcode_set and not cfg.iso_barcode_set.exists():
        rep.err(f"barcode set not found: {cfg.iso_barcode_set}")

    if cfg.iso_type == "hq":
        for p in cfg.iso:
            if p.suffix.lower() not in (".fa", ".fasta", ".fna", ".gz"):
                rep.warn(f"iso_type=hq expects FASTA, got {p.name}")
    rep.note(f"iso: {len(cfg.iso)} file(s), type={cfg.iso_type}")


# --------------------------------------------------------------------------
# RNA-Seq sample sheet
# --------------------------------------------------------------------------
@dataclass
class Sample:
    sample_id: str
    r1: str | None
    r2: str | None
    bam: str | None
    group: str | None

    @property
    def kind(self) -> str:
        if self.bam and self.r1:
            return "bam+fastq"
        return "bam" if self.bam else "fastq"


def check_samples(cfg: RunConfig, rep: Report) -> list[Sample]:
    if not cfg.rna:
        return []
    sheet = cfg.rna
    if not sheet.exists():
        rep.err(f"sample sheet not found: {sheet}")
        return []

    header, rows = _read_sheet(sheet)
    need = {"sample_id"}
    if not need.issubset(header):
        rep.err(f"sample sheet must have a header with 'sample_id' (got: {header})")
        return []
    if "r1" not in header and "bam" not in header:
        rep.err("sample sheet needs at least one of 'r1' (FASTQ) or 'bam' columns")
        return []
    if "group" not in header:
        rep.warn("no 'group' column: replicate filter and group means will be skipped")

    idx = {k: i for i, k in enumerate(header)}
    samples: list[Sample] = []
    seen: set[str] = set()
    for lineno, r in enumerate(rows, 2):
        def get(key: str) -> str | None:
            i = idx.get(key)
            if i is None or i >= len(r):
                return None
            v = r[i].strip()
            return None if v in ("", "-", ".") else v

        sid = get("sample_id")
        if not sid:
            rep.err(f"sample sheet line {lineno}: empty sample_id")
            continue
        if sid in seen:
            rep.err(f"sample sheet line {lineno}: duplicate sample_id {sid!r}")
            continue
        seen.add(sid)

        r1, r2, bam, group = get("r1"), get("r2"), get("bam"), get("group")
        if not r1 and not bam:
            rep.err(f"sample {sid}: neither FASTQ (r1) nor BAM given")
            continue
        if r2 and not r1:
            rep.err(f"sample {sid}: r2 given without r1")
            continue
        for label, path in (("r1", r1), ("r2", r2), ("bam", bam)):
            if path and not Path(path).exists():
                rep.err(f"sample {sid}: {label} not found: {path}")
        samples.append(Sample(sid, r1, r2, bam, group))

    kinds: dict[str, int] = {}
    for s in samples:
        kinds[s.kind] = kinds.get(s.kind, 0) + 1
    rep.facts["samples"] = {
        "n": len(samples),
        "kinds": kinds,
        "groups": sorted({s.group for s in samples if s.group}),
        "paired": sum(1 for s in samples if s.r2),
    }
    rep.note(f"samples: {len(samples)} ({kinds}), groups={rep.facts['samples']['groups']}")
    return samples


# --------------------------------------------------------------------------
# G3 + G4 BAM gates
# --------------------------------------------------------------------------
def check_bam(path: Path, lengths: dict[str, int], rep: Report, label: str) -> None:
    try:
        import pysam  # noqa: PLC0415  (optional at import time)
    except ImportError:
        rep.err("pysam is required for BAM checks (pip install pysam / conda install pysam)")
        return

    if not path.exists():
        return
    bai = Path(str(path) + ".bai")
    csi = Path(str(path) + ".csi")
    if not bai.exists() and not csi.exists():
        _status(f"index missing for {path.name}, building")
        try:
            pysam.index(str(path))
        except Exception as exc:  # pragma: no cover - depends on env
            rep.warn(f"{label}: cannot build BAM index ({exc})")

    try:
        af = pysam.AlignmentFile(str(path), "rb", check_sq=False)
    except Exception as exc:
        rep.err(f"{label}: cannot open BAM ({exc})")
        return

    sq = {ref: af.lengths[i] for i, ref in enumerate(af.references)}
    af.close()

    # G3a: must be aligned to this genome (soft: naming may be near-relative)
    missing = set(sq) - set(lengths)
    if lengths and missing and len(missing) == len(sq):
        rep.err(
            f"{label}: none of the {len(sq)} BAM sequences exist in the genome — "
            "this looks like a transcriptome alignment or a different assembly"
        )
        return
    if missing:
        rep.warn(f"{label}: {len(missing)} BAM sequence(s) absent from genome (near-relative assembly?)")

    spliced = 0
    mapq_zero = 0
    mapq_max = 0
    scanned = 0
    af = pysam.AlignmentFile(str(path), "rb", check_sq=False)
    for rec in af.fetch(until_eof=True):
        if rec.is_unmapped:
            continue
        scanned += 1
        if any(op == 3 for op, _ in rec.cigartuples or []):  # 3 = N
            spliced += 1
        if rec.mapping_quality == 0:
            mapq_zero += 1
        mapq_max = max(mapq_max, rec.mapping_quality)
        if scanned >= BAM_SCAN_READS:
            break
    af.close()

    if scanned and spliced == 0:
        if scanned < MIN_SPLICE_SCAN:
            rep.warn(
                f"{label}: 0 spliced read in the first {scanned} mapped reads, but that is below "
                f"{MIN_SPLICE_SCAN} — cannot prove the BAM is not splice-aware (test data / "
                "exon-only simulation?); please confirm the aligner preset yourself"
            )
        else:
            rep.err(
                f"{label}: no spliced read in the first {scanned} mapped reads — this BAM was not "
                "aligned splice-aware (or used a wrong preset); exon structure cannot be trusted"
            )
    if scanned >= MIN_SPLICE_SCAN:
        if mapq_zero == scanned:
            rep.warn(f"{label}: all {scanned} sampled reads have MAPQ=0 (multi-mapping or low quality)")
        elif mapq_zero == 0:
            rep.warn(f"{label}: no MAPQ=0 read in {scanned} sampled reads — possibly MAPQ-filtered "
                     "(quantification bias)")

    # G4: unmapped retention via samtools flagstat
    proc = _run(["samtools", "flagstat", str(path)])
    total = mapped = dup = 0
    if proc.returncode == 0:
        for line in proc.stdout.splitlines():
            m = re.match(r"^(\d+)\s+\+\s+\d+\s+in total", line)
            if m:
                total = int(m.group(1))
            m = re.match(r"^(\d+)\s+\+\s+\d+\s+mapped\s+\(", line)
            if m:
                mapped = int(m.group(1))
            m = re.match(r"^(\d+)\s+\+\s+\d+\s+duplicates", line)
            if m:
                dup = int(m.group(1))
        if total and mapped:
            unmapped = total - mapped
            if unmapped == 0:
                rep.warn(f"{label}: no unmapped read retained -> quantification denominator is truncated")
            rep.facts.setdefault("bam_qc", []).append(
                {"file": str(path), "total": total, "mapped": mapped,
                 "unmapped": unmapped, "duplicates": dup,
                 "sampled": scanned, "spliced_sampled": spliced, "mapq_max": mapq_max}
            )
    else:
        rep.warn(f"{label}: samtools flagstat failed, QC numbers unavailable")


# --------------------------------------------------------------------------
# optional GTF/GFF
# --------------------------------------------------------------------------
def check_gtf(cfg: RunConfig, lengths: dict[str, int], rep: Report) -> None:
    if not cfg.gtf:
        rep.warn(
            "no --gtf: strand inheritance for unstranded RNA-Seq falls back to splice motifs, "
            "class_code and genic/intergenic labels unavailable"
        )
        return
    if not cfg.gtf.exists():
        rep.err(f"annotation not found: {cfg.gtf}")
        return
    exon_lines = 0
    chroms: set[str] = set()
    with cfg.gtf.open(encoding="utf-8-sig", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i > 2_000_000:
                break
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 9:
                chroms.add(parts[0])
                if parts[2] == "exon":
                    exon_lines += 1
    if exon_lines == 0:
        rep.err(f"{cfg.gtf.name}: no exon feature found (is this a GTF/GFF3?)")
        return
    missing = chroms - set(lengths)
    if lengths and missing and len(missing) == len(chroms):
        rep.err(f"{cfg.gtf.name}: chromosome naming does not match genome at all")
    elif missing:
        rep.warn(f"{cfg.gtf.name}: {len(missing)} chromosome(s) absent from genome")
    rep.facts["gtf"] = {"path": str(cfg.gtf), "exon_features": exon_lines, "n_chrom": len(chroms)}
    rep.note(f"annotation: {exon_lines} exon features")


# --------------------------------------------------------------------------
# mandatory external dependency
# --------------------------------------------------------------------------
def check_sqanti3(cfg: RunConfig, rep: Report) -> None:
    qc = cfg.sqanti3_dir / "sqanti3_qc.py"
    if not qc.exists():
        rep.err(f"SQANTI3 not found at {cfg.sqanti3_dir} (expected {qc}); it is a mandatory dependency")
        return
    version = "unknown"
    proc = _run(["git", "-C", str(cfg.sqanti3_dir), "describe", "--tags"])
    if proc.returncode == 0:
        version = proc.stdout.strip()
    rep.facts["sqanti3"] = {"path": str(cfg.sqanti3_dir), "version": version, "pin": cfg.sqanti3_pin}
    if version not in ("unknown", cfg.sqanti3_pin) and not version.startswith(cfg.sqanti3_pin):
        rep.warn(f"SQANTI3 version {version} differs from pinned {cfg.sqanti3_pin}")
    rep.note(f"SQANTI3: {version}")


# --------------------------------------------------------------------------
# orchestrator
# --------------------------------------------------------------------------
def validate(cfg: RunConfig) -> Report:
    rep = Report()
    lengths = ensure_fai(cfg, rep)
    hairpins = read_hairpins(cfg, lengths, rep) if lengths else []
    check_iso_inputs(cfg, rep)
    samples = check_samples(cfg, rep)
    for s in samples:
        if s.bam:
            check_bam(Path(s.bam), lengths, rep, f"sample {s.sample_id}")
    # NOTE: Iso-Seq inputs are raw/clustered reads, not alignments, so the
    # splice-awareness gate must not run on them.
    check_gtf(cfg, lengths, rep)
    check_sqanti3(cfg, rep)

    if cfg.quant and cfg.mode == "iso":
        rep.note("quant requested but no RNA-Seq input: expression matrix will be skipped")

    rep.facts["run"] = {
        "mode": cfg.mode,
        "iso_type": cfg.iso_type,
        "quant_backend": cfg.resolve_quant_backend() if cfg.quant else "off",
        "hairpins": len(hairpins),
        "samples": len(samples),
    }
    return rep
