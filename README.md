# pripipe

**pripipe** annotates plant **pri-miRNA** transcripts from Iso-Seq and/or RNA-Seq data.
For every pre-miRNA in a supplied BED file it reports which assembled transcripts
contain the hairpin, how those transcripts are structured, and — when RNA-Seq is
available — how strongly each locus is expressed.

```
Iso-Seq (subreads / CCS / FLNC / hq) ─┐
                                      ├─→ transcript catalog ─→ anchoring ─→ pri-miRNA loci
RNA-Seq (FASTQ or BAM) ───────────────┘            ↑
                                          pre-miRNA BED6
```

## Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Inputs](#inputs)
- [Common options](#common-options)
- [Outputs](#outputs)
- [Resuming and reusing runs](#resuming-and-reusing-runs)
- [Reproducibility](#reproducibility)
- [Runtime](#runtime)
- [License](#license)

## Requirements

Python ≥ 3.10 with `pysam`. The heavy lifting is done by external tools, which must
be on `PATH`:

| Stage | Tools |
|---|---|
| Iso-Seq | `ccs`, `lima`, `isoseq3`, `minimap2`, `samtools` |
| RNA-Seq | `fastp`, `hisat2`, `samtools` |
| Assembly | `stringtie` ≥ 3.0 |
| Catalog | `SQANTI3`, `gffread`, `gffcompare` |
| Collapse | `isoseq3 collapse` or Cupcake (`collapse_isoforms_by_sam.py`) |
| Quantification | `kallisto` (FASTQ) or `stringtie -e -B` (BAM) |
| Intervals | `bedtools` |

`SQANTI3` is a mandatory dependency; it is called from its source directory
(`--sqanti3-dir`) rather than from `PATH`.

A tested environment specification is provided in `env/`:

```
env/environment.yml      # rebuild by name
env/spec-explicit.txt    # exact pins, recommended for reproducibility
env/verify_env.sh        # self-check, prints ALL PASS when everything is wired up
```

## Installation

```bash
git clone https://github.com/mercedes-cykt/pripipe.git
cd pripipe

conda env create -f env/environment.yml
conda activate pri-miRNA-pipe
bash env/verify_env.sh

pip install -e .          # provides the `pripipe` command
pripipe --version
```

Activate the environment before running. External tools invoked through
`/bin/sh` (for example `gffread` inside SQANTI3) are only visible once the
environment is active.

## Quick start

A small end-to-end dataset ships in `testdata/` (a 2.4 Mb genome slice):

```bash
bash ci/run_ci.sh
```

A typical run on real data:

```bash
pripipe run \
  --genome  Sta_v4.fasta \
  --mirna-bed mirna_pre.bed \
  --gtf     Sta.gtf \
  --iso     flnc.bam --iso-type flnc \
  --rna     samples.tsv \
  --out     results/ \
  --threads 32
```

Inputs are validated first; the run stops with a non-zero exit code if a gate
fails. Use `pripipe validate` to check inputs without running anything.

## Inputs

### Required

| Input | Description |
|---|---|
| `--genome` | Reference genome FASTA. A `.fai` index is created automatically when missing. |
| `--mirna-bed` | pre-miRNA loci in BED6: `chrom start end name score strand`. Coordinates are 0-based half-open; `name` must be unique (it is the downstream key) and `strand` must be `+` or `-`. |
| `--iso` and/or `--rna` | At least one sequencing input. |

### Optional

| Input | Effect when omitted |
|---|---|
| `--gtf` | No `class_code`, no `genic`/`intergenic` annotation, weaker strand inheritance for unstranded libraries. A warning is emitted. |

### Iso-Seq

`--iso-type` selects the entry point:

| `--iso-type` | Steps performed | Extra requirement |
|---|---|---|
| `subreads` (default) | ccs → lima → refine → cluster2 | `--iso-primers` |
| `ccs` | lima → refine → cluster2 | `--iso-primers` |
| `flnc` | cluster2 | — |
| `hq` | alignment only (already clustered) | — |

All paths continue with `minimap2 -ax splice:hq`, redundancy removal and
`gffread -T`. Multiple `--iso` files are concatenated, so multi-cell libraries
can be passed directly.

### RNA-Seq

`--rna` takes a TSV sample sheet with a header:

```
sample_id   r1          r2          bam         group
flowerAA    f1.fq.gz    f2.fq.gz    -           flower
rootAA      -           -           rootAA.bam  root
mix01       m1.fq.gz    m2.fq.gz    mix01.bam   mix
```

- at least one of `r1` and `bam` is required per row;
- supplying both means the BAM is used for assembly and the FASTQ for quantification;
- single-end runs put `-` in `r2`;
- `group` is optional and enables per-group means.

## Common options

| Option | Default | Description |
|---|---|---|
| `--iso-type` | `subreads` | Iso-Seq entry point (see above) |
| `--collapse-strategy` | `cupcake` | Redundancy removal: `cupcake` or `isoseq3` |
| `--iso-singletons` / `--no-iso-singletons` | on | Keep FLNC reads that cannot be clustered. Higher recall with singletons on; stricter consensus with them off. |
| `--quant-backend` | `auto` | `auto` selects `kallisto` when FASTQ is present and `stringtie -e -B` otherwise |
| `--no-quant` | — | Annotation only, skip quantification |
| `--extend-hairpin` | 0 | Grow each hairpin by N bp before anchoring (boundary tolerance) |
| `--min-flank` | 50 | Flag transcripts with less than N bp flanking the hairpin (flagged, never dropped) |
| `--allow-partial` | off | Also report transcripts that cover the hairpin only partially |
| `--threads` | 8 | Threads passed to the external tools |
| `--from` / `--to` | — | Restrict the run to a range of stages |
| `--existing-catalog` | — | Reuse a finished SQANTI3 output directory |

`pripipe run --help` lists every option with its default.

Stages, in order: `iso → rna → assemble → catalog → anchor → quant → report`.
Stages that do not apply to the given inputs are skipped automatically.

## Outputs

```
results/
├── catalog/     transcripts.gtf, meta.json, fl_counts.tsv
├── anchor/      anchor_pairs.tsv, anchor_loci.tsv, intermediates.tsv, anchor_stats.json
├── quant/       per-sample quantification and matrices
├── pri/
│   ├── pri_mirna.tsv              one row per pri-miRNA locus
│   ├── pri_pairs.tsv              one row per hairpin x transcript pair
│   ├── pri_expression_tpm.tsv
│   └── pri_expression_counts.tsv
├── qc/          validation.json, flagstat_summary.tsv
└── manifest.json
```

`pri_mirna.tsv` key columns:

| Column | Meaning |
|---|---|
| `locus_id` | pre-miRNA name from the BED file |
| `n_transcripts`, `n_spliced`, `n_intronic` | supporting evidence |
| `anchor_type` | dominant structure of the locus |
| `evidence_level` | `L1_iso`, `L2_rnaseq` or `L3_annotation` |
| `genomic_env` | `genic-exonic`, `genic-intronic`, `intergenic` or `unknown` |
| `class_code` | SQANTI3 structural category |
| `n_intermediates` | processing intermediates observed at the locus |
| `pri_sum_tpm`, `pri_max_tpm`, `tpm_<group>` | locus-level expression (RNA-Seq runs only) |

`intermediates.tsv` lists transcripts whose end stops at a hairpin boundary,
i.e. 5' or 3' remnants of processing.

## Resuming and reusing runs

Each stage writes `<stage>/<stage>.done`; rerunning skips completed stages.

```bash
# recompute anchoring and reporting only
pripipe run ... --from catalog --to report

# reuse an existing SQANTI3 run
pripipe run ... --existing-catalog sqanti3_out --catalog-sources iso
```

## Reproducibility

`manifest.json` records the input paths and sizes, every parameter, the stage
statistics, the validation warnings and the tool versions (including
`git describe` of the SQANTI3 checkout). Together with the versioned environment
in `env/` it is sufficient to reproduce a run.

`pri_mirna.tsv` carries a `schema_version` column so that tables from different
releases can be told apart.

## Runtime

| Dataset | Resources | Wall time |
|---|---|---|
| `testdata/` (2.4 Mb slice, 115 transcripts, 1,208 read pairs) | 8 threads | ~20 s |
| Iso-Seq FLNC library, 46k FLNC reads, with Cupcake collapse | 24 threads | ~13 min |

Iso-Seq runs that start from raw subreads are dominated by `ccs`.

## License

GNU General Public License v3.0. See [LICENSE](LICENSE).

Copyright (C) 2026 mercedes-cykt &lt;mrcdes@163.com&gt;
