# pripipe

从 **Iso-Seq 和/或 RNA-Seq** 出发，输出植物 **pri-miRNA 注释**：哪些转录本是某个 pre-miRNA 的初级转录本、属于哪种结构类型、支持证据有多强；有 RNA-Seq 时同时给出位点级表达。

```
Iso-Seq ─┐
         ├─→ 转录本 catalog ─→ 挂靠(pre-miRNA × 转录本) ─→ pri-miRNA 注释 + 表达
RNA-Seq ─┘                            ↑
                              pre-miRNA BED6
```

## 它做什么 / 不做什么

**做**：Iso-Seq 上游处理、RNA-Seq 比对与组装、转录本 catalog、pre-miRNA 挂靠（spliced / intronic / 中间体）、位点级表达定量、可复现的 manifest。

**不做**：miRNA 位点鉴定（不调用 miRDeep-P2，位点由输入 BED 决定）、PEI/isomiR/进化等下游分析。miRNA 注释的质量决定本工具的上限。

## 安装

```bash
# 方式 A：按环境名重建（可能解出更新版本）
conda env create -f env/environment.yml

# 方式 B：逐包精确复现（推荐）
conda create -n pri-miRNA-pipe --file env/spec-explicit.txt

conda activate pri-miRNA-pipe
bash env/verify_env.sh        # 期望输出 ALL PASS
```

关键依赖：`ccs / lima / isoseq3 / minimap2 / samtools / hisat2 / fastp / stringtie>=3 / gffcompare / gffread / kallisto / bedtools / pysam`，外加 **SQANTI3**（源码目录，必选）与 **Cupcake**（`collapse_isoforms_by_sam.py`，`--collapse-strategy cupcake` 时需要）。版本与路径见 `env/README.md`。

> **必须 `conda activate`**，或至少用该环境的解释器启动：SQANTI3 内部用 `/bin/sh` 调 `gffread`，只改 PATH 不激活会报 `gffread: not found`。

## 快速开始

```bash
pripipe run \
  --genome Sta_v4.fasta \
  --mirna-bed mirna_pre.bed \
  --gtf Sta.gtf \
  --iso flnc.bam --iso-type flnc \
  --rna samples.tsv \
  --out results/ --threads 32
```

先跑自带的迷你数据集（约 20 秒）：

```bash
bash ci/run_ci.sh
```

## 输入

### 必需

| 输入 | 说明 |
|---|---|
| `--genome` FASTA | 参考基因组；`.fai` 缺失会 `samtools faidx` 自动生成 |
| `--mirna-bed` BED6 | pre-miRNA 位点：`chrom start end name score strand`。`name` 必须唯一（下游主键），`strand` 必须是 `+/-`；坐标 0-based 半开 |
| `--iso` 或 `--rna` | 至少一个测序数据来源 |

### 建议

| 输入 | 缺失后果 |
|---|---|
| `--gtf` | 无 `class_code`、无 `genic/intergenic` 标注、非链特异数据的链继承变弱（会告警） |

### Iso-Seq：`--iso-type` 决定从哪里开始

| `--iso-type` | 工具要做的 | 额外要求 |
|---|---|---|
| `subreads`（默认） | ccs → lima → refine → cluster2 | **`--iso-primers`** |
| `ccs` | lima → refine → cluster2 | **`--iso-primers`** |
| `flnc` | cluster2 | — |
| `hq` | 直接比对（已聚类/去冗余） | — |

之后统一：`minimap2 -ax splice:hq` → collapse（`cupcake` 或 `isoseq3`）→ `gffread -T`。

多 SMRT Cell 直接给多个 `--iso` 文件，会先 `samtools cat` 合并。

### RNA-Seq：`--rna` 样本清单

必须有表头，列名固定：

```
sample_id   r1              r2              bam         group
flowerAA    f1.fq.gz        f2.fq.gz        -           flower
rootAA      -               -               rootAA.bam  root
mix01       m1.fq.gz        m2.fq.gz        mix01.bam   mix
```

- `r1` 与 `bam` 至少给一个；**两者都给** = BAM 用于组装、FASTQ 用于 kallisto 定量；
- 单端：`r2` 写 `-`；
- `group` 缺失会告警，只输出 per-sample 列。

### 判据参数（都有默认值）

| 参数 | 默认 | 含义 |
|---|---|---|
| `--extend-hairpin` | 0 | 挂靠前把发夹扩 N bp（边界容差） |
| `--min-flank` | 50 | 侧翼短于 N bp 只**标记**，不丢弃 |
| `--allow-partial` | off | 是否报告只部分覆盖发夹的转录本 |
| `--multi-hairpin` | per-hairpin | 一个转录本含多个发夹时怎么记 |

## 输出

```
results/
├── iso/          # 上游转录本模型（work/ 里保留中间文件）
├── rna/          # BAM、strandness.txt
├── assemble/     # 逐样本 StringTie + merged
├── catalog/      # transcripts.gtf + meta.json（provenance）+ fl_counts.tsv
├── anchor/       # anchor_pairs.tsv / anchor_loci.tsv / intermediates.tsv / stats.json
├── quant/        # 位点 × 样本 TPM / counts
├── pri/
│   ├── pri_mirna.tsv           # 主表：一行一个 pri-miRNA 位点
│   ├── pri_pairs.tsv           # 证据明细：一行一个 发夹×转录本 对
│   ├── pri_expression_tpm.tsv
│   └── pri_expression_counts.tsv
├── qc/           # flagstat_summary.tsv、validation.json
└── manifest.json # 输入/参数/工具版本/阶段统计 —— 复现的最小闭环
```

`pri_mirna.tsv` 主要列：

| 列 | 含义 |
|---|---|
| `locus_id` | pre-miRNA 名（来自 BED） |
| `n_transcripts` / `n_spliced` / `n_intronic` | 支持证据的构成 |
| `anchor_type` | 该位点的主导类型 |
| `evidence_level` | `L1_iso`（有 Iso-Seq）/ `L2_rnaseq` / `L3_annotation` |
| `genomic_env` | `genic-exonic` / `genic-intronic` / `intergenic`（需 `--gtf`） |
| `class_code` | SQANTI3 的 structural_category |
| `n_intermediates` | 加工中间体条数 |
| `pri_sum_tpm` / `pri_max_tpm` / `tpm_<group>` | 位点级表达（有 RNA-Seq 时） |

## 两套基线（同一份数据、同参数，只差 singleton）

| 模式 | 参数 | 位点 | 用途 |
|---|---|---|---|
| 宽松（默认） | — | 69 | 发现新 pri，高召回 |
| 严格 | `--no-iso-singletons` | 38 | 与历史结果对齐、投稿复现 |

`--iso-singletons`（默认开）会保留无法聚类的单条 FLNC；它们既贡献低丰度 pri，也带来"只有 1 条转录本支持"的低证据位点。需要严格口径时关掉即可。

## 断点续跑

每个阶段成功后写 `<stage>/<stage>.done`，重跑自动跳过：

```bash
pripipe run ... --from catalog --to report   # 只重算挂靠与报表
pripipe run ... --existing-catalog sqanti3_out --catalog-sources iso
```

## 常见坑（都是实测踩过的）

1. **`subreads`/`ccs` 忘了 `--iso-primers`** —— validate 阶段直接报错，不用等到 lima。
2. **`lima` 输出不是裸前缀** —— 它会写 `<prefix>.<primer5p>--<primer3p>.bam`，工具已自动收集/合并。
3. **`isoseq3 refine` 是三参数** —— `<demux.bam> <primers.fasta> <flnc.bam>`。
4. **`cluster2` 默认丢 singleton** —— 低丰度 isoform 会被砍掉，用 `--iso-singletons` 控制。
5. **Cupcake 的比对参数是 `-s SAM`** —— 不是 `-b BAM`，工具内部自动转 SAM。
6. **BAM 必须剪接感知** —— CIGAR 无 `N` 会被判错（`subreads`/`ccs`/`flnc` 是原始 reads，不做此检查）。
7. **BED 带 UTF-8 BOM** —— Windows 生成的文件常见，工具已按 `utf-8-sig` 读取。
8. **`/tmp` 禁用** —— 临时文件走 `--tmp-dir`。

## 开发

```
pripipe/
├── cli.py        # 命令行
├── config.py     # 参数、阶段表、断点
├── models.py     # Hairpin / Transcript / AnchorPair / Intermediate
├── validate.py   # 七道输入闸门
├── stages/       # iso, rna, assemble, catalog, anchor, quant, report
├── utils.py      # 子进程 PATH 注入
├── tests/        # 回归夹具（anchor 复现 41/60/9/51）
├── testdata/     # 2.4 Mb 切片 + 4 发夹 + 115 转录本 + 1208 对 reads
└── ci/run_ci.sh  # 端到端冒烟（~20 s，7 项断言）
```

回归测试：

```bash
python3 -m pripipe.tests.repro_anchor --catalog <corrected.gtf> --bed <mirna.bed> \
    --span-raw <span_raw.tsv> --spliced-raw <spliced_raw.tsv> --out check/
# 期望：RESULT: MATCH（41 位点 / 60 对 / spliced 9 / intronic 51）
```

## 许可

GNU General Public License v3.0（见 `LICENSE`）。Copyright (C) 2026 mercedes-cykt <mrcdes@163.com>。

商用或闭源集成需另行联系作者授权；学术使用请引用本仓库。

## 状态

**v1.0.1**。四个 `--iso-type` 分支都已在真实数据上跑通；schema 版本 `1.0`。详见 `SPEC.md`（设计契约）与 `env/README.md`（环境）。
