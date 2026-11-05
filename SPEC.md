# pripipe 设计规格（v0.1，冻结中）

植物 pri-miRNA 获取工具：**输入转录本 + pre-miRNA 坐标 → 输出哪些转录本是 pri-miRNA**。
不做 miRNA 鉴定、不做 PEI/isomiR 等下游分析；有 RNA-Seq 时输出 pri 表达。

## 1. 输入契约

| 参数 | 必需性 | 默认 | 说明 |
|---|---|---|---|
| `--genome` | 必需 | — | FASTA；`.fai` 缺失自动 `samtools faidx` |
| `--mirna-bed` | 必需 | — | 前体 BED6；`name` 唯一、`strand` ∈ `+/-` |
| `--iso` | 与 `--rna` 至少一个 | — | 可多值（多 SMRT Cell） |
| `--iso-type` | 否 | `subreads` | `subreads\|ccs\|flnc\|hq` |
| `--iso-primers` | `subreads`/`ccs` 时必需 | — | lima 用；缺失在 validate 阶段直接失败 |
| `--iso-barcode-set` | 否 | — | refine 用（多重混样） |
| `--iso-min-rq` | 否 | 0.9 | ccs |
| `--rna` | 与 `--iso` 至少一个 | — | 样本清单 TSV |
| `--gtf` | 否 | — | GTF/GFF3；缺失降级 + 警告 |
| `--strandedness` | 否 | `auto` | `auto\|rf\|fr\|unstranded` |
| `--quant` / `--no-quant` | 否 | on | |
| `--quant-backend` | 否 | `auto` | `auto\|kallisto\|stringtie` |
| `--extend-hairpin` | 否 | 0 | bp |
| `--min-flank` | 否 | 50 | bp |
| `--allow-partial` | 否 | false | |
| `--multi-hairpin` | 否 | `per-hairpin` | `per-hairpin\|merge` |
| `--threads` / `--tmp-dir` / `--out` | 否 | 8 / `./tmp` / `./pripipe_out` | 禁止 `/tmp` |

模式由输入自动推断：仅 `--iso` → `iso`；仅 `--rna` → `rna`；都给 → `both`。

### 样本清单（`--rna`）

必须有表头；列名固定：

```
sample_id   r1          r2          bam         group
flowerAA    f1.fq.gz    f2.fq.gz    -           flower
rootAA      -           -           rootAA.bam  root
```

- `r1` 与 `bam` 至少给一个；两者都给 = BAM 用于组装、FASTQ 用于 kallisto 定量；
- SE 时 `r2` 写 `-`；
- `group` 缺失 → 告警继续（只出 per-sample）。

## 2. 阶段流水

| 阶段 | 说明 |
|---|---|
| `iso` | 按 `--iso-type`：ccs → lima → refine → cluster2 → minimap2 → isoseq3 collapse → gffread -T |
| `rna` | FASTQ → hisat2（或直接用 BAM） |
| `assemble` | StringTie 逐样本 + merge（有 GTF 则 `-G`） |
| `catalog` | **SQANTI3 必选**（源码 v5.5.4）→ 归一化 catalog + 并集去冗余 |
| `anchor` | hairpin × 转录本双路挂靠 + 中间体 + 证据分级 |
| `quant` | `auto`：有 FASTQ → kallisto；仅 BAM → `stringtie -e -B` |
| `report` | `pri_mirna.tsv` + 表达矩阵 + manifest |

每阶段产出 `<stage>.done`，支持 `--from/--to` 与断点续跑。

## 3. 输出（schema 1.0）

```
out/
├── catalog/transcripts.gtf
├── pri/pri_mirna.tsv
├── pri/pri_expression_tpm.tsv
├── pri/pri_expression_counts.tsv
├── qc/flagstat_summary.tsv
├── qc/validation.json
└── manifest.json
```

`pri_mirna.tsv` 列：`schema_version, locus_id, hairpin_id, chrom, hp_start, hp_end, strand,
hp_length, transcript_id, transcript_source, transcript_start, transcript_end, exon_count,
anchor_type, flank_5p, flank_3p, evidence_level, n_isoforms, n_samples_support,
genomic_env, class_code, expr_backend, pri_sum_tpm, pri_max_tpm, pri_sum_counts, n_fl`

## 4. 校验闸门（`pripipe validate`）

| 编号 | 检查 | 判定 |
|---|---|---|
| G0 | 基因组 + `.fai` | 自动生成；失败 = error |
| G1 | 染色体命名（BED/GTF/BAM vs fai） | 全不匹配 = error；部分缺失 = warning |
| G2 | BED6 合法（0-based、边界、strand、ID 唯一、长度 50–1500） | 非法 = error；长度越界 = warning |
| G3 | BAM 剪接性（CIGAR 含 `N`） | 抽样 20k reads 中 0 条 = error |
| G4 | BAM 完整性（保留 unmapped、MAPQ 过滤痕迹） | warning |
| G5 | hairpin 相互重叠 | warning（多顺聚族合法） |
| G6 | `subreads`/`ccs` 必须有 primers | error |
| G7 | SQANTI3 存在 + 版本记录 | error / warning |

## 5. 依赖契约

- **SQANTI3 必选**：源码 `v5.5.4`（git `v5.5.4-19-g4f2bf55`），路径 `/nfs/users/zhaojw/software/SQANTI3`；manifest 记录 `git describe`；
- 统一环境 `pri-miRNA-pipe`（`env/environment.yml` + `env/spec-explicit.txt`）；
- 定量后端：`kallisto 0.51.1`、`stringtie 3.0.3`（`-e -B`）。

## 6. 判据默认值（2026-11-04 定稿）

| 判据 | 取值 | 语义 |
|---|---|---|
| hairpin 长度合理区间 | `[50, 500]` nt | 越界只 warning，不阻断 |
| hairpin `name` 重复 | **error** | 名字是主键，不自动加后缀（否则与原注释对不上） |
| BAM 剪接性抽样 | 前 20,000 条 mapped reads | 0 条含 `N` 时：**≥5000 条 sampled → error**；**<5000 → warning**（样本量不足以区分"没做剪接比对"与"数据本身无跨剪接 reads"，如合成/测试数据） |
| MAPQ 过滤启发式 | 仅当 sampled ≥ 5000 时启用 | 全 0 或全非 0 各一条 warning |
| 染色体命名不全匹配 | warning | 支持近缘物种映射 |
| unmapped reads = 0 | warning | 定量分母截断 |
