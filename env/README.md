# pri-miRNA-pipe 统一运行环境

> 建立：2026-11-04 ｜ 用途：Scaevola taccada（草海桐）pri-miRNA 全长转录组课题全流程
> 位置：`/nfs/users/zhaojw/miniconda3/envs/pri-miRNA-pipe`（3.1 GB）
> 本目录：`/nfs/users/zhaojw/project/Scaevola_taccada_miRNA/env/`

## 1. 启用法

```bash
source /nfs/users/zhaojw/miniconda3/etc/profile.d/conda.sh
conda activate pri-miRNA-pipe
```

**必须 activate，不要只改 PATH**：SQANTI3 内部用 `/bin/sh` 调 `gffread`，只把 `bin` 加进 PATH 会报 `gffread: not found`。

## 2. 包含工具

| 环节 | 工具 | 版本 |
|---|---|---|
| Iso-Seq 上游 | ccs | 6.4.0 |
| | lima | 26.2.1 |
| | isoseq3 | 26.2.0 |
| | minimap2 | 2.28-r1209 |
| 短读 | hisat2 | 2.2.1 |
| | fastp | 0.23.4 |
| | samtools | 1.21 |
| 组装 | stringtie | 3.0.3（支持 `--mix`） |
| | gffcompare | 0.12.10 |
| | gffread | 0.12.7 |
| 定量 | kallisto | 0.51.1 |
| 区间 | bedtools | 2.31.1 |
| 转录本 QC | SQANTI3（源码 `../..//software/SQANTI3`，v5.5.4） | — |
| 统计/出图 | R 4.3.3（data.table / ggplot2 / Gviz）/ pandoc / gmap | — |
| 基础 | python 3.11.13（numpy 1.26.4 / pandas 2.2.3 / pysam 0.22.1 / bx-python 0.11.0） | — |

**不含**：miRDeep-P2（继续用 `mirdp2` env）、Trinity（按课题决定不使用）。

## 3. 重建方式

```bash
# 方式 A：按环境名（可能解出更新版本）
conda env create -f environment.yml

# 方式 B：逐包精确版本（推荐复现）
conda create -n pri-miRNA-pipe --file spec-explicit.txt
```

## 4. 自检

```bash
bash verify_env.sh
```

预期：全部二进制打印版本；`stringtie --mix` 计数 ≥1；SQANTI3 `sqanti3_qc.py --help` 正常。

## 5. 已知差异

- `lima 26.2.1` / `isoseq3 26.2.0` 比 2026-09 原始环境（2.13.0 / 4.3.0）新。现有 Iso-Seq 下游产物（46,313 HQ → 26,430 转录本）已冻结，不受影响；若从原始 subreads 重跑上游，FLNC 数字可能略有差异。
- `NOISeq` 未安装（SQANTI3 v5.5.4 不使用；6.x 才需要）。
