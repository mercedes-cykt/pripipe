# testdata（从真实 Scaevola 课题裁剪）

窗口 `Chr02:19,354,885-21,754,984`（约 2.4 Mb），围绕 `Sta-MIRN163a`。

| 文件 | 内容 | 来源 |
|---|---|---|
| `genome.fa` | 2.4 Mb 基因组切片（header 改名为 Chr02） | `genome/Sta_v4.fasta` |
| `mirna_pre.bed` | 4 个 pre-miRNA（坐标已平移） | `annotation-basic-info.bed` |
| `ref.gtf` | 窗口内 880 个注释特征（坐标已平移） | `genome/Sta_annotation_20260906.gtf` |
| `iso.fa` | 115 条窗口内 Iso-Seq 转录本（`--iso-type hq` 用） | `sqanti3_out/scaevola_corrected.fasta` |
| `rna_1/2.fq.gz` | 1,208 对窗口内 reads | `sr_bams/flowerAA.sorted.bam` |
| `samples.tsv` | 单样本清单 | — |
| `expected.json` | CI 的基线数字 | `ci/run_ci.sh` 实测 |

所有坐标都做了同一位移，因此 BED / GTF / FASTA 与切片基因组自洽。
重新生成：`python3 build_testdata.py`（脚本需项目原始路径）。
