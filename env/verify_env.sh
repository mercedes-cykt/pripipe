#!/bin/bash
# verify_env.sh — pri-miRNA-pipe 环境自检
set -u
V=/nfs/users/zhaojw/miniconda3/envs/pri-miRNA-pipe
B=$V/bin
export PATH=$B:$PATH
fail=0

chk() {  # $1=label  $2=command
  out=$(eval "$2" 2>&1 | head -1)
  if [ -z "$out" ]; then echo "  [FAIL] $1"; fail=$((fail+1)); else echo "  [ OK ] $1 : $out"; fi
}

echo "== 二进制 =="
chk ccs        "$B/ccs --version"
chk lima       "$B/lima --version"
chk isoseq3    "$B/isoseq3 --version"
chk minimap2   "$B/minimap2 --version"
chk samtools   "$B/samtools --version"
chk hisat2     "$B/hisat2 --version"
chk fastp      "$B/fastp --version"
chk stringtie  "$B/stringtie --version"
chk gffcompare "$B/gffcompare --version"
chk gffread    "$B/gffread --version"
chk kallisto   "$B/kallisto version"
chk bedtools   "$B/bedtools --version"
chk python3    "$B/python3 --version"
chk R          "$B/Rscript -e 'cat(as.character(getRversion()))'"

echo "== 功能 =="
n=$($B/stringtie --help 2>&1 | grep -c -- '--mix')
[ "$n" -ge 1 ] && echo "  [ OK ] stringtie --mix 支持" || { echo "  [FAIL] stringtie 无 --mix"; fail=$((fail+1)); }

$B/python3 -c "import numpy,pandas,pysam,bx;print('py-deps',numpy.__version__,pandas.__version__)" \
  && echo "  [ OK ] python 依赖" || { echo "  [FAIL] python 依赖"; fail=$((fail+1)); }

$B/Rscript -e 'stopifnot(requireNamespace("data.table",quietly=TRUE),requireNamespace("ggplot2",quietly=TRUE));cat("R-deps OK\n")' \
  || { echo "  [FAIL] R 依赖"; fail=$((fail+1)); }

$B/python3 /nfs/users/zhaojw/software/SQANTI3/sqanti3_qc.py --help >/dev/null 2>&1 \
  && echo "  [ OK ] sqanti3_qc.py --help" || { echo "  [FAIL] sqanti3_qc.py"; fail=$((fail+1)); }

echo "== 结果 =="
[ $fail -eq 0 ] && echo "ALL PASS" || echo "FAILED: $fail"
exit $fail
