#!/usr/bin/env bash
# verify_env.sh — pripipe environment self-check.
#
# Portable: uses the active conda environment when there is one, so the same
# script runs on a laptop, on the lab server and on the CI runner.
#
#   PRI_ENV_PREFIX   explicit prefix to check (wins over everything)
#   CONDA_PREFIX     used automatically when a conda env is activated
#   SQANTI3_DIR      SQANTI3 checkout; skipped when neither the variable nor
#                    the lab default exists
set -uo pipefail

if [ -n "${PRI_ENV_PREFIX:-}" ]; then
  V=$PRI_ENV_PREFIX
elif [ -n "${CONDA_PREFIX:-}" ] && [ -x "$CONDA_PREFIX/bin/python3" ]; then
  V=$CONDA_PREFIX
else
  V=/nfs/users/zhaojw/miniconda3/envs/pri-miRNA-pipe
fi
B=$V/bin

if [ ! -d "$B" ]; then
  echo "environment prefix not found: $V" >&2
  echo "set PRI_ENV_PREFIX, or activate the conda environment first" >&2
  exit 1
fi
export PATH=$B:$PATH
fail=0

# $1 = label, $2 = command. Success means exit status 0 *and* some output;
# a "No such file or directory" message must never count as a pass.
chk() {
  local out status
  out=$(eval "$2" 2>&1); status=$?
  # first non-blank line: several tools print a leading blank line
  out=$(printf '%s\n' "$out" | grep -m1 -v '^[[:space:]]*$' || true)
  if [ "$status" -ne 0 ] || [ -z "$out" ]; then
    echo "  [FAIL] $1"
    echo "        cmd: $2"
    echo "        exit: $status, output: ${out:-<empty>}"
    fail=$((fail + 1))
  else
    echo "  [ OK ] $1 : $out"
  fi
}

echo "== prefix =="
echo "  $V"

echo "== binaries =="
chk ccs        "ccs --version"
chk lima       "lima --version"
chk isoseq3    "isoseq3 --version"
chk minimap2   "minimap2 --version"
chk samtools   "samtools --version"
chk hisat2     "hisat2 --version"
chk fastp      "fastp --version"
chk stringtie  "stringtie --version"
chk gffcompare "gffcompare --version"
chk gffread    "gffread --version"
chk kallisto   "kallisto version"
chk bedtools   "bedtools --version"
chk python3    "python3 --version"
chk R          "Rscript -e 'cat(as.character(getRversion()))'"

echo "== features =="
n=$(stringtie --help 2>&1 | grep -c -- '--mix' || true)
if [ "${n:-0}" -ge 1 ]; then
  echo "  [ OK ] stringtie --mix"
else
  echo "  [FAIL] stringtie lacks --mix"
  fail=$((fail + 1))
fi

chk "python deps" "python3 -c 'import numpy, pandas, pysam, bx; print(\"py-deps\", numpy.__version__, pandas.__version__)'"

chk "R deps" "Rscript -e 'stopifnot(requireNamespace(\"data.table\", quietly=TRUE), requireNamespace(\"ggplot2\", quietly=TRUE)); cat(\"R-deps OK\")'"

SQ3=${SQANTI3_DIR:-}
if [ -z "$SQ3" ] && [ -f /nfs/users/zhaojw/software/SQANTI3/sqanti3_qc.py ]; then
  SQ3=/nfs/users/zhaojw/software/SQANTI3
fi
if [ -n "$SQ3" ] && [ -f "$SQ3/sqanti3_qc.py" ]; then
  if python3 "$SQ3/sqanti3_qc.py" --version >/dev/null 2>&1; then
    echo "  [ OK ] sqanti3_qc.py"
  else
    echo "  [FAIL] sqanti3_qc.py (python3 $SQ3/sqanti3_qc.py --version)"
    fail=$((fail + 1))
  fi
elif [ -n "$SQ3" ]; then
  echo "  [FAIL] SQANTI3_DIR set but $SQ3/sqanti3_qc.py is missing"
  fail=$((fail + 1))
else
  echo "  [SKIP] SQANTI3 (set SQANTI3_DIR to check it)"
fi

echo "== result =="
if [ "$fail" -eq 0 ]; then
  echo "ALL PASS"
else
  echo "FAILED: $fail"
fi
exit "$fail"
