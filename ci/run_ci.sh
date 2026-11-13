#!/usr/bin/env bash
# End-to-end smoke test on testdata/ (about a minute once the environment exists).
#
# Environment variables:
#   PYTHON        interpreter to use            (default: python3)
#   SQANTI3_DIR   SQANTI3 source checkout       (default: the lab path)
#   CUPCAKE_DIR   directory with Cupcake's python + collapse_isoforms_by_sam.py
#
# Failures are echoed as GitHub "::error::" annotations so that the reason shows
# up in the run summary without needing to download the raw log.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-$ROOT/ci_out}
PY=${PYTHON:-python3}
SQ3=${SQANTI3_DIR:-/nfs/users/zhaojw/software/SQANTI3}
CUPCAKE=${CUPCAKE_DIR:-/nfs/users/zhaojw/miniconda3/envs/isoseq/bin}
LOG=${TMPDIR:-/tmp}/pripipe_e2e.log

cd "$ROOT"
rm -rf "$OUT"

echo "=== pripipe CI: $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
echo "python : $("$PY" -c 'import sys; print(sys.executable)' 2>/dev/null || echo "$PY")"
echo "sqanti3: $SQ3"
echo "cupcake: $CUPCAKE"

fail() {
  echo "::error title=pripipe e2e::$1"
  exit 1
}

set +e
"$PY" -m pripipe.cli run \
  --genome testdata/genome.fa \
  --mirna-bed testdata/mirna_pre.bed \
  --gtf testdata/ref.gtf \
  --iso testdata/iso.fa --iso-type hq \
  --rna testdata/samples.tsv \
  --sqanti3-dir "$SQ3" \
  --cupcake-dir "$CUPCAKE" \
  --out "$OUT" --threads 8 >"$LOG" 2>&1
rc=$?
set -e

tail -12 "$LOG"
if [ "$rc" -ne 0 ]; then
  echo "--- last lines of $LOG ---"
  tail -25 "$LOG" || true
  while IFS= read -r line; do
    echo "::error::$line"
  done < <(grep -E 'Error|error|Traceback|raise |Exception|not found|No such' "$LOG" | tail -5 || true)
  fail "pipeline exited with $rc"
fi

echo
echo "=== assertions ==="
set +e
"$PY" - "$OUT" <<'PY'
import json, sys
from pathlib import Path

out = Path(sys.argv[1])
pri = out / "pri/pri_mirna.tsv"
if not pri.exists():
    print(f"::error title=pripipe e2e::missing {pri}")
    sys.exit(1)
rows = pri.read_text().splitlines()
header = rows[0].split("\t")
data = [r.split("\t") for r in rows[1:]]
stats = json.loads((out / "anchor/anchor_stats.json").read_text())
quant = json.loads((out / "quant/quant_stats.json").read_text())

checks = {
    "anchor loci > 0": stats["n_loci"] > 0,
    "anchor pairs > 0": stats["n_pairs"] > 0,
    "pri table rows == loci": len(data) == stats["n_loci"],
    "schema_version present": "schema_version" in header,
    "expression columns present": any(h.startswith("tpm_") for h in header),
    "counts matrix exists": (out / "pri/pri_expression_counts.tsv").exists(),
    "manifest exists": (out / "manifest.json").exists(),
}
print(f"loci={stats['n_loci']} pairs={stats['n_pairs']} "
      f"spliced={stats['spliced']} intronic={stats['intronic']} "
      f"backend={quant['backend']} columns={len(header)}")
for name, ok in checks.items():
    print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    if not ok:
        print(f"::error title=pripipe e2e::check failed: {name}")
sys.exit(0 if all(checks.values()) else 1)
PY
rc=$?
set -e
[ "$rc" -eq 0 ] || fail "assertions failed with $rc"

echo
echo "=== e2e OK ==="
