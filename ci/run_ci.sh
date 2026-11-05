#!/bin/bash
# ci/run_ci.sh — end-to-end smoke test on testdata/ (~2-3 min)
set -euo pipefail
export PATH=/nfs/users/zhaojw/miniconda3/envs/pri-miRNA-pipe/bin:$PATH
PY=/nfs/users/zhaojw/miniconda3/envs/pri-miRNA-pipe/bin/python3
ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-$ROOT/ci_out}
cd "$ROOT"
rm -rf "$OUT"

echo "=== pripipe CI: $(date) ==="
$PY -m pripipe.cli run \
  --genome testdata/genome.fa \
  --mirna-bed testdata/mirna_pre.bed \
  --gtf testdata/ref.gtf \
  --iso testdata/iso.fa --iso-type hq \
  --rna testdata/samples.tsv \
  --out "$OUT" --threads 8 2>&1 | tail -12

echo
echo "=== 断言 ==="
$PY - "$OUT" <<'PY'
import json, sys
from pathlib import Path

out = Path(sys.argv[1])
pri = out / "pri/pri_mirna.tsv"
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
sys.exit(0 if all(checks.values()) else 1)
PY
