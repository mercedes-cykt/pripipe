#!/usr/bin/env python3
"""Regression harness for the anchoring stage.

Reproduces the 2026-09 numbers on Scaevola taccada:
    expected: 41 loci / 60 pairs / spliced 9 / intronic 51
    intermediates (strict): 6 isoforms of one locus

Usage:
    python -m pripipe.tests.repro_anchor \
        --catalog  sqanti3_out/scaevola_corrected.gtf \
        --bed      annotation-basic-info.bed \
        --span-raw pri_mirna/span_raw.tsv \
        --spliced-raw pri_mirna/spliced_raw.tsv \
        --out      /path/to/check
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from ..stages.anchor import anchor_pri, load_hairpins_bed, parse_gtf, write_outputs


def load_bedtools_pairs(path: Path) -> set[tuple[str, str]]:
    """bedtools -wa -wb output: A(6) + B(6) -> (hairpin name, transcript id)."""
    out: set[tuple[str, str]] = set()
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 10:
            continue
        out.add((f[3], f[9]))
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", required=True, type=Path)
    ap.add_argument("--bed", required=True, type=Path)
    ap.add_argument("--span-raw", type=Path)
    ap.add_argument("--spliced-raw", type=Path)
    ap.add_argument("--out", type=Path, default=Path("anchor_check"))
    args = ap.parse_args(argv)

    pins = load_hairpins_bed(args.bed)
    txs = parse_gtf(args.catalog)
    print(f"hairpins: {len(pins)}   transcripts: {len(txs)}")

    pairs, intermeds, loci, st = anchor_pri(pins, txs, min_flank=0,
                                            allow_partial=False, window=50)
    write_outputs(args.out, pairs, intermeds, loci, st)

    print("\n== computed ==")
    for k, v in st.items():
        print(f"  {k:22s} {v}")

    ok = True
    if args.span_raw and args.span_raw.exists():
        expected = load_bedtools_pairs(args.span_raw)
        got = {(p.hairpin.name, p.tx.tid) for p in pairs}
        print(f"\n== vs span_raw: expected {len(expected)} pairs, got {len(got)}")
        missing, extra = expected - got, got - expected
        print(f"  missing {len(missing)}: {sorted(missing)[:5]}")
        print(f"  extra   {len(extra)}: {sorted(extra)[:5]}")
        ok &= not missing and not extra

    if args.spliced_raw and args.spliced_raw.exists():
        expected = load_bedtools_pairs(args.spliced_raw)
        got = {(p.hairpin.name, p.tx.tid) for p in pairs if p.anchor_type == "spliced"}
        print(f"\n== vs spliced_raw: expected {len(expected)}, got {len(got)}")
        missing, extra = expected - got, got - expected
        print(f"  missing {len(missing)}: {sorted(missing)[:5]}")
        print(f"  extra   {len(extra)}: {sorted(extra)[:5]}")
        ok &= not missing and not extra

    print("\nRESULT:", "MATCH" if ok else "MISMATCH")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
