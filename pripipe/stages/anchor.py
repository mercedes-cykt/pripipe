#!/usr/bin/env python3
"""Hairpin x transcript anchoring: the core of pripipe.

Two questions are answered per pre-miRNA locus:

1. which transcripts fully contain the hairpin on the same strand
   (*pri-miRNA* candidates, split into `spliced` vs `intronic`), and
2. which transcripts stop right at the hairpin boundary
   (processing intermediates: 5' / 3' remnants).

The union of the Iso-Seq and RNA-Seq catalogs is used, so a candidate can be
supported by either source or both.
"""

from __future__ import annotations

import bisect
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

from ..config import RunConfig
from ..models import AnchorPair, Hairpin, Intermediate, Transcript

TX_ID_RE = re.compile(r'transcript_id\s*[=]?\s*"?([^";]+)"?')
SRC_RE = re.compile(r'pri_source\s+"([^"]+)"')
CLUSTER_TOL = 3    # bp tolerance when clustering transcript endpoints
MIN_CLUSTER = 2    # isoforms required before an intermediate is called


def _status(msg: str) -> None:
    sys.stderr.write(f"[anchor] {msg}\n")


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------
def parse_gtf(path: Path, source: str = "RNA-Seq") -> dict[str, Transcript]:
    """Parse a GTF/GFF3 into transcript models (exons 0-based half-open).

    A `pri_source` attribute (written by the catalog stage) overrides `source`
    so that a merged catalog keeps its per-transcript provenance.
    """
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    meta: dict[str, tuple[str, str, str]] = {}
    with path.open(encoding="utf-8-sig", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = TX_ID_RE.search(f[8])
            if not m:
                continue
            tid = m.group(1)
            start, end = int(f[3]) - 1, int(f[4])       # 1-based -> 0-based
            exons[tid].append((start, end))
            src = SRC_RE.search(f[8])
            meta.setdefault(tid, (f[0], f[6], src.group(1) if src else source))

    txs: dict[str, Transcript] = {}
    for tid, spans in exons.items():
        chrom, strand, tx_source = meta[tid]
        txs[tid] = Transcript(tid=tid, chrom=chrom, strand=strand,
                              exons=sorted(spans), source=tx_source)
    return txs


def load_hairpins_bed(path: Path) -> list[Hairpin]:
    pins: list[Hairpin] = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        pins.append(Hairpin(chrom=f[0], start=int(f[1]), end=int(f[2]),
                            name=f[3], score=f[4] if len(f) > 4 else "0",
                            strand=f[5] if len(f) > 5 else "+"))
    return pins


def merge_catalogs(catalogs: list[tuple[str, dict[str, Transcript]]]) -> dict[str, Transcript]:
    """Union of catalogs; identical models are merged and keep every source."""
    merged: dict[tuple, Transcript] = {}
    for source, txs in catalogs:
        for tx in txs.values():
            key = tx.chain_key()
            hit = merged.get(key)
            if hit is None:
                merged[key] = tx
                continue
            if source not in hit.source.split("+"):
                hit.source = "+".join(sorted(set(hit.source.split("+")) | {source},
                                             key=lambda s: (s != "Iso-Seq", s)))
            hit.aliases.append(tx.tid)
    out: dict[str, Transcript] = {}
    for tx in merged.values():
        out[tx.tid] = tx
    return out


def build_index(txs: dict[str, Transcript]):
    """(chrom, strand) -> (sorted starts, transcript list)."""
    idx: dict[tuple[str, str], tuple[list[int], list[Transcript]]] = {}
    buckets: dict[tuple[str, str], list[Transcript]] = defaultdict(list)
    for tx in txs.values():
        buckets[(tx.chrom, tx.strand)].append(tx)
    for key, lst in buckets.items():
        lst.sort(key=lambda t: t.start)
        idx[key] = ([t.start for t in lst], lst)
    return idx


# --------------------------------------------------------------------------
# anchoring
# --------------------------------------------------------------------------
def _flanks(hp: Hairpin, tx: Transcript) -> tuple[int, int]:
    if tx.strand == "+":
        return hp.start - tx.start, tx.end - hp.end
    return tx.end - hp.end, hp.start - tx.start


def find_pairs(pins: list[Hairpin], txs: dict[str, Transcript],
               min_flank: int, allow_partial: bool) -> list[AnchorPair]:
    idx = build_index(txs)
    pairs: list[AnchorPair] = []
    for hp in pins:
        key = (hp.chrom, hp.strand)
        if key not in idx:
            continue
        starts, lst = idx[key]
        hi = bisect.bisect_right(starts, hp.end)
        for tx in lst[:hi]:
            if tx.end < hp.start:
                continue
            contained = tx.contains(hp.start, hp.end)
            exon_idx = tx.exon_containing(hp.start, hp.end) if contained else None
            if contained:
                kind = "spliced" if exon_idx is not None else "intronic"
            elif allow_partial:
                kind = "partial"
            else:
                continue
            f5, f3 = _flanks(hp, tx)
            pairs.append(AnchorPair(hairpin=hp, tx=tx, anchor_type=kind,
                                    flank_5p=f5, flank_3p=f3, exon_index=exon_idx))
    if min_flank > 0:
        short = sum(1 for p in pairs if min(p.flank_5p, p.flank_3p) < min_flank)
        if short:
            _status(f"{short} pair(s) have a flank shorter than {min_flank} bp (flagged, not dropped)")
    return pairs


def find_intermediates(pins: list[Hairpin], txs: dict[str, Transcript],
                       pairs: list[AnchorPair], window: int) -> list[Intermediate]:
    """Transcripts whose end stops inside / next to a hairpin boundary."""
    covered = {(p.hairpin.name, p.tx.tid) for p in pairs}
    idx = build_index(txs)
    raw: list[Intermediate] = []
    for hp in pins:
        key = (hp.chrom, hp.strand)
        if key not in idx:
            continue
        starts, lst = idx[key]
        hi = bisect.bisect_right(starts, hp.end + window)
        for tx in lst[:hi]:
            if tx.end < hp.start - window or (hp.name, tx.tid) in covered:
                continue
            for kind, endpoint in (("5prime_remnant", tx.end), ("3prime_remnant", tx.start)):
                if hp.start - window <= endpoint <= hp.end + window:
                    inside = hp.start < endpoint < hp.end
                    # signed distance to the nearest hairpin boundary:
                    # positive = endpoint sits inside the hairpin, negative = outside
                    if inside:
                        offset = min(endpoint - hp.start, hp.end - endpoint)
                    else:
                        offset = -min(abs(endpoint - hp.start), abs(endpoint - hp.end))
                    raw.append(Intermediate(hairpin=hp, tx=tx, kind=kind,
                                            endpoint=endpoint, offset=offset,
                                            strict=inside))
    # cluster endpoints: a real processing site should be shared by >=2 isoforms
    groups: dict[tuple[str, str], list[Intermediate]] = defaultdict(list)
    for item in raw:
        groups[(item.hairpin.name, item.kind)].append(item)
    out: list[Intermediate] = []
    for items in groups.values():
        items.sort(key=lambda x: x.endpoint)
        cluster: list[Intermediate] = [items[0]]
        for prev, cur in zip(items, items[1:]):
            if cur.endpoint - prev.endpoint <= CLUSTER_TOL:
                cluster.append(cur)
            else:
                _flush(cluster, out)
                cluster = [cur]
        _flush(cluster, out)
    return out


def _flush(cluster: list[Intermediate], out: list[Intermediate]) -> None:
    if len(cluster) >= MIN_CLUSTER:
        for item in cluster:
            item.n_cluster = len(cluster)
        out.extend(cluster)


# --------------------------------------------------------------------------
# aggregation + output
# --------------------------------------------------------------------------
def locus_rows(pairs: list[AnchorPair], intermeds: list[Intermediate]) -> list[dict]:
    per_locus: dict[str, dict] = {}
    tx_to_hp: dict[str, set[str]] = defaultdict(set)
    for p in pairs:
        row = per_locus.setdefault(p.hairpin.name, {
            "hairpin_id": p.hairpin.name, "chrom": p.hairpin.chrom,
            "hp_start": p.hairpin.start, "hp_end": p.hairpin.end,
            "strand": p.hairpin.strand, "n_pairs": 0, "n_tx": set(),
            "spliced": 0, "intronic": 0, "partial": 0,
            "sources": set(),
            "max_flank_5p": 0, "max_flank_3p": 0, "n_intermediates": 0,
        })
        row["n_pairs"] += 1
        row["n_tx"].add(p.tx.tid)
        row[p.anchor_type] += 1
        row["sources"].update(p.tx.source.split("+"))
        row["max_flank_5p"] = max(row["max_flank_5p"], p.flank_5p)
        row["max_flank_3p"] = max(row["max_flank_3p"], p.flank_3p)
        tx_to_hp[p.tx.tid].add(p.hairpin.name)
    for item in intermeds:
        row = per_locus.get(item.hairpin.name)
        if row:
            row["n_intermediates"] += 1
    rows = []
    for row in per_locus.values():
        row["n_tx"] = len(row["n_tx"])
        row["sources"] = "+".join(sorted(row["sources"]))
        rows.append(row)
    rows.sort(key=lambda r: (r["chrom"], r["hp_start"]))
    return rows


def stats(pairs: list[AnchorPair], intermeds: list[Intermediate]) -> dict:
    kinds = defaultdict(int)
    for p in pairs:
        kinds[p.anchor_type] += 1
    poly = defaultdict(set)
    for p in pairs:
        poly[p.tx.tid].add(p.hairpin.name)
    return {
        "n_loci": len({p.hairpin.name for p in pairs}),
        "n_pairs": len(pairs),
        "spliced": kinds["spliced"],
        "intronic": kinds["intronic"],
        "partial": kinds["partial"],
        "n_polycistronic_tx": sum(1 for v in poly.values() if len(v) >= 2),
        "n_intermediates": len(intermeds),
        "intermediate_strict": sum(1 for i in intermeds if i.strict),
        "intermediate_loci": len({i.hairpin.name for i in intermeds}),
    }


PAIR_COLS = ["hairpin_id", "chrom", "hp_start", "hp_end", "strand", "hp_length",
             "transcript_id", "transcript_source", "tx_start", "tx_end", "n_exons",
             "anchor_type", "exon_index", "flank_5p", "flank_3p"]


def write_outputs(out: Path, pairs: list[AnchorPair], intermeds: list[Intermediate],
                  loci: list[dict], st: dict) -> None:
    out.mkdir(parents=True, exist_ok=True)
    with (out / "anchor_pairs.tsv").open("w", encoding="utf-8") as fh:
        fh.write("\t".join(PAIR_COLS) + "\n")
        for p in sorted(pairs, key=lambda x: (x.hairpin.chrom, x.hairpin.start, x.tx.tid)):
            fh.write("\t".join(str(v) for v in [
                p.hairpin.name, p.hairpin.chrom, p.hairpin.start, p.hairpin.end,
                p.hairpin.strand, p.hairpin.length, p.tx.tid, p.tx.source,
                p.tx.start, p.tx.end, p.tx.n_exons, p.anchor_type,
                p.exon_index if p.exon_index is not None else "NA",
                p.flank_5p, p.flank_3p]) + "\n")

    with (out / "anchor_loci.tsv").open("w", encoding="utf-8") as fh:
        cols = ["hairpin_id", "chrom", "hp_start", "hp_end", "strand", "n_pairs",
                "n_tx", "spliced", "intronic", "partial", "sources",
                "max_flank_5p", "max_flank_3p", "n_intermediates"]
        fh.write("\t".join(cols) + "\n")
        for r in loci:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    with (out / "intermediates.tsv").open("w", encoding="utf-8") as fh:
        cols = ["hairpin_id", "hp_start", "hp_end", "strand", "transcript_id",
                "tx_start", "tx_end", "kind", "offset_to_boundary", "strict", "n_cluster"]
        fh.write("\t".join(cols) + "\n")
        for i in sorted(intermeds, key=lambda x: (x.hairpin.name, x.kind, x.endpoint)):
            fh.write("\t".join(str(v) for v in [
                i.hairpin.name, i.hairpin.start, i.hairpin.end, i.hairpin.strand,
                i.tx.tid, i.tx.start, i.tx.end, i.kind, i.offset,
                "yes" if i.strict else "no", i.n_cluster]) + "\n")

    (out / "anchor_stats.json").write_text(json.dumps(st, indent=2), encoding="utf-8")


# --------------------------------------------------------------------------
# stage entry
# --------------------------------------------------------------------------
def anchor_pri(pins: list[Hairpin], txs: dict[str, Transcript],
               min_flank: int = 0, allow_partial: bool = False,
               window: int = 50) -> tuple[list[AnchorPair], list[Intermediate], list[dict], dict]:
    pairs = find_pairs(pins, txs, min_flank, allow_partial)
    intermeds = find_intermediates(pins, txs, pairs, window)
    loci = locus_rows(pairs, intermeds)
    return pairs, intermeds, loci, stats(pairs, intermeds)


def catalog_sources(cat_dir: Path, cfg: RunConfig) -> list[str]:
    """Provenance of the transcript catalog.

    The catalog stage writes ``meta.json`` with the real sources; without it we
    fall back to whatever inputs were given on this run.
    """
    meta = cat_dir / "meta.json"
    if meta.exists():
        try:
            src = json.loads(meta.read_text(encoding="utf-8")).get("sources")
            if src:
                return list(src)
        except json.JSONDecodeError:
            pass
    if cfg.iso and cfg.rna:
        return ["Iso-Seq", "RNA-Seq"]
    return ["Iso-Seq"] if cfg.iso else ["RNA-Seq"]


def stage_anchor(cfg: RunConfig, rep) -> None:
    cat_dir = cfg.stage_dir("catalog")
    gtf = cat_dir / "transcripts.gtf"
    if not gtf.exists():
        sys.stderr.write(f"[anchor] catalog not found: {gtf}\n")
        raise SystemExit(4)
    pins = load_hairpins_bed(cfg.mirna_bed)
    # Anchor on the *filtered* union catalog. The per-source files next to it are
    # pre-filter copies kept for provenance and must not be used here, otherwise
    # models removed by the catalog filters would still be reported as evidence.
    sources = catalog_sources(cat_dir, cfg)
    txs = parse_gtf(gtf, source="+".join(sources))
    pairs, intermeds, loci, st = anchor_pri(
        pins, txs, min_flank=0, allow_partial=cfg.allow_partial,
        window=cfg.extend_hairpin or 50)
    write_outputs(cfg.stage_dir("anchor"), pairs, intermeds, loci, st)
    rep.facts["anchor"] = st
    _status(f"{st['n_loci']} loci / {st['n_pairs']} pairs "
            f"(spliced {st['spliced']}, intronic {st['intronic']}), "
            f"{st['n_intermediates']} intermediates")
