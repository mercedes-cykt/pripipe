"""Shared data model for pripipe."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Hairpin:
    """A pre-miRNA locus (BED6, 0-based half-open)."""

    chrom: str
    start: int
    end: int
    name: str
    score: str = "0"
    strand: str = "+"

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class Transcript:
    """A transcript model (exons are 0-based half-open, sorted)."""

    tid: str
    chrom: str
    strand: str
    exons: list[tuple[int, int]]
    source: str = "RNA-Seq"
    aliases: list[str] = field(default_factory=list)

    @property
    def start(self) -> int:
        return self.exons[0][0]

    @property
    def end(self) -> int:
        return self.exons[-1][1]

    @property
    def n_exons(self) -> int:
        return len(self.exons)

    def contains(self, start: int, end: int) -> bool:
        return self.start <= start and end <= self.end

    def exon_containing(self, start: int, end: int) -> int | None:
        """Index of the exon that fully contains [start, end), else None."""
        for i, (s, e) in enumerate(self.exons):
            if s <= start and end <= e:
                return i
        return None

    def chain_key(self) -> tuple:
        """Identity used to merge the same model coming from two catalogs."""
        return (self.chrom, self.strand, tuple(self.exons))


@dataclass
class AnchorPair:
    hairpin: Hairpin
    tx: Transcript
    anchor_type: str          # spliced | intronic | partial
    flank_5p: int
    flank_3p: int
    exon_index: int | None = None

    @property
    def flank_ok(self) -> bool:
        return min(self.flank_5p, self.flank_3p) >= 0


@dataclass
class Intermediate:
    hairpin: Hairpin
    tx: Transcript
    kind: str                 # 5prime_remnant | 3prime_remnant
    endpoint: int             # transcript end that sits at the cut site
    offset: int               # >0 inside hairpin, <0 outside
    strict: bool              # endpoint strictly inside the hairpin
    n_cluster: int = 1
