"""Unit tests that need no external tools (run with `python -m unittest`)."""

from __future__ import annotations

import shutil
import unittest
from pathlib import Path

from pripipe.cli import main
from pripipe.config import RunConfig
from pripipe.models import Hairpin
from pripipe.stages.anchor import anchor_pri, parse_gtf
from pripipe.stages.catalog import apply_filters, dedup_structural
from pripipe.models import Transcript
from pripipe.validate import Report, read_hairpins

SCRATCH = Path(__file__).resolve().parent / "_scratch"


class ScratchTest(unittest.TestCase):
    """Each test gets its own directory under pripipe/tests/_scratch."""

    def setUp(self) -> None:
        self.tmp = SCRATCH / self._testMethodName
        shutil.rmtree(self.tmp, ignore_errors=True)
        self.tmp.mkdir(parents=True, exist_ok=True)

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp, ignore_errors=True)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def gtf_line(chrom, feature, start, end, strand, tid, exon_number=None):
    attrs = f'gene_id "{tid}"; transcript_id "{tid}";'
    if exon_number is not None:
        attrs += f' exon_number "{exon_number}";'
    return f"{chrom}\tsrc\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n"


class AnchorTest(ScratchTest):
    def setUp(self) -> None:
        super().setUp()
        self.gtf = self.tmp / "cat.gtf"
        body = ""
        # T1: single exon 0-500  -> hp1 (100-200) is spliced
        body += gtf_line("chr1", "transcript", 1, 500, "+", "T1")
        body += gtf_line("chr1", "exon", 1, 500, "+", "T1", 1)
        # T2: exons 0-150 + 250-500 -> hp1 spans a junction (intronic), hp2 (300-350) spliced
        body += gtf_line("chr1", "transcript", 1, 500, "+", "T2")
        body += gtf_line("chr1", "exon", 1, 150, "+", "T2", 1)
        body += gtf_line("chr1", "exon", 250, 500, "+", "T2", 2)
        # T3a/T3b: start inside hp1 -> shared 3' remnant (intermediate)
        body += gtf_line("chr1", "transcript", 150, 400, "+", "T3a")
        body += gtf_line("chr1", "exon", 150, 400, "+", "T3a", 1)
        body += gtf_line("chr1", "transcript", 150, 420, "+", "T3b")
        body += gtf_line("chr1", "exon", 150, 420, "+", "T3b", 1)
        write(self.gtf, body)
        self.hp1 = Hairpin("chr1", 100, 200, "hp1", "0", "+")
        self.hp2 = Hairpin("chr1", 300, 350, "hp2", "0", "+")

    def test_parse_gtf(self) -> None:
        txs = parse_gtf(self.gtf)
        self.assertEqual(len(txs), 4)
        self.assertEqual(txs["T2"].n_exons, 2)
        self.assertEqual((txs["T1"].start, txs["T1"].end), (0, 500))

    def test_anchor_types(self) -> None:
        pairs, _, loci, st = anchor_pri([self.hp1, self.hp2], parse_gtf(self.gtf),
                                        min_flank=0, allow_partial=False, window=50)
        got = {(p.hairpin.name, p.tx.tid): p.anchor_type for p in pairs}
        # hp1 (100-200): T1 spliced, T2 spans the junction -> intronic
        self.assertEqual(got[("hp1", "T1")], "spliced")
        self.assertEqual(got[("hp1", "T2")], "intronic")
        # hp2 (300-350): contained by all four transcripts
        self.assertEqual(got[("hp2", "T1")], "spliced")
        self.assertEqual(got[("hp2", "T3a")], "spliced")
        self.assertEqual(st["n_loci"], 2)
        self.assertEqual(st["n_pairs"], 6)

    def test_intermediate_cluster(self) -> None:
        _, intermeds, _, st = anchor_pri([self.hp1], parse_gtf(self.gtf),
                                         min_flank=0, allow_partial=False, window=50)
        kinds = {(i.tx.tid, i.kind) for i in intermeds}
        self.assertIn(("T3a", "3prime_remnant"), kinds)
        self.assertIn(("T3b", "3prime_remnant"), kinds)
        self.assertEqual(st["n_intermediates"], 2)
        # GTF is 1-based: transcript start 150 -> 0-based 149, i.e. 49 bp inside hp1
        offsets = {i.offset for i in intermeds}
        self.assertEqual(len(offsets), 1)
        self.assertTrue(49 <= offsets.pop() <= 50)

    def test_no_partial_by_default(self) -> None:
        pairs, _, _, _ = anchor_pri([self.hp1], parse_gtf(self.gtf),
                                    min_flank=0, allow_partial=False, window=50)
        self.assertNotIn("T3a", {p.tx.tid for p in pairs})


class ValidateTest(ScratchTest):
    def setUp(self) -> None:
        super().setUp()
        write(self.tmp / "g.fa", ">chr1\n" + "A" * 1000 + "\n")
        write(self.tmp / "g.fa.fai", "chr1\t1000\t6\t1000\t1001\n")
        stub = self.tmp / "sq3"
        stub.mkdir()
        write(stub / "sqanti3_qc.py", "# stub\n")
        self.sq3 = stub

    def _cfg(self, bed: Path) -> RunConfig:
        return RunConfig(genome=self.tmp / "g.fa", mirna_bed=bed, sqanti3_dir=self.sq3)

    def _pins(self, bed: Path):
        rep = Report()
        pins = read_hairpins(self._cfg(bed), {"chr1": 1000}, rep)
        return pins, rep

    def test_good_bed(self) -> None:
        bed = write(self.tmp / "ok.bed", "chr1\t100\t200\thp1\t0\t+\nchr1\t300\t400\thp2\t0\t-\n")
        pins, rep = self._pins(bed)
        self.assertEqual(len(pins), 2)
        self.assertEqual(rep.errors, [])

    def test_bad_strand(self) -> None:
        bed = write(self.tmp / "bad.bed", "chr1\t100\t200\thp1\t0\t.\n")
        _, rep = self._pins(bed)
        self.assertTrue(any("strand" in e for e in rep.errors))

    def test_duplicate_name(self) -> None:
        bed = write(self.tmp / "dup.bed",
                    "chr1\t100\t200\thp1\t0\t+\nchr1\t300\t400\thp1\t0\t+\n")
        _, rep = self._pins(bed)
        self.assertTrue(any("duplicate" in e for e in rep.errors))

    def test_out_of_bounds(self) -> None:
        bed = write(self.tmp / "oob.bed", "chr1\t900\t1200\thp1\t0\t+\n")
        _, rep = self._pins(bed)
        self.assertTrue(any("exceeds" in e for e in rep.errors))


class DedupTest(unittest.TestCase):
    """Redundancy removal must never drop a transcript without recording it."""

    def _txs(self):
        return {
            f"T{i}": Transcript(tid=f"T{i}", chrom="chr1", strand="+",
                                exons=[(i * 20, 500 + i * 20)])
            for i in range(6)   # TSS 0,20,40,...,PAS 500,520,...
        }

    def test_no_silent_loss(self) -> None:
        txs = self._txs()
        kept, collapsed = dedup_structural(txs, tss_tol=50, pas_tol=50)
        self.assertEqual(len(kept) + len(collapsed), len(txs))
        self.assertLessEqual(len(kept), len(txs))
        for member, rep in collapsed:
            self.assertIn(rep, kept)
            self.assertNotIn(member, kept)

    def test_tolerance_splits_clusters(self) -> None:
        txs = self._txs()
        kept_tight, collapsed_tight = dedup_structural(txs, tss_tol=5, pas_tol=5)
        self.assertEqual(len(collapsed_tight), 0)
        self.assertEqual(len(kept_tight), len(txs))

    def test_longest_wins(self) -> None:
        txs = {
            "short": Transcript("short", "chr1", "+", [(0, 500)]),
            "long": Transcript("long", "chr1", "+", [(0, 530)]),
        }
        kept, collapsed = dedup_structural(txs, 50, 50)
        self.assertIn("long", kept)
        self.assertEqual(collapsed, [("short", "long")])

    def test_filter_reasons(self) -> None:
        from pripipe.config import RunConfig
        cfg = RunConfig(genome=Path("g.fa"), mirna_bed=Path("b.bed"))
        txs = {
            "tiny": Transcript("tiny", "chr1", "+", [(0, 100)]),
            "huge": Transcript("huge", "chr1", "+", [(0, 200_000)]),
            "ok": Transcript("ok", "chr1", "+", [(0, 400), (600, 900)]),
        }
        kept, dropped, flagged, stats = apply_filters(txs, {}, cfg)
        reasons = {tid: reason for tid, reason, _ in dropped}
        self.assertEqual(reasons["tiny"], "shorter_than_min")
        self.assertEqual(reasons["huge"], "longer_than_max")
        self.assertIn("ok", kept)
        self.assertEqual(stats["n_in"], 3)


class CliTest(unittest.TestCase):
    def test_version(self) -> None:
        with self.assertRaises(SystemExit) as ctx:
            main(["--version"])
        self.assertEqual(ctx.exception.code, 0)

    def test_help(self) -> None:
        with self.assertRaises(SystemExit) as ctx:
            main(["--help"])
        self.assertEqual(ctx.exception.code, 0)


if __name__ == "__main__":
    unittest.main()
