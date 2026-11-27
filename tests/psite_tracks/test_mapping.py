"""Run with: python3 -m unittest discover -s tests/psite_tracks"""
import csv
import gzip
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[2]
RSCRIPT = os.environ.get("RSCRIPT", "Rscript")


class PsiteTracksTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        gtf = []
        for transcript, strand, exons in [
            ("txP", "+", [(101, 110), (201, 210)]),
            ("txPcopy", "+", [(101, 110), (201, 210)]),
            ("txPalt", "+", [(51, 55), (101, 110), (201, 210)]),
            ("txM", "-", [(301, 310), (401, 410)]),
            ("txMalt", "-", [(301, 310), (401, 410), (451, 455)]),
        ]:
            for rank, (start, end) in enumerate(exons, 1):
                gtf.append('chrT\ttest\texon\t{}\t{}\t.\t{}\t.\tgene_id "{}"; transcript_id "{}"; exon_number "{}";'.format(
                    start, end, strand, transcript, transcript, rank))
        (self.root / "genes.gtf").write_text("\n".join(gtf) + "\n")
        (self.root / "genome.fai").write_text("chrT\t1000\t0\t60\t61\n")
        (self.root / "transcripts.tsv").write_text("transcript_id\ntxP\ntxPcopy\ntxPalt\ntxM\ntxMalt\n")
        self.run_r("prepare_psite_tracks.R", "-g", "genes.gtf", "-f", "genome.fai", "-t", "transcripts.tsv")

    def run_r(self, script, *args):
        result = subprocess.run([RSCRIPT, str(REPO / "bin" / script)] + list(args),
                                cwd=str(self.root), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        self.assertEqual(result.returncode, 0, result.stdout)

    def map_sites(self, rows):
        with gzip.open(self.root / "sample.psite.tsv.gz", "wt") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["read_name", "transcript", "psite"])
            writer.writerows(rows)
        self.run_r("get_psite_tracks.R", "-p", "sample.psite.tsv.gz", "-a", "psite_exons.rds",
                   "-s", "genome.chrom.sizes", "-o", "sample")

    def track(self, strand):
        with (self.root / ("sample.psites." + strand + ".unsorted.bedgraph")).open() as handle:
            return {tuple(row[:3]): int(row[3]) for row in csv.reader(handle, delimiter="\t")}

    def test_read_identity_strands_and_exon_boundaries(self):
        self.map_sites([
            ("a", "txP", 1), ("b", "txP", 1), ("a", "txPcopy", 1),
            ("c", "txP", 10), ("d", "txP", 11), ("e", "txP", 20),
            ("m1", "txM", 1), ("m10", "txM", 10), ("m11", "txM", 11), ("m20", "txM", 20),
            ("conflict", "txP", 1), ("conflict", "txP", 11), ("missing", "unknown", 5),
            ("0001", "txP", 1), ("m20", "txM", 20),
        ])
        self.assertEqual(self.track("forward"), {("chrT", "100", "101"): 3,
            ("chrT", "109", "110"): 1, ("chrT", "200", "201"): 1, ("chrT", "209", "210"): 1})
        self.assertEqual(self.track("reverse"), {("chrT", "409", "410"): 1,
            ("chrT", "400", "401"): 1, ("chrT", "309", "310"): 1, ("chrT", "300", "301"): 1})
        with (self.root / "sample.psites.mapping_stats.tsv").open() as handle:
            stats = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({key: int(value) for key, value in stats.items() if key != "sample"}, {
            "input_records": 15, "input_read_names": 12, "mapped_records": 14, "unmapped_records": 1,
            "duplicate_assignments": 2, "ambiguous_reads": 1, "forward_psites": 6, "reverse_psites": 4})
        with gzip.open(self.root / "sample.codon_coverage_psite.bed.gz", "rt") as handle:
            self.assertEqual(sum(int(row[4]) for row in csv.reader(handle, delimiter="\t")), 10)

    def test_high_coverage_retains_distinct_reads(self):
        self.map_sites([("read_" + str(i), "txP", 1) for i in range(2000)])
        self.assertEqual(self.track("forward"), {("chrT", "100", "101"): 2000})
        self.assertFalse((self.root / "sample.psites.reverse.unsorted.bedgraph").exists())

    def test_all_ambiguous_or_unmapped(self):
        self.map_sites([("conflict", "txP", 1), ("conflict", "txM", 1), ("missing", "unknown", 1)])
        self.assertFalse(list(self.root.glob("*.unsorted.bedgraph")))
        with gzip.open(self.root / "sample.codon_coverage_psite.bed.gz", "rt") as handle:
            self.assertEqual(handle.read(), "")

    def test_isoforms_collapse_after_projection_and_bigwig_round_trip(self):
        self.map_sites([
            ("a", "txP", 1), ("a", "txPalt", 6), ("b", "txPcopy", 1),
            ("junction", "txP", 11), ("junction", "txPalt", 16),
            ("m", "txM", 1), ("m", "txMalt", 6),
            ("mj", "txM", 11), ("mj", "txMalt", 16),
            ("ambiguous", "txP", 1), ("ambiguous", "txPalt", 16),
        ])
        self.assertEqual(self.track("forward"), {("chrT", "100", "101"): 2, ("chrT", "200", "201"): 1})
        self.assertEqual(self.track("reverse"), {("chrT", "409", "410"): 1, ("chrT", "309", "310"): 1})
        with (self.root / "sample.psites.mapping_stats.tsv").open() as handle:
            stats = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(int(stats["duplicate_assignments"]), 4)
        self.assertEqual(int(stats["ambiguous_reads"]), 1)
        self.assertEqual(int(stats["forward_psites"]) + int(stats["reverse_psites"]), 5)
        code = '''
        library(rtracklayer)
        for (direction in c("forward", "reverse")) {
          path <- paste0("sample.psites.", direction, ".unsorted.bedgraph")
          expected <- sort(import(path, format = "bedGraph"))
          seqlengths(expected) <- c(chrT = 1000)
          export(expected, paste0(direction, ".bw"), format = "BigWig")
          actual <- import(paste0(direction, ".bw"))
          stopifnot(identical(start(expected), start(actual)), identical(end(expected), end(actual)),
                    identical(as.numeric(score(expected)), as.numeric(score(actual))))
        }
        '''
        result = subprocess.run([RSCRIPT, "-e", code], cwd=self.root, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
