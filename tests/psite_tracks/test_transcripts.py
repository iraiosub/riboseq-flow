"""Transcript-reference selection; run with the other tests in tests/psite_tracks."""
import csv
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[2]
RSCRIPT = os.environ.get("RSCRIPT", "Rscript")


class TranscriptSelectionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.tmp.name)
        rows = []
        models = [
            ("g1", "txA", "+", "protein_coding", "basic", [(101, 160), (201, 260)], [(119, 160), (201, 218)]),
            ("g1", "txB", "+", "protein_coding", "basic", [(101, 160), (301, 330)], [(119, 160), (301, 318)]),
            ("g1", "txNMD", "+", "nonsense_mediated_decay", "basic", [(101, 160), (401, 460)], [(119, 160), (401, 430)]),
            ("g2", "txPartial", "+", "protein_coding", "cds_start_NF", [(501, 620)], [(501, 560)]),
            ("g3", "txNoncoding", "+", "lncRNA", "basic", [(701, 790)], []),
            ("g4", "txMinus", "-", "protein_coding", "basic", [(901, 960), (1001, 1060)], [(943, 960), (1001, 1042)]),
        ]
        for gene, transcript, strand, biotype, tag, exons, cds in models:
            attributes = 'gene_id "{}"; transcript_id "{}"; gene_name "{}"; gene_type "protein_coding"; transcript_type "{}"; tag "{}";'.format(gene, transcript, gene, biotype, tag)
            for feature, ranges in (("exon", exons), ("CDS", cds)):
                for start, end in ranges:
                    rows.append('\t'.join(map(str, ("chrT", "test", feature, start, end, ".", strand, "0" if feature == "CDS" else ".", attributes))))
        cls.gtf = cls.root / "genes.gtf"
        cls.gtf.write_text('\n'.join(rows) + '\n')
        cls.tables = {}
        for mode in ("default", "representative", "all_cds"):
            work = cls.root / mode
            work.mkdir()
            command = [RSCRIPT, str(REPO / "bin/get_transcript_info.R"), "-g", str(cls.gtf), "-o", "Homo sapiens"]
            if mode != "default":
                command += ["--transcripts", mode]
            result = subprocess.run(command, cwd=work, capture_output=True, text=True)
            if result.returncode:
                raise AssertionError(result.stdout + result.stderr)
            cls.tables[mode] = next(work.glob("*.transcript_info.tsv"))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def read_table(self, mode):
        with self.tables[mode].open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_default_preserves_representatives_and_coordinates(self):
        self.assertEqual(self.tables["default"].read_bytes(), self.tables["representative"].read_bytes())
        rows = self.read_table("default")
        self.assertEqual({row["transcript_id"] for row in rows}, {"txA", "txMinus"})
        for row in rows:
            self.assertEqual((row["cds_start"], row["cds_end"], row["cds_length"], row["tx_length"]), ("19", "78", "60", "120"))
        self.assertEqual(self.tables["default"].name, "genes.longest_cds.transcript_info.tsv")

    def test_all_cds_includes_alternative_nmd_and_partial_cds(self):
        rows = self.read_table("all_cds")
        self.assertEqual({row["transcript_id"] for row in rows}, {"txA", "txB", "txNMD", "txPartial", "txMinus"})
        self.assertEqual(len(rows), 5)
        self.assertEqual(self.tables["all_cds"].name, "genes.all_cds.transcript_info.tsv")
        self.assertEqual(sum(row["gene_id"] == "g1" for row in rows), 3)
        self.assertEqual(next(row for row in rows if row["transcript_id"] == "txNMD")["cds_length"], "72")
        gtf = (self.root / "all_cds/genes.all_cds_transcripts.gtf").read_text()
        self.assertIn('transcript_id "txNMD"', gtf)
        self.assertNotIn('transcript_id "txNoncoding"', gtf)

    def test_all_cds_does_not_require_biotype_or_gene_name_attributes(self):
        work = self.root / "bare"
        work.mkdir()
        attributes = 'gene_id "gene"; transcript_id "txBare";'
        gtf = work / "bare.gtf"
        gtf.write_text('chrT\ttest\texon\t1\t120\t.\t+\t.\t' + attributes + '\n'
                       'chrT\ttest\tCDS\t19\t78\t.\t+\t0\t' + attributes + '\n')
        result = subprocess.run([RSCRIPT, str(REPO / "bin/get_transcript_info.R"), "-g", str(gtf),
                                 "-o", "Homo sapiens", "--transcripts", "all_cds"], cwd=work, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        with (work / "bare.all_cds.transcript_info.tsv").open() as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(row["gene_name"], "gene")
        self.assertEqual(row["transcript_id"], "txBare")


if __name__ == "__main__":
    unittest.main()
