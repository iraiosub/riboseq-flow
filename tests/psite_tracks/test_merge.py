"""Run with: python3 -m unittest discover -s tests/psite_tracks -p test_merge.py"""
import csv
import gzip
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[2]
RSCRIPT = os.environ.get("RSCRIPT", "Rscript")
COVERAGE_FILES = {
    "codon_coverage_psite.tsv.gz": ["transcript", "start", "end", "from_cds_start", "from_cds_stop", "region"],
    "codon_coverage_rpf.tsv.gz": ["transcript", "start", "end", "from_cds_start", "from_cds_stop", "region"],
    "cds_coverage_psite.tsv.gz": ["transcript", "length_cds", "length_selection"],
    "cds_plus42nt_minus27nt_coverage_psite.tsv.gz": ["transcript", "length_cds", "length_selection"],
}


class MergePsiteTablesTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="psite merge 'quote' ")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.tables = self.root / "tables"

    def write_table(self, path, header, rows):
        path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(path, "wt") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(header)
            writer.writerows(rows)

    def write_sample(self, directory, sample, rows):
        directory = self.tables / directory
        self.write_table(directory / "psite_offset.tsv.gz",
                         ["sample", "length", "corrected_offset_from_5"],
                         [(sample, 29, 12), (sample, 28, 12)])
        for filename, keys in COVERAGE_FILES.items():
            values = [3, 5, 0, -8, "cds"] if "start" in keys else [90, "all"]
            self.write_table(directory / filename, keys + [sample],
                             [(transcript, *values, count) for transcript, count in rows])

    def run_merge(self):
        result = subprocess.run([RSCRIPT, "--vanilla", str(REPO / "bin/merge_psite_tables.R"), str(self.tables)],
                                cwd=self.root, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertEqual({path.name for path in self.root.glob("*.tsv.gz")},
                         set(COVERAGE_FILES) | {"psite_offset.tsv.gz"})

    def read_table(self, filename):
        with gzip.open(self.root / filename, "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            return reader.fieldnames, list(reader)

    def test_merge_gzipped_samples_preserves_keys_counts_and_order(self):
        self.write_sample("01", "002", [("0001", 7), ("0003", 11)])
        self.write_sample("02", "001", [("0001", 5), ("0002", 3)])
        self.run_merge()
        header, offsets = self.read_table("psite_offset.tsv.gz")
        self.assertEqual([(row["sample"], row["length"]) for row in offsets],
                         [("001", "28"), ("001", "29"), ("002", "28"), ("002", "29")])
        self.assertTrue(all(row["corrected_offset_from_5"] == "12" for row in offsets))
        for filename, keys in COVERAGE_FILES.items():
            header, coverage = self.read_table(filename)
            self.assertEqual(header, keys + ["001", "002"])
            self.assertEqual({row["transcript"]: (row["001"], row["002"]) for row in coverage},
                             {"0001": ("5", "7"), "0002": ("3", "0"), "0003": ("0", "11")})
            values = ["3", "5", "0", "-8", "cds"] if "start" in keys else ["90", "all"]
            self.assertTrue(all([row[key] for key in keys[1:]] == values for row in coverage))

    def test_single_sample_reads_gzip_and_preserves_values(self):
        self.write_sample("01", "001", [("0001", 5)])
        self.run_merge()
        for filename, keys in COVERAGE_FILES.items():
            header, coverage = self.read_table(filename)
            self.assertEqual(header, keys + ["001"])
            self.assertEqual([(row["transcript"], row["001"]) for row in coverage], [("0001", "5")])


if __name__ == "__main__":
    unittest.main()
