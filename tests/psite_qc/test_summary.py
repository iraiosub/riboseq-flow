"""Run with: python3 -m unittest discover -s tests/psite_qc"""
import csv
import gzip
from pathlib import Path
import json
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "psite_multiqc.py"


class PsiteQCTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.a = self.write_table("a.psite.tsv.gz", "sampleA", [
            (27, -26), (27, -25), (27, -1), (31, 0), (31, 0),
            (31, 1), (31, 2), (31, 40), (31, 41), (31, 60),
            (31, 100), (31, 101), (31, 125), (31, 126),
        ])
        self.b = self.write_table("b.psite.tsv", "sample<B>&", [(33, 200)])

    def write_table(self, filename, sample, sites):
        path = self.root / filename
        opener = gzip.open if filename.endswith(".gz") else open
        with opener(path, "wt") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "length", "psite_from_start", "psite_from_stop", "psite_region"])
            for length, start in sites:
                region = "5utr" if start < 0 else "cds" if start <= 100 else "3utr"
                writer.writerow([sample, length, start, start - 100, region])
        return path

    def run_summary(self, *paths):
        subprocess.run([sys.executable, str(SCRIPT)] + [str(path) for path in paths],
                       cwd=str(self.root), check=True)

    def read_plot(self, name):
        return json.loads((self.root / ("psite_" + name + "_mqc.json")).read_text())

    def test_region_and_frame_counts_and_retained_lengths(self):
        self.run_summary(self.b, self.a)
        regions = self.read_plot("regions")
        self.assertEqual(regions["data"]["sampleA"], {"5' UTR": 3, "CDS": 8, "3' UTR": 3})
        self.assertEqual(regions["data"]["sample<B>&"], {"5' UTR": 0, "CDS": 0, "3' UTR": 1})
        frames = self.read_plot("frames")["data"]
        self.assertEqual(frames["sampleA | 5' UTR"], {"Frame 0": 0, "Frame 1": 1, "Frame 2": 2})
        self.assertEqual(frames["sampleA | CDS"], {"Frame 0": 3, "Frame 1": 3, "Frame 2": 2})
        self.assertEqual(frames["sampleA | 3' UTR"], {"Frame 0": 1, "Frame 1": 0, "Frame 2": 2})
        self.assertNotIn("sample<B>& | CDS", frames)
        self.assertIn("<strong>sampleA</strong>: 27, 31 nt", regions["extra"])
        self.assertIn("<strong>sample&lt;B&gt;&amp;</strong>: 33 nt", regions["extra"])
        self.assertEqual(regions["parent_name"], "riboWaltz: P-site QC")
        self.assertIn("riboWaltz", regions["section_href"])
        self.assertEqual(regions["categories"], {"5' UTR": {}, "CDS": {}, "3' UTR": {}})
        self.assertEqual(self.read_plot("frames")["categories"],
                         {"Frame 0": {}, "Frame 1": {}, "Frame 2": {}})

    def test_profile_boundaries_zero_filling_and_normalisation(self):
        self.run_summary(self.a, self.b)
        start = self.read_plot("start")["data"]
        stop = self.read_plot("stop")["data"]
        self.assertEqual(set(start), {"sampleA"})
        self.assertEqual(set(stop), {"sampleA"})
        start = {int(position): value for position, value in start["sampleA"].items()}
        stop = {int(position): value for position, value in stop["sampleA"].items()}
        self.assertEqual(set(start), set(range(-25, 41)))
        self.assertEqual(set(stop), set(range(-40, 26)))
        self.assertAlmostEqual(sum(start.values()), 1)
        self.assertAlmostEqual(sum(stop.values()), 1)
        self.assertAlmostEqual(start[0], 2 / 7)
        self.assertAlmostEqual(start[-25], 1 / 7)
        self.assertAlmostEqual(start[40], 1 / 7)
        self.assertEqual(start[4], 0)
        self.assertEqual(stop[-40], 0.25)
        self.assertEqual(stop[25], 0.25)

    def test_input_order_does_not_change_report_data(self):
        self.run_summary(self.a, self.b)
        first = {path.name: path.read_bytes() for path in self.root.glob("*_mqc.json")}
        self.run_summary(self.b, self.a)
        self.assertEqual(first, {path.name: path.read_bytes() for path in self.root.glob("*_mqc.json")})

    def test_no_p_sites_produces_no_sections(self):
        empty = self.write_table("empty.psite.tsv.gz", "empty", [])
        self.run_summary(empty)
        self.assertEqual(list(self.root.glob("*_mqc.json")), [])


if __name__ == "__main__":
    unittest.main()
