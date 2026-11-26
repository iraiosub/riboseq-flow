"""Run with: python3 -m unittest discover -s tests/rust_qc"""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "rust_multiqc.py"


class RustQCTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.write_table("qc.tsv", ["name", "most_abundant_length", "most_periodic_length"],
                         [("A", 28, 29), ("B", 31, 27)])
        self.write_table("counts.tsv", ["sample", "length", "psites_used"],
                         [("A", 27, 1000), ("A", 28, 999), ("A", 29, 1500),
                          ("B", 27, 2000), ("B", 31, 3000)])
        self.write_table("profiles.tsv", ["sample", "length", "offset_codon", "kld"],
                         [(sample, length, offset, length / 1000 + offset / 100000)
                          for sample, length in (("A", 27), ("A", 28), ("A", 29), ("B", 27), ("B", 31))
                          for offset in range(10, -11, -1)])

    def write_table(self, name, header, rows):
        with (self.root / name).open("w") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(header)
            writer.writerows(rows)

    def run_summary(self):
        subprocess.run([sys.executable, str(SCRIPT), "--qc-summary", "qc.tsv",
                        "--profiles", "profiles.tsv", "--lengths", "counts.tsv"],
                       cwd=self.root, check=True)

    def read_section(self, name):
        return json.loads((self.root / ("rust_" + name + "_mqc.json")).read_text())

    def test_support_boundary_and_lossless_heatmap(self):
        self.run_summary()
        section = self.read_section("all_lengths")
        self.assertEqual(section["ycats"], ["S1 | 27 nt | A", "S1 | 29 nt | A", "S2 | 27 nt | B", "S2 | 31 nt | B"])
        self.assertEqual(section["xcats"], list(map(str, range(-10, 11))))
        self.assertEqual(len(section["data"]), 4)
        self.assertEqual(len(section["data"][0]), 21)
        self.assertAlmostEqual(section["data"][0][0], 0.0269)
        self.assertAlmostEqual(section["data"][0][-1], 0.0271)
        self.assertNotIn("extra", section)

    def test_qc_selections_are_preserved_without_substitution(self):
        self.run_summary()
        frequent = self.read_section("frequent")
        self.assertEqual(set(frequent["data"]), {"B | 31 nt"})
        self.assertIn("A: 28 nt", frequent["description"])
        periodic = self.read_section("periodic")
        self.assertEqual(set(periodic["data"]), {"A | 29 nt", "B | 27 nt"})
        self.assertAlmostEqual(periodic["data"]["A | 29 nt"]["0"], 0.029)

    def test_no_supported_lengths_keeps_explanation_and_publication(self):
        self.write_table("counts.tsv", ["sample", "length", "psites_used"], [("A", 27, 999)])
        self.run_summary()
        section = self.read_section("support")
        self.assertEqual(section["plot_type"], "html")
        self.assertIn("1,000", section["data"])
        self.assertIn("10.1038/ncomms12915", section["parent_description"])
        self.assertFalse((self.root / "rust_all_lengths_mqc.json").exists())

    def test_missing_periodic_selection_has_no_invented_curve(self):
        self.write_table("qc.tsv", ["name", "most_abundant_length", "most_periodic_length"],
                         [("B", 31, "NA")])
        self.run_summary()
        section = self.read_section("periodic")
        self.assertEqual(section["plot_type"], "html")
        self.assertIn("B: no selection", section["data"])


if __name__ == "__main__":
    unittest.main()
