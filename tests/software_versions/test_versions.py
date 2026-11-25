"""Run with: python3 -m unittest discover -s tests/software_versions"""

import importlib.util
from pathlib import Path
import subprocess
import sys
import shutil
import tempfile
import unittest

import yaml


SCRIPT = Path(__file__).resolve().parents[2] / "bin/collect_software_versions.py"
spec = importlib.util.spec_from_file_location("collect_versions", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class SoftwareVersionsTest(unittest.TestCase):
    @unittest.skipUnless(shutil.which("Rscript"), "Rscript is required to check R metadata")
    def test_r_runtime_and_package_versions(self):
        with tempfile.TemporaryDirectory() as directory:
            script = SCRIPT.with_name("software_versions.R")
            subprocess.run(["Rscript", str(script), "RIBOSEQ:PROBE", "stats"], cwd=directory, check=True)
            versions = yaml.safe_load((Path(directory) / "versions.yml").read_text())["RIBOSEQ:PROBE"]
            self.assertIsInstance(versions["R"], str)
            self.assertEqual(versions["stats"], versions["R"])

    def test_repeated_processes_keep_all_versions_without_duplicates(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = []
            for index, version in enumerate(("1.10", "1.18", "1.10")):
                path = Path(directory) / str(index) / "versions.yml"
                path.parent.mkdir()
                path.write_text('"RIBOSEQ:MAP":\n    samtools: ' + version + '\n    star: 2.7.10a\n')
                paths.append(path)
            premap = Path(directory) / "premap.yml"
            premap.write_text('"RIBOSEQ:PREMAP":\n    samtools: 1.10\n')
            paths.append(premap)
            versions = module.collect_versions(paths, "23.10.1", "1.3.0")
            self.assertEqual(versions["MAP"]["samtools"], ["1.10", "1.18"])
            self.assertEqual(versions["MAP"]["star"], "2.7.10a")
            self.assertEqual(versions["PREMAP"]["samtools"], "1.10")
            self.assertEqual(set(versions), {"MAP", "PREMAP", "Workflow"})
            self.assertEqual(versions["Workflow"]["Nextflow"], "23.10.1")
            self.assertEqual(versions["Workflow"]["riboseq-flow"], "1.3.0")

    def test_cli_writes_one_yaml_with_string_versions(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "versions.yml"
            source.write_text('"RIBOSEQ:FASTQC":\n    fastqc: 0.11.9\n'
                              '"RIBOSEQ:MAP":\n    samtools: 1.10\n'
                              '"RIBOSEQ:GENOMIC_PSITE_TRACKS:PSITE_BEDGRAPHTOBIGWIG":\n    bedGraphToBigWig: 469\n')
            subprocess.run([sys.executable, str(SCRIPT), "23.10.1", "1.3.0 / main / abc1234", str(source)],
                           cwd=root, check=True)
            versions = yaml.safe_load((root / "software_versions.yml").read_text())
            self.assertEqual(versions["MAP"]["samtools"], "1.10")
            self.assertEqual(versions["FASTQC"]["fastqc"], "0.11.9")
            self.assertEqual(versions["PSITE_BEDGRAPHTOBIGWIG"]["bedGraphToBigWig"], "469")
            self.assertEqual(versions["Workflow"]["riboseq-flow"], "1.3.0 / main / abc1234")
            self.assertEqual(sorted(path.name for path in root.iterdir()),
                             ["software_versions.yml", "versions.yml"])


if __name__ == "__main__":
    unittest.main()
