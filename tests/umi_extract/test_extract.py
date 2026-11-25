"""Run with: python3 -m unittest discover -s tests/umi_extract -v"""

import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


MODULE = Path(__file__).resolve().parents[2] / "modules/local/umitools.nf"
RECORDS = [("read1", "ACGT", "ATGCTAGCTAGCTAGCTAGC", "GA"),
           ("read2", "TTAA", "GCTAGCTAGCTAGCTAGCTAGCT", "CT")]


@unittest.skipUnless(shutil.which("nextflow") and shutil.which("umi_tools"),
                     "Nextflow and UMI-tools are required")
class UmiExtractTest(unittest.TestCase):
    def check_extract(self, pattern, method="regex"):
        with tempfile.TemporaryDirectory(prefix="umi-extract-") as directory:
            root = Path(directory)
            reads = root / "reads with space.fastq"
            reads.write_text("".join("@" + name + "\n" + left + insert + right + "\n+\n" +
                                     "I" * len(left + insert + right) + "\n"
                                     for name, left, insert, right in RECORDS))
            (root / "main.nf").write_text(
                "nextflow.enable.dsl=2\n"
                "include { UMITOOLS_EXTRACT } from '" + MODULE.as_posix() + "'\n"
                "workflow {\n"
                "    UMITOOLS_EXTRACT(Channel.of(tuple('sample', file(params.fastq))))\n"
                "}\n")
            (root / "nextflow.config").write_text(
                "process.executor = 'local'\ndocker.enabled = false\nsingularity.enabled = false\n")
            (root / "params.json").write_text(json.dumps({
                "fastq": str(reads), "umi_pattern": pattern, "umi_extract_method": method}))
            result = subprocess.run(
                ["nextflow", "run", "main.nf", "-ansi-log", "false", "-params-file", "params.json"],
                cwd=root, env=dict(os.environ, NXF_DISABLE_CHECK_LATEST="true",
                                   UMI_QUOTING_SENTINEL="unexpected-expansion"),
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            self.assertEqual(result.returncode, 0, result.stdout)
            outputs = list((root / "work").glob("*/*/sample.umi_extract.fastq.gz"))
            self.assertEqual(len(outputs), 1)
            with gzip.open(outputs[0], "rt") as handle:
                lines = handle.read().splitlines()
            expected = []
            for name, left, insert, right in RECORDS:
                umi = left + right if method == "regex" else left
                sequence = insert if method == "regex" else insert + right
                expected.extend(["@" + name + "_" + umi, sequence, "+", "I" * len(sequence)])
            self.assertEqual(lines, expected)
            self.assertIn(pattern, outputs[0].with_name("sample.umi_extract.log").read_text())
            self.assertTrue(outputs[0].with_name("versions.yml").is_file())

    def test_string_pattern(self):
        self.check_extract("NNNN", "string")

    def test_split_regex(self):
        self.check_extract("^(?P<umi_1>.{4}).+(?P<umi_2>.{2})$")

    def test_regex_with_apostrophe(self):
        self.check_extract("^(?P<umi_1>.{4})(?#5' UMI).+(?P<umi_2>.{2})$")

    def test_regex_with_repeated_apostrophes(self):
        self.check_extract("^(?P<umi_1>.{4})(?#5'' and '3' UMI).+(?P<umi_2>.{2})$")

    def test_regex_with_shell_characters(self):
        comment = "don't expand \"$UMI_QUOTING_SENTINEL\" or `printf ALTERED`; & | * ? [ACGT] {1,2}"
        self.check_extract("^(?P<umi_1>.{4})(?#" + comment + ").+(?P<umi_2>.{2})$")

    def test_regex_with_backslashes(self):
        self.check_extract(r"\A(?P<umi_1>[ACGT]{4})[ACGT]+(?P<umi_2>[ACGT]{2})\Z")

    def test_regex_with_whitespace(self):
        self.check_extract("^(?P<umi_1>.{4})(?#5'\tUMI\nand 3' UMI).+(?P<umi_2>.{2})$")


if __name__ == "__main__":
    unittest.main()
