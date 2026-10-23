"""Real featureCounts regressions for combined and empty sample columns."""

import csv
import importlib.util
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
MISSING = shutil.which("featureCounts") is None or importlib.util.find_spec("pysam") is None


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t"))


@unittest.skipIf(MISSING, "Requires featureCounts and pysam")
class CountFeaturesTests(unittest.TestCase):
    def make_fixture(self, directory):
        annotation = directory / "reference.gff3"
        annotation.write_text(
            "##gff-version 3\n"
            "chr\ttest\tCDS\t200\t300\t.\t+\t0\tgene_id=cds_gene\n"
            "chr\ttest\trRNA\t600\t700\t.\t+\t.\tgene_id=rrna_gene\n")
        header = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr\tLN:1000\n"
        sequence, quality = "ACGT" * 7 + "AC", "I" * 30

        def mapped(name, position):
            return f"{name}\t0\tchr\t{position}\t60\t30M\t*\t0\t0\t{sequence}\t{quality}\n"

        empty, first, second = (directory / name for name in ("empty.sam", "first sample.sam", "second.sam"))
        empty.write_text(header)
        first.write_text(header + mapped("cds", 201) + mapped("rrna", 601))
        second.write_text(header + mapped("cds", 250) +
                          f"unmapped\t4\t*\t0\t0\t*\t*\t0\t0\t{sequence}\t{quality}\n")
        return annotation, [empty, first, second]

    def run_helper(self, directory, feature_type):
        annotation, sams = self.make_fixture(directory)
        output = directory / "combined.txt"
        result = subprocess.run([sys.executable, str(REPO / "count_features.py"), "--threads", "1",
                                 "--annotation", str(annotation), "--gene-id", "gene_id",
                                 "--feature-type", feature_type, "--allow-overlap", "--output", str(output),
                                 *map(str, sams)], text=True, capture_output=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertFalse(list(directory.glob(".featurecounts-*")))
        return sams, rows(output), rows(Path(str(output) + ".summary"))

    def test_independent_samples_preserve_empty_and_nonzero_columns(self):
        with tempfile.TemporaryDirectory() as temp:
            sams, table, summary = self.run_helper(Path(temp), "CDS")
            self.assertEqual(table[0], ["Geneid", "Chr", "Start", "End", "Strand", "Length", *map(str, sams)])
            self.assertEqual(table[1], ["cds_gene", "chr", "200", "300", "+", "101", "0", "1", "1"])
            self.assertEqual(summary[0], ["Status", *map(str, sams)])
            counts = {row[0]: list(map(int, row[1:])) for row in summary[1:]}
            self.assertEqual(counts["Assigned"], [0, 1, 1])
            self.assertEqual(counts["Unassigned_NoFeatures"], [0, 1, 0])
            self.assertEqual(counts["Unassigned_Unmapped"], [0, 0, 1])

    def test_absent_feature_type_writes_zero_matrix_and_balanced_summary(self):
        with tempfile.TemporaryDirectory() as temp:
            sams, table, summary = self.run_helper(Path(temp), "tRNA")
            self.assertEqual(len(table), 1)
            self.assertEqual(table[0][6:], list(map(str, sams)))
            counts = {row[0]: list(map(int, row[1:])) for row in summary[1:]}
            self.assertEqual(counts["Assigned"], [0, 0, 0])
            self.assertEqual(counts["Unassigned_NoFeatures"], [0, 2, 1])
            self.assertEqual(counts["Unassigned_Unmapped"], [0, 0, 1])
            self.assertEqual([sum(values[i] for values in counts.values()) for i in range(3)], [0, 2, 2])


if __name__ == "__main__":
    unittest.main()
