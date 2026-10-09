"""Regression tests for the post-trimming BBSplit read-length boundary."""

import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SPLITTER = REPO / "split_trimmed_fastq.py"


def fastq_record(name, length, *, plus="+", newline="\n", final_newline=True):
    sequence = ("ACGTN" * ((length + 4) // 5))[:length]
    qualities = "".join(chr(33 + position % 60) for position in range(length))
    lines = ["@" + name, sequence, plus, qualities]
    return (newline.join(lines) + (newline if final_newline else "")).encode("ascii")


def read_fastq_bytes(path):
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as handle:
            return handle.read()
    return path.read_bytes()


class LengthFilterTests(unittest.TestCase):
    def run_split(self, directory, data, *, compressed_input=False,
                  compressed_eligible=False, compressed_oversized=False):
        input_path = directory / ("trimmed.fastq.gz" if compressed_input else "trimmed.fastq")
        eligible = directory / ("eligible.fastq.gz" if compressed_eligible else "eligible.fastq")
        oversized = directory / ("over_6000bp.fastq.gz" if compressed_oversized else "over_6000bp.fastq")
        report = directory / "sample_22bp.length_filter.tsv"
        if compressed_input:
            with gzip.open(input_path, "wb") as handle:
                handle.write(data)
        else:
            input_path.write_bytes(data)
        result = subprocess.run(
            [sys.executable, str(SPLITTER), str(input_path), str(eligible), str(oversized), str(report)],
            capture_output=True, text=True,
        )
        return result, eligible, oversized, report

    def metrics(self, path):
        lines = path.read_text().splitlines()
        self.assertEqual(lines[0], "metric\tcount")
        return {key: int(value) for key, value in
                (line.split("\t") for line in lines[1:])}

    def test_boundary_and_exact_record_preservation_for_plain_and_gzip(self):
        below = fastq_record("duplicate header with spaces and metadata", 5999,
                             plus="+duplicate header with spaces and metadata")
        boundary = fastq_record("duplicate header with spaces and metadata", 6000)
        above = fastq_record("over limit with metadata", 6001)
        short = fastq_record("short", 22)
        # Interleave dispositions to catch dropped/reordered records and avoid
        # an assertion based solely on output file size or header identifiers.
        data = below + above + boundary + short
        for compressed_input in (False, True):
            for compressed_eligible in (False, True):
                for compressed_oversized in (False, True):
                    with self.subTest(input_gzip=compressed_input, eligible_gzip=compressed_eligible,
                                      oversized_gzip=compressed_oversized), tempfile.TemporaryDirectory() as temp:
                        result, eligible, oversized, report = self.run_split(
                            Path(temp), data, compressed_input=compressed_input,
                            compressed_eligible=compressed_eligible,
                            compressed_oversized=compressed_oversized,
                        )
                        self.assertEqual(result.returncode, 0, result.stderr)
                        self.assertEqual(read_fastq_bytes(eligible), below + boundary + short)
                        self.assertEqual(read_fastq_bytes(oversized), above)
                        self.assertEqual(self.metrics(report), {
                            "input": 4, "eligible": 3, "oversized": 1, "max_read_length": 6000,
                        })

    def test_line_endings_and_unterminated_final_quality_are_preserved(self):
        first = fastq_record("crlf comment", 25, newline="\r\n")
        last = fastq_record("last read", 6001, final_newline=False)
        with tempfile.TemporaryDirectory() as temp:
            result, eligible, oversized, _ = self.run_split(Path(temp), first + last)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(read_fastq_bytes(eligible), first)
            self.assertEqual(read_fastq_bytes(oversized), last)

    def test_empty_and_all_excluded_inputs_create_both_outputs(self):
        records = fastq_record("long", 6001) + fastq_record("much longer", 10000)
        for data, count in ((b"", 0), (records, 2)):
            for compressed in (False, True):
                with self.subTest(records=count, gzip=compressed), tempfile.TemporaryDirectory() as temp:
                    result, eligible, oversized, report = self.run_split(
                        Path(temp), data, compressed_input=compressed,
                        compressed_eligible=compressed, compressed_oversized=compressed,
                    )
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertEqual(read_fastq_bytes(eligible), b"")
                    self.assertEqual(read_fastq_bytes(oversized), data)
                    self.assertEqual(self.metrics(report), {
                        "input": count, "eligible": 0, "oversized": count, "max_read_length": 6000,
                    })

    def test_malformed_fastq_is_rejected(self):
        malformed = {
            "missing header marker": b"read\nACGT\n+\n!!!!\n",
            "missing separator marker": b"@read\nACGT\nseparator\n!!!!\n",
            "quality too short": b"@read\nACGT\n+\n!!!\n",
            "quality too long": b"@read\nACGT\n+\n!!!!!\n",
            "truncated sequence": b"@read\n",
            "truncated separator": b"@read\nACGT\n",
            "truncated qualities": b"@read\nACGT\n+\n",
        }
        for description, data in malformed.items():
            with self.subTest(description=description), tempfile.TemporaryDirectory() as temp:
                result, _, _, _ = self.run_split(Path(temp), data)
                self.assertNotEqual(result.returncode, 0)
                self.assertTrue(result.stderr.strip(), "Invalid FASTQ should include an error diagnostic")


if __name__ == "__main__":
    unittest.main()
