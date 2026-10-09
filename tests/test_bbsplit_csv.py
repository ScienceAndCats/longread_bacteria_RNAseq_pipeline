"""Check read conservation and oversized-read reporting in both summary CSVs."""

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
import read_breakdown_csv


ELIGIBLE = "Reads eligible for alignment (<=6000 bp)"
OVERSIZED = "Reads over 6000 bp (not processed)"


def write_metrics(path, metrics):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("metric\tcount\n" + "".join(f"{key}\t{value}\n" for key, value in metrics.items()))


def fixture(directory, *, host=False, sample="sample", minimum=22, empty=False):
    stem = f"{sample}_{minimum}bp"
    passing, eligible, oversized = (0, 0, 0) if empty else (10, 8, 2)
    total = 0 if empty else 13
    decoy, bacteria, host_aligned = (0, 0, 0) if empty else (2, 3, 2 if host else 0)
    bacteria_rrna = 0 if empty else 1
    host_rrna = 0 if empty else 1
    (directory / f"{stem}.cutadapt_log.txt").write_text(
        f"Total reads processed: {total}\nReads written (passing filters): {passing}\n"
    )
    write_metrics(directory / f"{stem}.length_filter.tsv", {
        "input": passing, "eligible": eligible, "oversized": oversized, "max_read_length": 6000,
    })
    alignment_dir = directory / "bbsplit_alignments"
    decoy_path = alignment_dir / "decoy" / f"{stem}.mapped_to_other_bugs.bbsplit.txt"
    bacteria_path = alignment_dir / "bacteria" / f"BACTERIA_{stem}.bbsplit.txt"
    write_metrics(decoy_path, {"input": eligible, "aligned": decoy, "unmapped": eligible - decoy})
    write_metrics(bacteria_path, {
        "input": eligible - decoy, "aligned": bacteria, "unmapped": eligible - decoy - bacteria,
    })
    (directory / "featurecounts_BACTERIA_rRNA.txt.summary").write_text(
        f"Status\t{alignment_dir / 'bacteria' / f'BACTERIA_{stem}.sam'}\nAssigned\t{bacteria_rrna}\n"
    )
    if host:
        host_dir = alignment_dir / "host"
        host_path = host_dir / f"HOST_{stem}.bbsplit.txt"
        host_input = eligible - decoy - bacteria
        write_metrics(host_path, {"input": host_input, "aligned": host_aligned, "unmapped": host_input - host_aligned})
        (host_dir / "featurecounts_HOST_rRNA.txt.summary").write_text(
            f"Status\t{host_dir / f'HOST_{stem}.sam'}\nAssigned\t{host_rrna}\n"
        )
    # The legacy main summary also uses CDS feature assignments and coverage.
    (directory / "featurecounts_BACTERIA_summary.txt").write_text(
        "# featureCounts synthetic fixture\n"
        f"Geneid\tChr\tStart\tEnd\tStrand\tLength\t{alignment_dir / 'bacteria' / f'BACTERIA_{stem}.sam'}\n"
        f"gene_a\tchr\t1\t25\t+\t25\t{0 if empty else 2}\n"
        "gene_b\tchr\t26\t50\t+\t25\t0\n"
    )
    coverage = alignment_dir / "bacteria" / f"BACTERIA_{stem}_coverage.txt"
    coverage.write_text("" if empty else
                        "coverage\nreference\nlength\npositions\nPercent covered: 20.0\nMean coverage: 3.5\n")
    return {
        "sample": sample, "total": total, "passing": passing, "eligible": eligible,
        "oversized": oversized, "decoy": decoy, "bacteria": bacteria,
        "host": host_aligned, "leftover": eligible - decoy - bacteria - host_aligned,
        "length": directory / f"{stem}.length_filter.tsv", "decoy_path": decoy_path,
        "bacteria_path": bacteria_path,
    }


class ReadBreakdownTests(unittest.TestCase):
    def test_conservation_with_and_without_host_and_configurable_minimum(self):
        for host in (False, True):
            for minimum in (20, 22):
                with self.subTest(host=host, minimum=minimum), tempfile.TemporaryDirectory() as temp:
                    directory = Path(temp)
                    expected = fixture(directory, host=host, minimum=minimum)
                    sample, counts, sources = read_breakdown_csv.build_table(directory)[0]
                    self.assertEqual(sample, expected["sample"])
                    for attribute, key in (("Total reads", "total"),
                                           ("Reads passing cutadapt length filter", "passing"),
                                           (ELIGIBLE, "eligible"), (OVERSIZED, "oversized"),
                                           ("Decoy aligned", "decoy"), ("Bacteria aligned", "bacteria"),
                                           ("Host aligned", "host"), ("Leftover reads", "leftover")):
                        self.assertEqual(counts[attribute], expected[key], attribute)
                    self.assertEqual(counts["Reads passing cutadapt length filter"], counts[ELIGIBLE] + counts[OVERSIZED])
                    self.assertEqual(counts[ELIGIBLE], sum(counts[key] for key in
                                     ("Decoy aligned", "Bacteria aligned", "Host aligned", "Leftover reads")))
                    self.assertEqual(counts["Bacteria aligned"], counts["Bacteria non-rRNA"] + counts["Bacteria rRNA"])
                    self.assertEqual(counts["Host aligned"], counts["Host non-rRNA"] + counts["Host rRNA"])
                    self.assertIn("length_filter.tsv", sources[OVERSIZED])
                    self.assertIn("bbsplit_alignments", sources["Bacteria aligned"])

    def test_empty_and_all_excluded_reads_are_reported_without_leftovers(self):
        for empty in (False, True):
            with self.subTest(empty=empty), tempfile.TemporaryDirectory() as temp:
                directory = Path(temp)
                expected = fixture(directory, empty=empty)
                if not empty:
                    write_metrics(expected["length"], {"input": 10, "eligible": 0, "oversized": 10, "max_read_length": 6000})
                    for report in (expected["decoy_path"], expected["bacteria_path"]):
                        write_metrics(report, {"input": 0, "aligned": 0, "unmapped": 0})
                    (directory / "featurecounts_BACTERIA_rRNA.txt.summary").write_text(
                        "Status\tBACTERIA_sample_22bp.sam\nAssigned\t0\n"
                    )
                _, counts, _ = read_breakdown_csv.build_table(directory)[0]
                self.assertEqual(counts[ELIGIBLE], 0)
                self.assertEqual(counts[OVERSIZED], 0 if empty else 10)
                self.assertEqual(counts["Leftover reads"], 0)

    def test_inconsistent_counts_are_rejected(self):
        invalid_metrics = {
            "length input differs from trim": ("length", {"input": 9, "eligible": 7, "oversized": 2, "max_read_length": 6000}),
            "length counts do not conserve": ("length", {"input": 10, "eligible": 8, "oversized": 1, "max_read_length": 6000}),
            "wrong length limit": ("length", {"input": 10, "eligible": 8, "oversized": 2, "max_read_length": 6001}),
            "oversized reaches decoy stage": ("decoy_path", {"input": 10, "aligned": 2, "unmapped": 8}),
            "bacteria receives wrong total": ("bacteria_path", {"input": 7, "aligned": 3, "unmapped": 4}),
            "stage does not conserve": ("bacteria_path", {"input": 6, "aligned": 3, "unmapped": 2}),
            "negative aligned count": ("bacteria_path", {"input": 6, "aligned": -1, "unmapped": 7}),
        }
        for description, (path_key, metrics) in invalid_metrics.items():
            with self.subTest(description=description), tempfile.TemporaryDirectory() as temp:
                directory = Path(temp)
                paths = fixture(directory)
                write_metrics(paths[path_key], metrics)
                with self.assertRaises(ValueError):
                    read_breakdown_csv.build_table(directory)

    def test_host_input_and_rrna_inconsistencies_are_rejected(self):
        for invalid in ("host input", "bacterial rrna", "host rrna"):
            with self.subTest(invalid=invalid), tempfile.TemporaryDirectory() as temp:
                directory = Path(temp)
                fixture(directory, host=True)
                if invalid == "host input":
                    write_metrics(directory / "bbsplit_alignments/host/HOST_sample_22bp.bbsplit.txt",
                                  {"input": 4, "aligned": 2, "unmapped": 2})
                elif invalid == "bacterial rrna":
                    (directory / "featurecounts_BACTERIA_rRNA.txt.summary").write_text(
                        "Status\tBACTERIA_sample_22bp.sam\nAssigned\t4\n"
                    )
                else:
                    (directory / "bbsplit_alignments/host/featurecounts_HOST_rRNA.txt.summary").write_text(
                        "Status\tHOST_sample_22bp.sam\nAssigned\t3\n"
                    )
                with self.assertRaises(ValueError):
                    read_breakdown_csv.build_table(directory)

    def test_breakdown_csv_quotes_sample_and_reports_percent_of_all_raw_reads(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            fixture(directory, sample="sample, comment")
            output = directory / "breakdown.csv"
            subprocess.run([sys.executable, str(REPO / "read_breakdown_csv.py"),
                            "--output-dir", str(directory), "--output", str(output)], check=True)
            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["Sample"], "sample, comment")
            self.assertEqual(rows[0][f"{OVERSIZED} (reads)"], "2")
            self.assertEqual(rows[0][f"{OVERSIZED} (% of total reads)"], "15.38%")


class MainSummaryTests(unittest.TestCase):
    def test_actual_empty_cutadapt_and_coverage_reports_produce_zero_counts(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            fixture(directory, empty=True)
            # These are the observed Cutadapt 5.2 / samtools 1.24 outputs for
            # an empty FASTQ and a header-only, zero-mapped-read BAM.
            (directory / "sample_22bp.cutadapt_log.txt").write_text(
                "This is cutadapt 5.2\nProcessing single-end reads on 1 core ...\nNo reads processed!\n"
            )
            subprocess.run([sys.executable, str(REPO / "bacteria_with_decoys_csvConversion.py")],
                           cwd=directory, check=True)
            output = directory / f"BACTERIA_{directory.name}.csv"
            with output.open(newline="") as handle:
                row = {key.strip(): value for key, value in next(csv.DictReader(handle)).items()}
            for column in ("Raw Reads", "Reads Written", "Overall Allignment", "Percent Covered", "Mean Coverage",
                           "Reads Assigned to features", "Features with non-zero reads", ELIGIBLE, OVERSIZED):
                self.assertEqual(float(row[column]), 0, column)
            _, counts, _ = read_breakdown_csv.build_table(directory)[0]
            self.assertTrue(all(count == 0 for count in counts.values()))

    def test_empty_coverage_with_mapped_reads_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            fixture(directory)
            (directory / "bbsplit_alignments/bacteria/BACTERIA_sample_22bp_coverage.txt").write_text("")
            result = subprocess.run([sys.executable, str(REPO / "bacteria_with_decoys_csvConversion.py")],
                                    cwd=directory, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)

    def test_main_csv_includes_excluded_count_and_preserves_existing_metrics(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            fixture(directory, sample="sample, comment")
            subprocess.run([sys.executable, str(REPO / "bacteria_with_decoys_csvConversion.py")],
                           cwd=directory, check=True)
            output = directory / f"BACTERIA_{directory.name}.csv"
            with output.open(newline="") as handle:
                rows = [{key.strip(): value for key, value in row.items()} for row in csv.DictReader(handle)]
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertEqual(row["Sample"], "sample, comment")
            for attribute, value in (("Raw Reads", "13"), ("Reads Written", "10"),
                                     ("Overall Allignment", "3"), ("Percent Covered", "20.0"),
                                     ("Mean Coverage", "3.5"), ("Reads Assigned to features", "2"),
                                     ("Features with non-zero reads", "1"), (ELIGIBLE, "8"), (OVERSIZED, "2")):
                self.assertEqual(row[attribute], value, attribute)


if __name__ == "__main__":
    unittest.main()
