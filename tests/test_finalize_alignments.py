"""Exercise RNA/DNA SAM merging and chronological read conservation."""

import sys
import tempfile
import unittest
from pathlib import Path

try:
    import pysam
except ImportError:
    pysam = None

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
if pysam is not None:
    import finalize_sample_alignments as finalizer
    from read_breakdown_csv import bbsplit_counts


def write_sam(path, references, records):
    """Records are (query name, reference name or None, SAM flag)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    header = {"HD": {"VN": "1.6"},
              "SQ": [{"SN": name, "LN": length} for name, length in references]}
    with pysam.AlignmentFile(str(path), "w", header=header) as sam:
        for name, reference, flag in records:
            read = pysam.AlignedSegment(sam.header)
            read.query_name = name
            read.query_sequence = "ACGT" * 5
            read.query_qualities = pysam.qualitystring_to_array("I" * 20)
            read.flag = flag
            if reference is None:
                read.flag |= 4
                read.reference_id = -1
                read.reference_start = -1
            else:
                read.reference_id = sam.get_tid(reference)
                read.reference_start = 7
                read.mapping_quality = 60
                read.cigarstring = "20M"
            sam.write(read)


def stage_fixture(directory, sample="sample_20bp", *, counts=None):
    """Write all four chronological stages with independent RNA/DNA dictionaries."""
    decoy = directory / "decoy" / f"{sample}.mapped_to_other_bugs.bbsplit.txt"
    decoy.parent.mkdir(parents=True)
    finalizer.write_counts(decoy, 12, 2)
    stages = (("bacteria", "rna", 10, 2), ("host", "rna", 8, 2),
              ("bacteria", "dna", 6, 3), ("host", "dna", 3, 1))
    result = {}
    for organism, molecule, input_count, aligned in stages:
        prefix = organism.upper() + "_" + sample
        path = directory / organism / molecule / f"{prefix}.sam"
        ref = f"{organism}_{molecule}|reference"
        actual_aligned = (counts or {}).get((organism, molecule), aligned)
        write_sam(path, [(ref, 100)], [
            # Repeated IDs represent separate input records and must survive.
            ("duplicate-id", ref, 0) for _ in range(actual_aligned)
        ] + [("unmapped", None, 4) for _ in range(input_count - actual_aligned)])
        report = path.with_suffix(".bbsplit.txt")
        finalizer.write_counts(report, input_count, aligned)
        result[organism, molecule] = path, report
    return result


@unittest.skipIf(pysam is None, "pysam is required for SAM finalization tests")
class CombineMappedTests(unittest.TestCase):
    def test_remaps_different_reference_dictionaries_and_preserves_repeated_ids(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            rna, dna, output = (directory / name for name in ("rna.sam", "dna.sam", "merged.sam"))
            write_sam(rna, [("RNA|second", 100), ("RNA|first", 120)],
                      [("same-id", "RNA|first", 0), ("rna-second", "RNA|second", 16)])
            write_sam(dna, [("DNA|second", 200), ("DNA|first", 220)],
                      [("same-id", "DNA|first", 0), ("dna-second", "DNA|second", 16)])
            self.assertEqual(finalizer.combine_mapped([rna, dna], output), 4)
            with pysam.AlignmentFile(str(output), "r") as sam:
                records = list(sam.fetch(until_eof=True))
                self.assertEqual([(read.query_name, read.reference_name) for read in records], [
                    ("same-id", "RNA|first"), ("rna-second", "RNA|second"),
                    ("same-id", "DNA|first"), ("dna-second", "DNA|second"),
                ])
                self.assertEqual([read.reference_start for read in records], [7] * 4)
                self.assertEqual([read.is_reverse for read in records], [False, True, False, True])
                self.assertEqual(dict(zip(sam.references, sam.lengths)), {
                    "RNA|second": 100, "RNA|first": 120, "DNA|second": 200, "DNA|first": 220,
                })

    def test_only_primary_mapped_alignments_are_kept(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source, output = directory / "input.sam", directory / "merged.sam"
            write_sam(source, [("DNA|chr", 100)], [
                ("primary", "DNA|chr", 0), ("secondary", "DNA|chr", 256),
                ("supplementary", "DNA|chr", 2048), ("unmapped", None, 4),
            ])
            self.assertEqual(finalizer.combine_mapped([source], output), 1)
            with pysam.AlignmentFile(str(output), "r") as sam:
                self.assertEqual([read.query_name for read in sam.fetch(until_eof=True)], ["primary"])

    def test_conflicting_reference_lengths_cannot_share_an_output_dictionary(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            first, second, output = (directory / name for name in ("first.sam", "second.sam", "merged.sam"))
            write_sam(first, [("same", 100)], [])
            write_sam(second, [("same", 101)], [])
            with self.assertRaisesRegex(ValueError, "Conflicting reference"):
                finalizer.combine_mapped([first, second], output)


@unittest.skipIf(pysam is None, "pysam is required for SAM finalization tests")
class FinalizationTests(unittest.TestCase):
    def test_rna_first_chain_produces_balanced_virtual_organism_reports(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            stage_fixture(directory)
            finalizer.finalize(directory, "sample_20bp", host_enabled=True, bacteria_rna=True, host_rna=True)
            bacteria_report = directory / "bacteria/BACTERIA_sample_20bp.bbsplit.txt"
            host_report = directory / "host/HOST_sample_20bp.bbsplit.txt"
            self.assertEqual(bbsplit_counts(bacteria_report), (10, 5))
            self.assertEqual(bbsplit_counts(host_report), (5, 3))
            for organism, expected in (("bacteria", 5), ("host", 3)):
                with pysam.AlignmentFile(str(directory / organism / f"{organism.upper()}_sample_20bp.sam"), "r") as sam:
                    reads = list(sam.fetch(until_eof=True))
                self.assertEqual(len(reads), expected)
                self.assertTrue(all(read.query_name == "duplicate-id" for read in reads))
                self.assertEqual({read.reference_name for read in reads}, {
                    f"{organism}_rna|reference", f"{organism}_dna|reference",
                })

    def test_wrong_chronological_stage_input_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            paths = stage_fixture(directory)
            finalizer.write_counts(paths["bacteria", "dna"][1], 8, 3)
            with self.assertRaisesRegex(ValueError, "Stage chain is inconsistent.*bacteria dna"):
                finalizer.finalize(directory, "sample_20bp", host_enabled=True, bacteria_rna=True, host_rna=True)

    def test_mapped_sam_count_that_differs_from_stage_reports_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            stage_fixture(directory, counts={("bacteria", "rna"): 1})
            with self.assertRaisesRegex(ValueError, "Mapped SAM count differs"):
                finalizer.finalize(directory, "sample_20bp", host_enabled=True, bacteria_rna=True, host_rna=True)

    def test_opposite_rna_and_dna_count_errors_cannot_cancel(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            # Aggregate bacterial total remains five, but each molecule report
            # disagrees with its SAM: the RNA/DNA categories would be corrupted.
            stage_fixture(directory, counts={("bacteria", "rna"): 1, ("bacteria", "dna"): 4})
            with self.assertRaisesRegex(ValueError, "Mapped SAM count differs"):
                finalizer.finalize(directory, "sample_20bp", host_enabled=True, bacteria_rna=True, host_rna=True)

    def test_dna_only_empty_sample_still_produces_valid_header_and_zero_report(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            decoy = directory / "decoy/sample_20bp.mapped_to_other_bugs.bbsplit.txt"
            decoy.parent.mkdir(parents=True)
            finalizer.write_counts(decoy, 0, 0)
            dna = directory / "bacteria/dna/BACTERIA_sample_20bp.sam"
            write_sam(dna, [("DNA|chr", 100)], [])
            finalizer.write_counts(dna.with_suffix(".bbsplit.txt"), 0, 0)
            finalizer.finalize(directory, "sample_20bp")
            self.assertEqual(bbsplit_counts(directory / "bacteria/BACTERIA_sample_20bp.bbsplit.txt"), (0, 0))
            with pysam.AlignmentFile(str(directory / "bacteria/BACTERIA_sample_20bp.sam"), "r") as sam:
                self.assertEqual(sam.references, ("DNA|chr",))
                self.assertEqual(list(sam.fetch(until_eof=True)), [])

    def test_host_rna_requires_host_alignment(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(ValueError, "Host RNA cannot be enabled"):
                finalizer.finalize(Path(temporary), "sample_20bp", host_rna=True)


if __name__ == "__main__":
    unittest.main()
