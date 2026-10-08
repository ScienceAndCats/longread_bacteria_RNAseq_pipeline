"""Synthetic molecule/FASTQ tests against Cutadapt's actual matching engine."""
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import dnaio

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from trim_neb_adapters import NebTrimmer, reverse_complement

FIVE = 'GTTCAGAGTTCTACAGTCCGACGATC'
THREE = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
NP = 'TTTCTGTTGGTGCTGATATTGC'
TRUSEQ = 'TGGAATTCTCGGGTGCCAAGG'
INSERT = 'ACCTGACTAGTCAGGTCATGAC'
P5 = 'AATGATACGGCGACCACCGAGATCTACAC'
P7 = 'ATCTCGTATGCCGTCTTCTGCTTG'


def molecule(name, sequence):
    # Vary qualities to catch reversal/slicing mistakes, not just length changes.
    qualities = ''.join(chr(40 + i % 35) for i in range(len(sequence)))
    return dnaio.SequenceRecord(name, sequence, qualities)


class MatchingTests(unittest.TestCase):
    def setUp(self):
        self.trimmer = NebTrimmer(FIVE, THREE, NP, .15, 12, 10)

    def assert_insert(self, sequence, insert=INSERT, rc=False, trimmer=None):
        original = molecule('read identifier with comment', sequence)
        start = sequence.index(insert)
        expected_quality = original.qualities[start:start+len(insert)]
        if rc:
            original = original.reverse_complement()
        read, evidence = (trimmer or self.trimmer).trim(original)
        self.assertEqual(read.sequence, insert)
        self.assertEqual(read.qualities, expected_quality)
        self.assertEqual(read.name, original.name)
        self.assertEqual(evidence['reverse_complemented'], int(rc))
        return evidence

    def test_single_and_repeated_neb_copies_both_orientations(self):
        for five_count, three_count in ((1, 1), (2, 1), (3, 1), (1, 2), (1, 3), (2, 2), (3, 3)):
            for rc in (False, True):
                with self.subTest(five=five_count, three=three_count, reverse=rc):
                    evidence = self.assert_insert(FIVE*five_count + INSERT + THREE*three_count, rc=rc)
                    self.assertEqual(evidence['flanks'], 'both')
                    self.assertEqual(evidence['neb_5p_matches'], five_count)
                    self.assertEqual(evidence['neb_3p_matches'], three_count)
                    self.assertEqual(evidence['repeated_adapter_sequences'], int(max(five_count, three_count) > 1))
                    self.assertEqual(evidence['ambiguous_concatemer'], 0)

    def test_outer_pcr_indexes_p5_p7_and_nanopore(self):
        for rc in (False, True):
            sequence = NP + P5 + 'AACCTTGG' + FIVE*2 + INSERT + THREE*3 + P7 + reverse_complement(NP)
            self.assert_insert(sequence, rc=rc)

    def test_adapter_substitutions_insertions_deletions(self):
        # Mutations are away from the insert junction, whose bases remain known.
        variants = (
            (FIVE[:7] + 'A' + FIVE[8:], THREE[:9] + 'T' + THREE[10:]),
            (FIVE[:8] + 'A' + FIVE[8:], THREE[:10] + 'C' + THREE[10:]),
            (FIVE[:8] + FIVE[9:], THREE[:10] + THREE[11:]),
        )
        for front, back in variants:
            for rc in (False, True):
                with self.subTest(front=front, back=back, reverse=rc):
                    self.assert_insert(P5 + front + INSERT + back + P7, rc=rc)

    def test_terminal_partial_adapters_and_overlap(self):
        for sequence in (FIVE[-12:] + INSERT + THREE[:12], FIVE[-12:] + INSERT, INSERT + THREE[:12]):
            for rc in (False, True):
                with self.subTest(sequence=sequence, reverse=rc):
                    self.assert_insert(sequence, rc=rc)
        original = molecule('under_overlap', FIVE[-11:] + INSERT)
        trimmed, evidence = self.trimmer.trim(original)
        self.assertEqual(trimmed.sequence, original.sequence)
        self.assertEqual(evidence['flanks'], 'none')

    def test_one_flank_not_reported_as_both(self):
        for sequence, flank in ((FIVE + INSERT, '5p_only'), (INSERT + THREE, '3p_only')):
            for rc in (False, True):
                evidence = self.assert_insert(sequence, rc=rc)
                self.assertEqual(evidence['flanks'], flank)

    def test_terminal_nanopore_on_incomplete_libraries(self):
        for np in (NP, reverse_complement(NP)):
            for rc in (False, True):
                self.assert_insert(np + FIVE + INSERT, rc=rc)
                self.assert_insert(INSERT + THREE + np, rc=rc)

    def test_no_adapters_and_internal_motifs_preserved(self):
        for insert in (INSERT, INSERT + FIVE + INSERT, INSERT + THREE + INSERT,
                       INSERT + NP + INSERT, INSERT + FIVE + INSERT + THREE + INSERT):
            with self.subTest(insert=insert):
                if FIVE in insert and THREE in insert:
                    # With true outer flanks, an internal same-orientation pair
                    # must not be searched/trimmed again after bounding.
                    self.assert_insert(FIVE + insert + THREE, insert=insert)
                else:
                    original = molecule('unbounded', insert)
                    read, evidence = self.trimmer.trim(original)
                    self.assertEqual(read.sequence, insert)
                    self.assertEqual(read.qualities, original.qualities)
                    self.assertEqual(evidence['adapter_matches'], 0)
                    self.assert_insert(FIVE + insert + THREE, insert=insert)

    def test_insert_ending_in_nanopore_motif_not_retrimmed(self):
        insert = INSERT + NP
        self.assert_insert(FIVE + insert + THREE, insert=insert)

    def test_dimers(self):
        for sequence in (FIVE + THREE, FIVE*3 + THREE*3):
            for rc in (False, True):
                read, evidence = self.trimmer.trim(molecule('dimer', reverse_complement(sequence) if rc else sequence))
                self.assertEqual(read.sequence, '')
                self.assertEqual(evidence['flanks'], 'both')

    def test_concatemer_flagged_unsplit(self):
        insert = INSERT + THREE + FIVE + INSERT
        for rc in (False, True):
            evidence = self.assert_insert(FIVE + insert + THREE, insert=insert, rc=rc)
            self.assertEqual(evidence['ambiguous_concatemer'], 1)
            self.assertEqual(evidence['detectable_adapter_motif_remaining'], 1)
            original = molecule('unbounded_concatemer', reverse_complement(insert) if rc else insert)
            read, evidence = self.trimmer.trim(original)
            self.assertEqual(read.sequence, original.sequence)
            self.assertEqual(evidence['ambiguous_concatemer'], 1)

    def test_round_limit_is_enforced_and_reported(self):
        trimmer = NebTrimmer(FIVE, THREE, NP, .15, 12, 1)
        read, evidence = trimmer.trim(molecule('limited', FIVE*3 + INSERT + THREE*3))
        self.assertEqual(read.sequence, FIVE*2 + INSERT + THREE*2)
        self.assertEqual(evidence['trimming_rounds'], 1)
        self.assertEqual(evidence['round_limit_reached'], 1)
        sequence = NP + reverse_complement(FIVE + INSERT) + NP
        read, evidence = trimmer.trim(molecule('limited_unoriented', sequence))
        self.assertEqual(read.sequence, reverse_complement(FIVE + INSERT))
        self.assertEqual(evidence['reverse_complemented'], 0)
        self.assertEqual(evidence['round_limit_reached'], 1)

    def test_nanopore_optional_and_truseq_opt_in(self):
        trimmer = NebTrimmer(FIVE, THREE, NP, .15, 12, 10, trim_nanopore=False)
        read, evidence = trimmer.trim(molecule('np_only', NP + INSERT))
        self.assertEqual(read.sequence, NP + INSERT)
        self.assertEqual(evidence['adapter_matches'], 0)
        read, _ = self.trimmer.trim(molecule('truseq_default', INSERT + TRUSEQ))
        self.assertEqual(read.sequence, INSERT + TRUSEQ)
        enabled = NebTrimmer(FIVE, THREE, NP, .15, 12, 10, truseq=TRUSEQ)
        read, _ = enabled.trim(molecule('truseq_enabled', INSERT + TRUSEQ))
        self.assertEqual(read.sequence, INSERT)


class FastqIntegrationTests(unittest.TestCase):
    def run_trim(self, directory, records, mode='neb_e7330', compressed=False, minimum=22, extra=''):
        input_path = directory / ('input.fastq.gz' if compressed else 'input.fastq')
        with dnaio.open(input_path, mode='w') as writer:
            for record in records:
                writer.write(record)
        config = directory / 'config.env'
        config.write_text(f'ADAPTER_TRIMMING_MODE="{mode}"\nMIN_READ_LENGTH="{minimum}"\nTHREADS="1"\n' + extra)
        output = directory / ('output.fastq.gz' if compressed else 'output.fastq')
        report = directory / 'sample_22bp'
        subprocess.run(['bash', str(REPO / 'trim_fastq.sh'), str(config), str(input_path), str(output), str(report)], check=True)
        with dnaio.open(output) as reader:
            results = list(reader)
        return results, report

    def test_fastq_logs_stats_final_length_and_duplicate_identifiers(self):
        for compressed in (False, True):
            with self.subTest(gzipped=compressed), tempfile.TemporaryDirectory() as temp:
                directory = Path(temp)
                constructs = [FIVE+INSERT+THREE, reverse_complement(FIVE+INSERT+THREE),
                              FIVE*3+INSERT+THREE*3, FIVE+INSERT,
                              FIVE+INSERT[:21]+THREE, FIVE+INSERT+'ACGTA'+THREE,
                              FIVE+THREE, INSERT, INSERT+THREE+FIVE+INSERT]
                originals = [molecule('duplicate header with comment', seq) for seq in constructs]
                results, report = self.run_trim(directory, originals, compressed=compressed)
                self.assertEqual([r.sequence for r in results], [INSERT]*4 + [INSERT+'ACGTA', INSERT, constructs[-1]])
                self.assertTrue(all(r.name == 'duplicate header with comment' for r in results))
                self.assertTrue(all(len(r.sequence) == len(r.qualities) for r in results))
                with Path(str(report)+'.adapter_reads.csv').open() as handle:
                    details = list(csv.DictReader(handle))
                self.assertEqual(len(details), 9)
                self.assertEqual([int(r['retained']) for r in details], [1,1,1,1,0,1,0,1,1])
                with Path(str(report)+'.adapter_stats.csv').open() as handle:
                    stats = next(csv.DictReader(handle))
                expected = {'total_input_reads': 9, 'trimmed_reads': 7, 'reads_with_both_flanks': 6,
                            'reads_with_one_flank': 1, 'reverse_complemented_reads': 1,
                            'reads_with_repeated_adapters': 1, 'reads_requiring_multiple_rounds': 1,
                            'reads_discarded_by_length': 2, 'final_retained_reads': 7,
                            'ambiguous_concatemer_reads': 1, 'reads_with_adapters_detected': 8,
                            'reads_with_detectable_adapters_remaining': 1,
                            'reads_with_terminal_adapters_remaining': 0}
                for key, value in expected.items():
                    self.assertEqual(int(stats[key]), value, key)
                text = Path(str(report)+'.cutadapt_log.txt').read_text()
                self.assertIn('Total reads processed: 9', text)
                self.assertIn('Reads written (passing filters): 7', text)
                with Path(str(report)+'.read_lengths.csv').open() as handle:
                    lengths = {int(r['length']): int(r['retained_reads']) for r in csv.DictReader(handle)}
                self.assertEqual(sum(lengths.values()), 7)
                self.assertEqual(lengths[22], 5)

    def test_legacy_optional_illumina_filter_runs_last(self):
        with tempfile.TemporaryDirectory() as temp:
            short = INSERT[:21]
            records = [molecule('valid', INSERT+NP), molecule('short_after_illumina', short+TRUSEQ+NP)]
            results, report = self.run_trim(Path(temp), records, mode='legacy',
                                           extra='TRIM_ILLUMINA_ADAPTERS="true"\n')
            self.assertEqual([(r.name, r.sequence) for r in results], [('valid', INSERT)])
            self.assertIn('Total reads processed:                       2', Path(str(report)+'.cutadapt_log.txt').read_text())

    def test_empty_input_and_all_short_input(self):
        for records in ([], [molecule('dimer', FIVE+THREE)]):
            with tempfile.TemporaryDirectory() as temp:
                results, report = self.run_trim(Path(temp), records)
                self.assertEqual(results, [])
                with Path(str(report)+'.adapter_stats.csv').open() as handle:
                    stats = next(csv.DictReader(handle))
                self.assertEqual(int(stats['total_input_reads']), len(records))
                self.assertEqual(int(stats['final_retained_reads']), 0)

    def test_minimum_twenty_and_configurable_overlap_without_nanopore(self):
        with tempfile.TemporaryDirectory() as temp:
            records = [molecule('short', FIVE+INSERT[:19]+THREE),
                       molecule('minimum', FIVE+INSERT[:20]+THREE),
                       molecule('long', FIVE+INSERT+THREE)]
            results, _ = self.run_trim(Path(temp), records, minimum=20,
                                      extra='CUTADAPT_MIN_OVERLAP="24"\nTRIM_NANOPORE_ADAPTERS="false"\n')
            self.assertEqual([(r.name, r.sequence) for r in results],
                             [('minimum', INSERT[:20]), ('long', INSERT)])


if __name__ == '__main__':
    unittest.main()
