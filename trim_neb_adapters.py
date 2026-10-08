#!/usr/bin/env python3
"""Extract NEB E7330 inserts with Cutadapt, then clean only exposed read ends.

Use the same adapter parser, AdapterCutter (--times 2) and ReverseComplementer
(--revcomp) as Cutadapt's CLI. No extra dependencies beyond Cutadapt are needed.
"""

import argparse
import csv
from collections import Counter

import cutadapt
import dnaio
from cutadapt.info import ModificationInfo
from cutadapt.modifiers import AdapterCutter, ReverseComplementer
from cutadapt.parser import make_adapter


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


class NebTrimmer:
    def __init__(self, five, three, nanopore, error_rate, overlap, rounds,
                 trim_nanopore=True, truseq=None):
        parameters = {'max_errors': error_rate, 'min_overlap': overlap}
        self.rounds = rounds

        def adapter(name, sequence, kind):
            return make_adapter(sequence, kind, parameters, name=name)

        # A linked -g requires both flanks; the rightmost -a boundary keeps
        # internal 3' motifs and exposes repeated terminal copies for cleanup.
        self.pair = adapter('neb_pair', f'{five}...{three};rightmost', 'front')
        self.pair_cutter = ReverseComplementer(
            AdapterCutter([self.pair], times=1), rc_suffix=None)
        neb_terminal = [adapter('neb_5p', f'X{five}', 'front'),
                        adapter('neb_3p', f'{three}X', 'back')]
        self.neb_cutter = AdapterCutter(neb_terminal, times=2)
        self.partial_cutter = ReverseComplementer(self.neb_cutter, rc_suffix=None)
        extras = []
        if trim_nanopore:
            for name, sequence in (('nanopore', nanopore),
                                   ('nanopore_rc', reverse_complement(nanopore))):
                extras.extend([adapter(name + '_5p', f'X{sequence}', 'front'),
                               adapter(name + '_3p', f'{sequence}X', 'back')])
        optional = []
        if truseq:
            optional = [adapter('truseq_3p', f'{truseq}X', 'back'),
                        adapter('truseq_rc_5p', f'X{reverse_complement(truseq)}', 'front')]
        # Nanopore cleanup is only available for incompletely bounded reads.
        self.complete_cutter = AdapterCutter(neb_terminal + optional, times=2)
        self.incomplete_cutter = AdapterCutter(neb_terminal + extras + optional, times=2)
        self.outer_cutter = AdapterCutter(extras + optional, times=2) if extras or optional else None
        self.five = adapter('neb_5p', five, 'front')
        self.three = adapter('neb_3p', three, 'back')
        sequences = [five, three] + ([nanopore] if trim_nanopore else []) + ([truseq] if truseq else [])
        self.detectors = [adapter('detect', seq, 'back') for original in sequences
                          for seq in (original, reverse_complement(original))]

    def concatemer(self, sequence):
        """An internal 3' -> 5' junction suggests multiple inserts; never split."""
        back = self.three.match_to(sequence)
        if back is None or back.rstop == len(sequence):
            return False
        return self.five.match_to(sequence[back.rstop:]) is not None

    def trim(self, original):
        info = ModificationInfo(original)
        read = self.pair_cutter(original, info)
        oriented = bool(info.matches)
        was_rc = bool(info.is_rc)
        counts = Counter()
        rounds = 0
        if info.matches:
            counts.update({'neb_5p': 1, 'neb_3p': 1})
            rounds = 1
        ambiguous = self.concatemer(read.sequence) or (
            not oriented and self.concatemer(reverse_complement(read.sequence)))
        while rounds < self.rounds:
            info = ModificationInfo(read)
            if not oriented:
                candidate = self.partial_cutter(read, info)
                if info.matches:
                    oriented = True
                    was_rc = bool(info.is_rc)
                    read = candidate
                elif self.outer_cutter:
                    read = self.outer_cutter(read, info)
            else:
                cutter = self.complete_cutter if counts['neb_5p'] and counts['neb_3p'] else self.incomplete_cutter
                read = cutter(read, info)
            if not info.matches:
                break
            counts.update(match.adapter.name for match in info.matches)
            rounds += 1
        both = bool(counts['neb_5p'] and counts['neb_3p'])
        # Observed interior matches are diagnostic only, never instructions to
        # trim further. They can be biological motifs, not necessarily adapters.
        remaining = any(a.match_to(read.sequence) is not None for a in self.detectors)
        cutter = self.complete_cutter if both else self.incomplete_cutter
        terminal_remaining = cutter.adapters.match_to(read.sequence) is not None
        if not oriented:
            terminal_remaining = terminal_remaining or (
                self.neb_cutter.adapters.match_to(reverse_complement(read.sequence)) is not None)
        ambiguous = ambiguous or (both and self.concatemer(read.sequence))
        return read, {
            'flanks': 'both' if both else '5p_only' if counts['neb_5p'] else '3p_only' if counts['neb_3p'] else 'none',
            'reverse_complemented': int(was_rc),
            'neb_5p_matches': counts['neb_5p'], 'neb_3p_matches': counts['neb_3p'],
            'adapter_matches': sum(counts.values()),
            'repeated_neb_adapters': int(counts['neb_5p'] > 1 or counts['neb_3p'] > 1),
            'repeated_adapter_sequences': int(
                any(n > 1 for n in counts.values()) or
                counts['nanopore_5p'] + counts['nanopore_3p'] > 1 or
                counts['nanopore_rc_5p'] + counts['nanopore_rc_3p'] > 1),
            'trimming_rounds': rounds, 'ambiguous_concatemer': int(ambiguous),
            'detectable_adapter_motif_remaining': int(remaining),
            'terminal_adapter_remaining': int(terminal_remaining),
            'round_limit_reached': int(terminal_remaining and rounds == self.rounds),
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--report-prefix', required=True)
    parser.add_argument('--five', default='GTTCAGAGTTCTACAGTCCGACGATC')
    parser.add_argument('--three', default='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC')
    parser.add_argument('--nanopore', default='TTTCTGTTGGTGCTGATATTGC')
    parser.add_argument('--nanopore-rc')
    parser.add_argument('--truseq')
    parser.add_argument('--no-nanopore', action='store_true')
    parser.add_argument('--error-rate', type=float, default=.15)
    parser.add_argument('--overlap', type=int, default=12)
    parser.add_argument('--max-rounds', type=int, default=10)
    parser.add_argument('--minimum-length', type=int, default=20)
    args = parser.parse_args()
    if not 0 <= args.error_rate < 1 or args.overlap < 1 or args.max_rounds < 1 or args.minimum_length < 0:
        parser.error('require 0 <= error-rate < 1, positive overlap/max-rounds and nonnegative minimum-length')
    enabled_sequences = [args.five, args.three, args.truseq]
    if not args.no_nanopore:
        enabled_sequences.append(args.nanopore)
    for sequence in enabled_sequences:
        if sequence is not None and (not sequence or set(sequence) - set('ACGT')):
            parser.error('adapter sequences must be nonempty uppercase DNA (ACGT)')
        if sequence is not None and len(sequence) < args.overlap:
            parser.error('minimum overlap cannot exceed an enabled adapter length')
    if not args.no_nanopore and args.nanopore_rc and args.nanopore_rc != reverse_complement(args.nanopore):
        parser.error('nanopore-rc must be the reverse complement of nanopore')
    trimmer = NebTrimmer(args.five, args.three, args.nanopore, args.error_rate,
                         args.overlap, args.max_rounds, not args.no_nanopore, args.truseq)
    stats = Counter({name: 0 for name in (
        'total_input_reads', 'reads_with_adapters_detected', 'trimmed_reads',
        'reads_with_both_flanks', 'reads_with_one_flank', 'reverse_complemented_reads',
        'reads_with_repeated_adapters', 'reads_requiring_multiple_rounds',
        'reads_discarded_by_length', 'reads_with_detectable_adapters_remaining',
        'reads_with_terminal_adapters_remaining', 'reads_reaching_round_limit',
        'ambiguous_concatemer_reads', 'final_retained_reads',
        'retained_reads_with_detectable_adapters_remaining')})
    lengths = Counter()
    fields = ['read_number', 'read_header', 'input_length', 'output_length', 'retained',
              'flanks', 'reverse_complemented', 'neb_5p_matches', 'neb_3p_matches',
              'adapter_matches', 'repeated_neb_adapters', 'repeated_adapter_sequences', 'trimming_rounds',
              'ambiguous_concatemer', 'detectable_adapter_motif_remaining',
              'terminal_adapter_remaining', 'round_limit_reached']
    with dnaio.open(args.input) as source, dnaio.open(args.output, mode='w', fileformat='fastq') as output, \
            open(args.report_prefix + '.adapter_reads.csv', 'w', newline='') as detail:
        writer = csv.DictWriter(detail, fieldnames=fields)
        writer.writeheader()
        for number, original in enumerate(source, 1):
            if original.qualities is None:
                raise ValueError('FASTQ input with quality scores is required')
            read, evidence = trimmer.trim(original)
            retained = len(read) >= args.minimum_length
            stats['total_input_reads'] += 1
            stats['reads_with_adapters_detected'] += bool(evidence['adapter_matches'] or evidence['detectable_adapter_motif_remaining'])
            stats['trimmed_reads'] += len(read) != len(original)
            stats['reads_with_both_flanks'] += evidence['flanks'] == 'both'
            stats['reads_with_one_flank'] += evidence['flanks'] in ('5p_only', '3p_only')
            stats['reverse_complemented_reads'] += evidence['reverse_complemented']
            stats['reads_with_repeated_adapters'] += evidence['repeated_adapter_sequences']
            stats['reads_requiring_multiple_rounds'] += evidence['trimming_rounds'] > 1
            stats['reads_discarded_by_length'] += not retained
            stats['reads_with_detectable_adapters_remaining'] += evidence['detectable_adapter_motif_remaining']
            stats['reads_with_terminal_adapters_remaining'] += evidence['terminal_adapter_remaining']
            stats['reads_reaching_round_limit'] += evidence['round_limit_reached']
            stats['ambiguous_concatemer_reads'] += evidence['ambiguous_concatemer']
            stats['final_retained_reads'] += retained
            stats['retained_reads_with_detectable_adapters_remaining'] += retained and evidence['detectable_adapter_motif_remaining']
            writer.writerow(dict(read_number=number, read_header=original.name,
                                 input_length=len(original), output_length=len(read),
                                 retained=int(retained), **evidence))
            if retained:
                output.write(read)
                lengths[len(read)] += 1
    distribution = sorted(lengths.items())
    stats['retained_length_min'] = distribution[0][0] if distribution else 0
    stats['retained_length_max'] = distribution[-1][0] if distribution else 0
    stats['retained_length_mean'] = sum(n*c for n, c in distribution) / stats['final_retained_reads'] if distribution else 0
    with open(args.report_prefix + '.adapter_stats.csv', 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['cutadapt_version', 'statistics_source'] + list(stats))
        writer.writeheader()
        writer.writerow(dict(cutadapt_version=cutadapt.__version__, statistics_source='per-read Cutadapt match objects; adapter_reads.csv', **stats))
    with open(args.report_prefix + '.read_lengths.csv', 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['length', 'retained_reads'])
        writer.writerows(distribution)
    print(f'Cutadapt {cutadapt.__version__} matching engine; NEB E7330 mode')
    print('Linked -g 5P...3P;rightmost with --revcomp; terminal -g X5P / -a 3PX with --times 2')
    print(f'Error rate: {args.error_rate}; minimum overlap: {args.overlap}; maximum rounds: {args.max_rounds}')
    # These compatibility lines describe the entire operation, including the
    # final filter, rather than pretending they are a native one-pass report.
    print(f'Total reads processed: {stats["total_input_reads"]:,}')
    print(f'Reads written (passing filters): {stats["final_retained_reads"]:,}')
    for name, value in stats.items():
        print(f'{name}: {value}')


if __name__ == '__main__':
    main()
