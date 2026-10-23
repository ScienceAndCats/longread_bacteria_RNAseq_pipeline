#!/usr/bin/env python3
"""Validate the RNA-first stage chain and combine each organism's mapped reads."""

import argparse
import csv
from pathlib import Path

import pysam

from read_breakdown_csv import bbsplit_counts


def combine_mapped(inputs, destination):
    """Merge coordinate dictionaries without merging RNA and genomic coordinates."""
    sequences = {}
    for path in inputs:
        with pysam.AlignmentFile(str(path), "r") as source:
            for record in source.header.to_dict().get("SQ", []):
                name = record["SN"]
                if name in sequences and sequences[name] != record:
                    raise ValueError(f"Conflicting reference sequence {name} in {path}")
                sequences[name] = record
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": list(sequences.values()),
        "PG": [{"ID": "finalize_sample_alignments", "PN": "finalize_sample_alignments.py"}],
    })
    mapped = 0
    with pysam.AlignmentFile(str(destination), "w", header=header) as output:
        for path in inputs:
            with pysam.AlignmentFile(str(path), "r") as source:
                for read in source.fetch(until_eof=True):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    # SAM text resolves reference names against the combined dictionary.
                    output.write(pysam.AlignedSegment.fromstring(read.to_string(), header))
                    mapped += 1
    return mapped


def write_counts(path, input_count, aligned):
    if not 0 <= aligned <= input_count:
        raise ValueError(f"Invalid aggregate alignment count for {path}")
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerows((("metric", "count"), ("input", input_count),
                          ("aligned", aligned), ("unmapped", input_count - aligned)))


def finalize(alignment_dir, sample, *, host_enabled=False, bacteria_rna=False, host_rna=False):
    if host_rna and not host_enabled:
        raise ValueError("Host RNA cannot be enabled without host alignment")
    decoy_input, decoy_aligned = bbsplit_counts(
        alignment_dir / f"decoy/{sample}.mapped_to_other_bugs.bbsplit.txt")
    after_decoy = decoy_input - decoy_aligned
    remaining = after_decoy
    stages = []
    if bacteria_rna:
        stages.append(("bacteria", "rna"))
    if host_rna:
        stages.append(("host", "rna"))
    stages.append(("bacteria", "dna"))
    if host_enabled:
        stages.append(("host", "dna"))
    paths = {"bacteria": [], "host": []}
    aligned = {"bacteria": 0, "host": 0}
    for organism, molecule in stages:
        prefix = organism.upper() + "_" + sample
        stage_dir = alignment_dir / organism / molecule
        stage_input, stage_aligned = bbsplit_counts(stage_dir / f"{prefix}.bbsplit.txt")
        if stage_input != remaining:
            raise ValueError(f"Stage chain is inconsistent for {sample}: {organism} {molecule}")
        sam_path = stage_dir / f"{prefix}.sam"
        with pysam.AlignmentFile(str(sam_path), "r") as source:
            mapped = sum(1 for read in source.fetch(until_eof=True)
                         if not (read.is_unmapped or read.is_secondary or read.is_supplementary))
        if mapped != stage_aligned:
            raise ValueError(f"Mapped SAM count differs from stage report for {sam_path}")
        remaining -= stage_aligned
        paths[organism].append(sam_path)
        aligned[organism] += stage_aligned
    # Aggregate reports retain the historical organism-level accounting:
    # bacteria input = after decoy; host input = after decoy minus all bacteria.
    # The molecule-stage reports above record the actual chronological inputs.
    for organism in ("bacteria", "host") if host_enabled else ("bacteria",):
        prefix = organism.upper() + "_" + sample
        combined = alignment_dir / organism / f"{prefix}.sam"
        mapped = combine_mapped(paths[organism], combined)
        if mapped != aligned[organism]:
            raise ValueError(f"Mapped SAM count differs from stage reports for {combined}")
        logical_input = after_decoy if organism == "bacteria" else after_decoy - aligned["bacteria"]
        write_counts(combined.with_suffix(".bbsplit.txt"), logical_input, mapped)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment-dir", type=Path, required=True)
    parser.add_argument("--sample-prefix", required=True)
    parser.add_argument("--host-enabled", action="store_true")
    parser.add_argument("--bacteria-rna", action="store_true")
    parser.add_argument("--host-rna", action="store_true")
    args = parser.parse_args()
    finalize(args.alignment_dir, args.sample_prefix, host_enabled=args.host_enabled,
             bacteria_rna=args.bacteria_rna, host_rna=args.host_rna)


if __name__ == "__main__":
    main()
