#!/usr/bin/env python3
"""Separate adapter-trimmed reads above BBSplit's 6,000 bp alignment limit."""

import argparse
import csv
import gzip
from contextlib import ExitStack
from pathlib import Path


MAX_READ_LENGTH = 6000


def open_fastq(path, mode):
    opener = gzip.open if path.suffix == ".gz" else open
    return opener(path, mode, encoding="utf-8", newline="")


def split_fastq(input_path, eligible_path, oversized_path, report_path):
    paths = [Path(path) for path in (input_path, eligible_path, oversized_path, report_path)]
    if len({path.resolve() for path in paths}) != len(paths):
        raise ValueError("Input, output, and report paths must be distinct")
    input_path, eligible_path, oversized_path, report_path = paths
    counts = {"input": 0, "eligible": 0, "oversized": 0, "max_read_length": MAX_READ_LENGTH}
    with ExitStack() as stack:
        source = stack.enter_context(open_fastq(input_path, "rt"))
        eligible = stack.enter_context(open_fastq(eligible_path, "wt"))
        oversized = stack.enter_context(open_fastq(oversized_path, "wt"))
        while True:
            header = source.readline()
            if not header:
                break
            record = [header, source.readline(), source.readline(), source.readline()]
            number = counts["input"] + 1
            if any(not line for line in record):
                raise ValueError(f"Incomplete FASTQ record {number} in {input_path}")
            sequence = record[1].rstrip("\r\n")
            qualities = record[3].rstrip("\r\n")
            if not header.startswith("@") or not record[2].startswith("+"):
                raise ValueError(f"Invalid FASTQ record {number} in {input_path}")
            if len(sequence) != len(qualities):
                raise ValueError(f"Sequence/quality length mismatch in FASTQ record {number} in {input_path}")
            bucket = "oversized" if len(sequence) > MAX_READ_LENGTH else "eligible"
            (oversized if bucket == "oversized" else eligible).writelines(record)
            counts["input"] += 1
            counts[bucket] += 1

    with report_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "count"))
        writer.writerows(counts.items())
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_fastq", type=Path)
    parser.add_argument("eligible_fastq", type=Path)
    parser.add_argument("oversized_fastq", type=Path)
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    split_fastq(args.input_fastq, args.eligible_fastq, args.oversized_fastq, args.report)


if __name__ == "__main__":
    main()
