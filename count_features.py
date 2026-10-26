#!/usr/bin/env python3
"""Run featureCounts once per SAM and combine compatible sample columns.

Some featureCounts builds abort when processing several SAMs in one invocation.
Independent invocations also let empty samples retain their own zero columns.
"""

import argparse
import csv
import gzip
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path


FEATURE_HEADER = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
# Large genes can produce long semicolon-separated annotation coordinate fields.
CSV_FIELD_SIZE_LIMIT = 64 * 1024 * 1024
SUMMARY_STATUSES = (
    "Assigned", "Unassigned_Unmapped", "Unassigned_Read_Type", "Unassigned_Singleton",
    "Unassigned_MappingQuality", "Unassigned_Chimera", "Unassigned_FragmentLength",
    "Unassigned_Duplicate", "Unassigned_MultiMapping", "Unassigned_Secondary",
    "Unassigned_NonSplit", "Unassigned_NoFeatures", "Unassigned_Overlapping_Length",
    "Unassigned_Ambiguity",
)


def has_feature_type(annotation, feature_type):
    opener = gzip.open if annotation.name.endswith(".gz") else open
    requested = set(feature_type.split(","))
    with opener(annotation, "rt") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 9 and fields[2] in requested:
                return True
    return False


def read_single_table(path, *, summary=False):
    csv.field_size_limit(CSV_FIELD_SIZE_LIMIT)
    with path.open(newline="") as handle:
        rows = list(csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t"))
    width = 1 if summary else 6
    header = ["Status"] if summary else FEATURE_HEADER
    if not rows or rows[0][:width] != header or len(rows[0]) != width + 1:
        raise ValueError(f"Unexpected featureCounts header in {path}")
    result = {}
    for row in rows[1:]:
        if len(row) != width + 1 or not re.fullmatch(r"[0-9]+", row[-1]):
            raise ValueError(f"Unexpected featureCounts row in {path}: {row}")
        key = tuple(row[:-1])
        if key in result:
            raise ValueError(f"Duplicate featureCounts row in {path}: {key}")
        result[key] = row[-1]
    if summary and ("Assigned",) not in result:
        raise ValueError(f"Missing Assigned result in {path}")
    return result


def no_feature_summary(sam):
    import pysam

    counts = {key: 0 for key in SUMMARY_STATUSES}
    with pysam.AlignmentFile(str(sam), "r") as handle:
        for read in handle.fetch(until_eof=True):
            if read.is_unmapped:
                counts["Unassigned_Unmapped"] += 1
            elif read.is_secondary or read.is_supplementary:
                counts["Unassigned_Secondary"] += 1
            else:
                counts["Unassigned_NoFeatures"] += 1
    return {(key,): str(value) for key, value in counts.items()}


def merge_column(matrix, values, path):
    if matrix is None:
        return {key: [value] for key, value in values.items()}
    if matrix.keys() != values.keys():
        raise ValueError(f"Inconsistent featureCounts rows across samples: {path}")
    for key in matrix:
        matrix[key].append(values[key])
    return matrix


def count_features(*, threads, annotation, gene_id, feature_type, output, sams, allow_overlap=False):
    if not sams:
        raise ValueError("At least one SAM input is required")
    if len(set(map(str, sams))) != len(sams):
        raise ValueError("Duplicate SAM inputs would create duplicate sample columns")
    output.parent.mkdir(parents=True, exist_ok=True)
    present = has_feature_type(annotation, feature_type)
    if not present:
        print(f"WARNING: No {feature_type} annotations in {annotation}; writing zero feature counts", file=sys.stderr)
    table, summary = None, None
    with tempfile.TemporaryDirectory(prefix=".featurecounts-", dir=output.parent) as scratch:
        scratch = Path(scratch)
        for index, sam in enumerate(sams):
            if present:
                single = scratch / f"sample_{index}.txt"
                command = ["featureCounts", "-T", str(threads), "-a", str(annotation),
                           "-s", "1", "-g", gene_id, "-t", feature_type, "-o", str(single)]
                if allow_overlap:
                    command.append("-O")
                command.append(str(sam))
                subprocess.run(command, check=True)
                values = read_single_table(single)
                statuses = read_single_table(single.with_suffix(".txt.summary"), summary=True)
            else:
                values, statuses = {}, no_feature_summary(sam)
            table = merge_column(table, values, sam)
            summary = merge_column(summary, statuses, sam)

        table_path, summary_path = scratch / "combined.txt", scratch / "combined.txt.summary"
        for path, header, rows in ((table_path, FEATURE_HEADER, table),
                                   (summary_path, ["Status"], summary)):
            with path.open("w", newline="") as handle:
                if path == table_path:
                    handle.write("# Program:count_features.py; featureCounts per sample; stranded=1\n")
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow([*header, *map(str, sams)])
                writer.writerows((*key, *values) for key, values in rows.items())
        os.replace(table_path, output)
        os.replace(summary_path, Path(str(output) + ".summary"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--gene-id", required=True)
    parser.add_argument("--feature-type", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allow-overlap", action="store_true")
    parser.add_argument("sams", type=Path, nargs="+")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    count_features(threads=args.threads, annotation=args.annotation, gene_id=args.gene_id,
                   feature_type=args.feature_type, output=args.output, sams=args.sams,
                   allow_overlap=args.allow_overlap)


if __name__ == "__main__":
    main()
