#!/usr/bin/env python3
"""Combine trimming, BBSplit, coverage, and featureCounts results by sample."""

import csv
import re
from pathlib import Path

from read_breakdown_csv import bbsplit_counts, cutadapt_counts, length_filter_counts, sample_from_log


# Keep the established column labels, including their historical spelling/spacing.
HEADER = (
    "Sample",
    " Raw Reads",
    " Reads Written",
    " Overall Allignment",
    " Percent Covered",
    " Mean Coverage",
    " Reads Assigned to features",
    " Features with non-zero reads",
    "Reads eligible for alignment (<=6000 bp)",
    "Reads over 6000 bp (not processed)",
)


def feature_counts(path):
    """Return assigned-read and nonzero-feature totals keyed by sample."""
    with path.open(newline="") as handle:
        rows = list(csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t"))
    if not rows or len(rows[0]) < 7 or rows[0][0] != "Geneid":
        raise ValueError(f"Unexpected featureCounts table format: {path}")
    samples = []
    for alignment in rows[0][6:]:
        match = re.fullmatch(r"BACTERIA_(.+)_\d+bp\.sam", Path(alignment).name)
        if not match or match.group(1) in samples:
            raise ValueError(f"Unexpected alignment name in {path}: {alignment}")
        samples.append(match.group(1))
    counts = {sample: [0, 0] for sample in samples}
    for row in rows[1:]:
        if len(row) != len(rows[0]):
            raise ValueError(f"Incomplete featureCounts row in {path}")
        for sample, count in zip(samples, row[6:]):
            if not re.fullmatch(r"[0-9]+", count):
                raise ValueError(f"Unexpected featureCounts read count in {path}: {count}")
            reads = int(count)
            counts[sample][0] += reads
            counts[sample][1] += int(reads != 0)
    return counts


def coverage_counts(path, aligned=None):
    """Keep the coverage values reported by samtools' histogram summary."""
    lines = path.read_text().splitlines()
    # samtools omits histogram output when there are no mapped reads.
    if not lines and aligned == 0:
        return "0", "0"
    if len(lines) < 6 or ":" not in lines[4] or ":" not in lines[5]:
        raise ValueError(f"Unexpected coverage report format: {path}")
    percent = lines[4].rsplit(":", 1)[1].replace(" ", "")
    mean = lines[5].rsplit(":", 1)[1].replace(" ", "")
    if not percent or not mean:
        raise ValueError(f"Missing coverage values in {path}")
    return percent, mean


def build_rows(output_dir):
    cutadapt_logs = sorted(output_dir.glob("*cutadapt_log.txt"), key=lambda path: path.name.casefold())
    if not cutadapt_logs:
        raise ValueError(f"No cutadapt logs found in {output_dir}")
    assigned = feature_counts(output_dir / "featurecounts_BACTERIA_summary.txt")
    alignment_dir = output_dir / "bbsplit_alignments/bacteria"
    rows = []
    seen = set()
    for log_path in cutadapt_logs:
        sample = sample_from_log(log_path)
        if sample in seen:
            raise ValueError(f"Duplicate cutadapt sample: {sample}")
        seen.add(sample)
        length = re.search(r"_(\d+)bp\.cutadapt_log", log_path.name).group(1)
        prefix = f"BACTERIA_{sample}_{length}bp"
        total, passing = cutadapt_counts(log_path)
        length_input, eligible, oversized = length_filter_counts(output_dir / f"{sample}_{length}bp.length_filter.tsv")
        decoy_input, decoy_aligned = bbsplit_counts(
            output_dir / f"bbsplit_alignments/decoy/{sample}_{length}bp.mapped_to_other_bugs.bbsplit.txt"
        )
        bacteria_input, aligned = bbsplit_counts(alignment_dir / f"{prefix}.bbsplit.txt")
        if length_input != passing or decoy_input != eligible or bacteria_input != eligible - decoy_aligned:
            raise ValueError(f"Stage totals are inconsistent for sample {sample}")
        if sample not in assigned:
            raise ValueError(f"Missing featureCounts result for sample {sample}")
        percent, mean = coverage_counts(alignment_dir / f"{prefix}_coverage.txt", aligned)
        rows.append((sample, total, passing, aligned, percent, mean, *assigned[sample], eligible, oversized))
    return rows


def main():
    output_dir = Path.cwd()
    rows = build_rows(output_dir)
    # Retain the two coverage text files produced by the original converter.
    (output_dir / "BACTERIA_percentCoverage.txt").write_text("".join(f"{row[4]}\n" for row in rows))
    (output_dir / "BACTERIA_meanCoverage.txt").write_text("".join(f"{row[5]}\n" for row in rows))
    with (output_dir / f"BACTERIA_{output_dir.name}.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(HEADER)
        writer.writerows(rows)


if __name__ == "__main__":
    main()
