#!/usr/bin/env python3
"""Create a per-sample read disposition table from pipeline reports."""

import argparse
import csv
import re
from pathlib import Path


ATTRIBUTES = (
    "Total reads",
    "Reads passing cutadapt length filter",
    "Reads eligible for alignment (<=6000 bp)",
    "Reads over 6000 bp (not processed)",
    "Decoy aligned",
    "Bacteria aligned",
    "Bacteria non-rRNA",
    "Bacteria rRNA",
    "Host aligned",
    "Host non-rRNA",
    "Host rRNA",
    "Leftover reads",
)


def report_number(report, pattern, path):
    match = re.search(pattern, report, re.MULTILINE)
    if not match:
        raise ValueError(f"Could not find the required count in {path}")
    return int(match.group(1).replace(",", ""))


def cutadapt_counts(path):
    report = path.read_text()
    # Cutadapt replaces its numeric summary with this message for empty input.
    if re.search(r"^No reads processed!$", report, re.MULTILINE):
        return 0, 0
    total = report_number(
        report, r"Total reads processed:\s*([\d,]+)", path
    )
    passing = report_number(
        report, r"Reads written \(passing filters\):\s*([\d,]+)", path
    )
    if passing > total:
        raise ValueError(f"Passing reads exceed total reads in {path}")
    return total, passing


def metric_counts(path, required):
    """Read a complete pipeline metric/count report without ignoring bad rows."""
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or rows[0] != ["metric", "count"]:
        raise ValueError(f"Unexpected metric report header: {path}")
    metrics = {}
    for row in rows[1:]:
        if len(row) != 2 or row[0] in metrics or not re.fullmatch(r"[0-9]+", row[1]):
            raise ValueError(f"Unexpected metric report row in {path}: {row}")
        metrics[row[0]] = int(row[1])
    if set(metrics) != set(required):
        raise ValueError(f"Unexpected metrics in {path}: expected {', '.join(required)}")
    return metrics


def bbsplit_counts(path):
    """Return input and aligned counts from the pipeline BBSplit TSV report."""
    metrics = metric_counts(path, ("input", "aligned", "unmapped"))
    if metrics["input"] != metrics["aligned"] + metrics["unmapped"]:
        raise ValueError(f"Inconsistent BBSplit counts in {path}")
    return metrics["input"], metrics["aligned"]


def length_filter_counts(path):
    """Return counts for the post-trimming 6000 bp alignment limit."""
    metrics = metric_counts(path, ("input", "eligible", "oversized", "max_read_length"))
    if metrics["max_read_length"] != 6000:
        raise ValueError(f"Unexpected maximum read length in {path}: expected 6000")
    if metrics["input"] != metrics["eligible"] + metrics["oversized"]:
        raise ValueError(f"Inconsistent length-filter counts in {path}")
    return metrics["input"], metrics["eligible"], metrics["oversized"]


def sample_from_log(path):
    match = re.fullmatch(r"(.+)_\d+bp\.cutadapt_log\.txt", path.name)
    if not match:
        raise ValueError(f"Unexpected cutadapt log name: {path}")
    return match.group(1)


def assigned_by_sample(path, prefix):
    """Read the Assigned row of a featureCounts .summary file."""
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or not rows[0] or rows[0][0] != "Status":
        raise ValueError(f"Unexpected featureCounts summary format: {path}")
    assigned = next((row for row in rows[1:] if row and row[0] == "Assigned"), None)
    if assigned is None:
        raise ValueError(f"No Assigned row in {path}")
    if len(assigned) != len(rows[0]):
        raise ValueError(f"Incomplete Assigned row in {path}")
    result = {}
    for sam, count in zip(rows[0][1:], assigned[1:]):
        name = Path(sam).name
        match = re.fullmatch(rf"{re.escape(prefix)}_(.+)_\d+bp\.sam", name)
        if not match:
            raise ValueError(f"Unexpected alignment name in {path}: {sam}")
        if match.group(1) in result or not re.fullmatch(r"[0-9]+", count):
            raise ValueError(f"Unexpected sample count in {path}: {sam}")
        result[match.group(1)] = int(count)
    return result


def percentage(count, total):
    return f"{count * 100 / total:.2f}%" if total else "0.00%"


def build_table(output_dir):
    cutadapt_logs = sorted(output_dir.glob("*cutadapt_log.txt"), key=lambda p: p.name.casefold())
    if not cutadapt_logs:
        raise ValueError(f"No cutadapt logs found in {output_dir}")

    bacteria_rrna_path = output_dir / "featurecounts_BACTERIA_rRNA.txt.summary"
    bacteria_rrna = assigned_by_sample(bacteria_rrna_path, "BACTERIA")
    host_rrna_path = output_dir / "bbsplit_alignments/host/featurecounts_HOST_rRNA.txt.summary"
    host_enabled = host_rrna_path.exists()
    host_rrna = assigned_by_sample(host_rrna_path, "HOST") if host_enabled else {}

    rows = []
    for cutadapt_path in cutadapt_logs:
        sample = sample_from_log(cutadapt_path)
        total, passing = cutadapt_counts(cutadapt_path)
        # Allow a configured minimum length other than the historical 22 bp default.
        length = re.search(r"_(\d+)bp\.cutadapt_log", cutadapt_path.name).group(1)
        length_path = output_dir / f"{sample}_{length}bp.length_filter.tsv"
        length_input, eligible, oversized = length_filter_counts(length_path)
        decoy_path = output_dir / f"bbsplit_alignments/decoy/{sample}_{length}bp.mapped_to_other_bugs.bbsplit.txt"
        bacteria_path = output_dir / f"bbsplit_alignments/bacteria/BACTERIA_{sample}_{length}bp.bbsplit.txt"
        decoy_input, decoy_aligned = bbsplit_counts(decoy_path)
        bacteria_input, bacteria_aligned = bbsplit_counts(bacteria_path)
        if length_input != passing or decoy_input != eligible or bacteria_input != eligible - decoy_aligned:
            raise ValueError(f"Stage totals are inconsistent for sample {sample}")

        host_path = output_dir / f"bbsplit_alignments/host/HOST_{sample}_{length}bp.bbsplit.txt"
        if host_enabled:
            host_input, host_aligned = bbsplit_counts(host_path)
            if host_input != bacteria_input - bacteria_aligned:
                raise ValueError(f"Host stage total is inconsistent for sample {sample}")
        else:
            host_aligned = 0

        if sample not in bacteria_rrna or (host_enabled and sample not in host_rrna):
            raise ValueError(f"Missing rRNA featureCounts result for sample {sample}")
        bacterial_rrna = bacteria_rrna[sample]
        host_rrna_count = host_rrna[sample] if host_enabled else 0
        if bacterial_rrna > bacteria_aligned or host_rrna_count > host_aligned:
            raise ValueError(f"rRNA count exceeds aligned count for sample {sample}")
        counts = {
            "Total reads": total,
            "Reads passing cutadapt length filter": passing,
            "Reads eligible for alignment (<=6000 bp)": eligible,
            "Reads over 6000 bp (not processed)": oversized,
            "Decoy aligned": decoy_aligned,
            "Bacteria aligned": bacteria_aligned,
            "Bacteria non-rRNA": bacteria_aligned - bacterial_rrna,
            "Bacteria rRNA": bacterial_rrna,
            "Host aligned": host_aligned,
            "Host non-rRNA": host_aligned - host_rrna_count,
            "Host rRNA": host_rrna_count,
            "Leftover reads": eligible - decoy_aligned - bacteria_aligned - host_aligned,
        }
        sources = {
            "Total reads": f"{cutadapt_path.name}: Total reads processed",
            "Reads passing cutadapt length filter": f"{cutadapt_path.name}: Reads written (passing filters)",
            "Reads eligible for alignment (<=6000 bp)": f"{length_path.name}: eligible",
            "Reads over 6000 bp (not processed)": f"{length_path.name}: oversized",
            "Decoy aligned": f"{decoy_path.relative_to(output_dir)}: pipeline primary-alignment count",
            "Bacteria aligned": f"{bacteria_path.relative_to(output_dir)}: pipeline primary-alignment count",
            "Bacteria rRNA": f"{bacteria_rrna_path.relative_to(output_dir)}: Assigned",
            "Host aligned": (f"{host_path.relative_to(output_dir)}: pipeline primary-alignment count" if host_enabled else "Host alignment not configured"),
            "Host rRNA": (f"{host_rrna_path.relative_to(output_dir)}: Assigned" if host_enabled else "Host alignment not configured"),
        }
        sources["Bacteria non-rRNA"] = "Calculated: Bacteria aligned minus Bacteria rRNA"
        sources["Host non-rRNA"] = "Calculated: Host aligned minus Host rRNA"
        sources["Leftover reads"] = "Calculated: eligible reads minus decoy, bacteria, and host aligned"
        rows.append((sample, counts, sources))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    destination = args.output or output_dir / f"BACTERIA_{output_dir.name}_read_breakdown.csv"
    header = ["Sample"]
    for attribute in ATTRIBUTES:
        header.extend((f"{attribute} (reads)", f"{attribute} (% of total reads)", f"{attribute} source"))
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for sample, counts, sources in build_table(output_dir):
            row = [sample]
            for attribute in ATTRIBUTES:
                row.extend((counts[attribute], percentage(counts[attribute], counts["Total reads"]), sources[attribute]))
            writer.writerow(row)


if __name__ == "__main__":
    main()
