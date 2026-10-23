#!/usr/bin/env python3
"""Combine only the completed samples in a submitted run's manifest."""

import argparse
import csv
import json
import os
import tempfile
from pathlib import Path

from bacteria_with_decoys_csvConversion import HEADER, build_rows
from count_features import FEATURE_HEADER, merge_column, read_single_table
from read_breakdown_csv import ATTRIBUTES, build_table, percentage


def host_enabled(config):
    return bool(config.get("HOST_REFERENCE_DIR", config.get("HOST_BBSPLIT_REFERENCE", "")))


def combine_feature_tables(sources, destination):
    table, summary = None, None
    alignments = []
    for source in sources:
        with source.open(newline="") as handle:
            header = next(csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t"))
        if len(header) != 7:
            raise ValueError(f"Expected one sample in featureCounts table: {source}")
        alignments.append(header[6])
        table = merge_column(table, read_single_table(source), source)
        summary = merge_column(summary, read_single_table(Path(str(source) + ".summary"), summary=True), source)
    for path, header, rows in ((destination, FEATURE_HEADER, table),
                               (Path(str(destination) + ".summary"), ["Status"], summary)):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as handle:
            if path == destination:
                handle.write("# Program:aggregate_samples.py; all manifest samples\n")
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([*header, *alignments])
            writer.writerows((*key, *values) for key, values in rows.items())
    with destination.open(newline="") as source, destination.with_suffix(".csv").open("w", newline="") as output:
        csv.writer(output).writerows(csv.reader(source, delimiter="\t"))


def aggregate(manifest_path):
    manifest = json.loads(Path(manifest_path).read_text())
    root = Path(manifest["run_dir"])
    (root / "aggregation.completed.json").unlink(missing_ok=True)
    samples = manifest["samples"]
    if not samples:
        raise ValueError("Run manifest contains no samples")
    enabled = host_enabled(manifest["config"])
    main_rows, breakdown_rows, sample_dirs = [], [], []
    for sample in samples:
        directory = root / "samples" / sample["sample"]
        marker = directory / "completed.json"
        if not marker.is_file():
            raise ValueError(f"Sample {sample['sample']} has not completed; final summaries will not be written")
        status = json.loads(marker.read_text())
        expected = {"run_id": manifest["run_id"], "sample": sample["sample"],
                    "config_digest": manifest["config_digest"], "reads": sample["reads"]}
        if any(status.get(key) != value for key, value in expected.items()):
            raise ValueError(f"Completion marker does not match this run: {marker}")
        rows = build_rows(directory)
        breakdown = build_table(directory, host_enabled=enabled)
        if len(rows) != 1 or len(breakdown) != 1 or rows[0][0] != sample["sample"] or breakdown[0][0] != sample["sample"]:
            raise ValueError(f"Unexpected samples in output directory: {directory}")
        main_rows.extend(rows)
        name, counts, sources = breakdown[0]
        # Preserve usable provenance after moving a sample's row to the run CSV.
        sources = {key: (f"samples/{name}/{value}" if
                        (directory / value.split(": ", 1)[0]).is_file() else value)
                   for key, value in sources.items()}
        breakdown_rows.append((name, counts, sources))
        sample_dirs.append(directory)

    # Validate every sample before publishing any of the final tables.
    with tempfile.TemporaryDirectory(prefix=".aggregate-", dir=root) as temporary:
        stage = Path(temporary)
        prefix = f"BACTERIA_{root.name}"
        with (stage / f"{prefix}.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(HEADER)
            writer.writerows(main_rows)
        with (stage / f"{prefix}_read_breakdown.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            header = ["Sample"]
            for attribute in ATTRIBUTES:
                header.extend((f"{attribute} (reads)", f"{attribute} (% of total reads)", f"{attribute} source"))
            writer.writerow(header)
            for name, counts, sources in breakdown_rows:
                row = [name]
                for attribute in ATTRIBUTES:
                    row.extend((counts[attribute], percentage(counts[attribute], counts["Total reads"]), sources[attribute]))
                writer.writerow(row)
        (stage / "BACTERIA_percentCoverage.txt").write_text("".join(f"{row[4]}\n" for row in main_rows))
        (stage / "BACTERIA_meanCoverage.txt").write_text("".join(f"{row[5]}\n" for row in main_rows))
        for organism in ("BACTERIA", "HOST") if enabled else ("BACTERIA",):
            relative = Path() if organism == "BACTERIA" else Path("host")
            source_relative = Path() if organism == "BACTERIA" else Path("bbsplit_alignments/host")
            for kind in ("summary", "rRNA"):
                filename = f"featurecounts_{organism}_{kind}.txt"
                combine_feature_tables([directory / source_relative / filename for directory in sample_dirs],
                                       stage / relative / filename)
        files = [path for path in stage.rglob("*") if path.is_file()]
        for path in files:
            destination = root / path.relative_to(stage)
            destination.parent.mkdir(parents=True, exist_ok=True)
            os.replace(path, destination)
    completion = root / "aggregation.completed.json"
    temporary = completion.with_suffix(".tmp")
    temporary.write_text(json.dumps({"run_id": manifest["run_id"], "samples": len(samples),
                                     "config_digest": manifest["config_digest"]}, indent=2) + "\n")
    os.replace(temporary, completion)
    print(f"Aggregated {len(samples)} samples into {root}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    aggregate(args.manifest)


if __name__ == "__main__":
    main()
