#!/usr/bin/env python3
"""Discover and stage one Nanopore sample per FASTQ file."""

import gzip
import hashlib
import random
from pathlib import Path


FASTQ_SUFFIXES = (".fastq.gz", ".fastq", ".fq.gz", ".fq")


def fastq_stem(path):
    for suffix in FASTQ_SUFFIXES:
        if path.name.lower().endswith(suffix):
            return path.name[:-len(suffix)]
    raise ValueError(f"Not a FASTQ file: {path}")


def input_metadata(path):
    path = path.resolve()
    info = path.stat()
    return {"path": str(path), "size": info.st_size, "mtime_ns": info.st_mtime_ns}


def discover_samples(directory, pattern="*.fastq*"):
    directory = Path(directory).expanduser().resolve()
    if not directory.is_dir():
        raise ValueError(f"FASTQ_DIR is not a directory: {directory}")
    files = sorted((p for p in directory.glob(pattern) if p.is_file() and
                    p.name.lower().endswith(FASTQ_SUFFIXES)), key=lambda p: str(p).casefold())
    if not files:
        raise ValueError(f"No FASTQ files in {directory} matching {pattern}")
    if len({p.resolve() for p in files}) != len(files):
        raise ValueError("Several FASTQ paths refer to the same input file")
    result, seen = [], set()
    for index, path in enumerate(sorted(files, key=lambda path: fastq_stem(path).casefold())):
        name = fastq_stem(path)
        if not name or name in {".", ".."} or "\n" in name or "\r" in name:
            raise ValueError(f"Invalid sample name: {name!r}")
        if name.casefold() in seen:
            raise ValueError(f"Several inputs produce sample name {name!r}; rename them to avoid collisions")
        seen.add(name.casefold())
        result.append({"index": index, "sample": name, "reads": [input_metadata(path)]})
    return result


def check_inputs(sample):
    for expected in sample["reads"]:
        if input_metadata(Path(expected["path"])) != expected:
            raise ValueError(f"Input changed after submission: {expected['path']}")


def records(path):
    path = Path(path)
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    with opener(path, "rt") as source:
        number = 0
        while True:
            header = source.readline()
            if not header:
                break
            number += 1
            rest = [source.readline() for _ in range(3)]
            if any(line == "" for line in rest):
                raise ValueError(f"Incomplete FASTQ record {number} in {path}")
            sequence, plus, quality = [line.rstrip("\r\n") for line in rest]
            if not header.startswith("@") or not plus.startswith("+") or len(sequence) != len(quality):
                raise ValueError(f"Malformed FASTQ record {number} in {path}")
            yield (header.rstrip("\r\n"), sequence, plus, quality)


def stage_sample(sample, directory, *, sample_reads=False, sample_size=5000, seed="1"):
    """Link the original FASTQ or reservoir-sample intact Nanopore records."""
    check_inputs(sample)
    if len(sample["reads"]) != 1:
        raise ValueError("Each Nanopore sample requires exactly one FASTQ input")
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    path = Path(sample["reads"][0]["path"])
    suffix = ".fastq.gz" if sample_reads or path.name.lower().endswith(".gz") else ".fastq"
    destination = directory / (sample["sample"] + suffix)
    if destination == path:
        raise ValueError("Staging directory must not overwrite the original FASTQ")
    destination.unlink(missing_ok=True)
    if not sample_reads:
        destination.symlink_to(path)
        return destination
    if sample_size < 1:
        raise ValueError("SAMPLE_SIZE must be positive")
    fingerprint = hashlib.sha256(str(path).encode()).hexdigest()
    rng, reservoir = random.Random(f"{seed}:{fingerprint}"), []
    for count, record in enumerate(records(path)):
        if count < sample_size:
            reservoir.append(record)
        else:
            replacement = rng.randrange(count + 1)
            if replacement < sample_size:
                reservoir[replacement] = record
    with gzip.open(destination, "wt", compresslevel=1) as output:
        for record in reservoir:
            output.write("\n".join(record) + "\n")
    check_inputs(sample)
    return destination
