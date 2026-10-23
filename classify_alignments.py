#!/usr/bin/env python3
"""Count exclusive read classes across RNA-first and genomic alignments.

Each mapped primary alignment is counted once, without collapsing read names.
Only aligned CIGAR blocks on the annotated strand contribute feature overlaps.
The priority for overlapping classes is rRNA, CDS, then other RNA. Everything
else is intergenic/unannotated. RNA-reference transcript extents count as other
RNA outside their CDS; genomic gene/mRNA extents do not imply a CDS or RNA hit.
"""

import argparse
from bisect import bisect_left
from collections import defaultdict
import csv
import gzip
from pathlib import Path
from urllib.parse import unquote

import pysam


CATEGORIES = ("CDS", "rRNA", "other_RNA", "intergenic_unannotated")
PRIORITY = ("rRNA", "CDS", "other_RNA")


class Intervals:
    """Disjoint half-open intervals supporting logarithmic overlap queries."""

    def __init__(self, intervals):
        merged = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
            else:
                merged.append((start, end))
        self.starts = [start for start, _ in merged]
        self.ends = [end for _, end in merged]

    def overlaps(self, start, end):
        if end <= start:
            return False
        index = bisect_left(self.starts, end) - 1
        return index >= 0 and self.ends[index] > start


def feature_category(feature_type, rna_reference=False):
    low = feature_type.lower()
    if low == "cds":
        return "CDS"
    if low == "rrna":
        return "rRNA"
    if low in {"mrna", "transcript", "primary_transcript"}:
        return "other_RNA" if rna_reference else None
    if "rna" in low or low in {"ribozyme", "guide_rna"}:
        return "other_RNA"
    return None


def gff_rows(path):
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt") as source:
        for number, line in enumerate(source, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Expected nine GFF fields at {path}:{number}")
            yield number, fields


def feature_links(attributes):
    """Decode GFF IDs without splitting percent-encoded commas within an ID."""
    fields = dict(item.split("=", 1) for item in attributes.split(";") if "=" in item)
    ident = unquote(fields.get("ID", ""))
    parents = tuple(unquote(parent) for parent in fields.get("Parent", "").split(",") if parent)
    return ident, parents


def annotation_index(path, references, rna_reference=False):
    """Validate annotations and index genomic noncoding RNAs by their exons.

    Noncoding RNA features without exon descendants retain their own interval,
    as in bacterial annotations. Two streaming passes avoid retaining the many
    coding exon rows in a human annotation solely to find noncoding exons.
    """
    path = Path(path)
    intervals = defaultdict(list)
    noncoding = []
    parent_links = {}
    has_exons = False
    reference_names = {name: name for name in references}
    for number, fields in gff_rows(path):
        reference, _, kind, start, end, _, strand, _, attributes = fields
        if reference not in references:
            raise ValueError(f"GFF reference {reference!r} is absent from the alignment header: {path}:{number}")
        reference = reference_names[reference]
        try:
            start, end = int(start) - 1, int(end)
        except ValueError as exc:
            raise ValueError(f"Invalid GFF coordinates at {path}:{number}") from exc
        if start < 0 or end <= start or end > references[reference]:
            raise ValueError(f"GFF coordinates exceed reference bounds at {path}:{number}")
        if strand not in {"+", "-", ".", "?"}:
            raise ValueError(f"Invalid GFF strand at {path}:{number}")
        low = kind.lower()
        if not rna_reference and low not in {"exon", "cds"}:
            ident, parents = feature_links(attributes)
            if ident and parents:
                parent_links[ident] = parents
        else:
            ident = ""
        has_exons |= low == "exon"
        category = feature_category(kind, rna_reference)
        # Strandless features cannot support the pipeline's -s 1 assignment.
        if category and strand in {"+", "-"}:
            if not rna_reference and category in {"rRNA", "other_RNA"}:
                noncoding.append((reference, start, end, strand, category, ident))
            else:
                intervals[(reference, strand, category)].append((start, end))

    target_ids = {feature[-1] for feature in noncoding if feature[-1]}
    exons_by_target = defaultdict(list)
    ancestors_cache = {}

    def noncoding_ancestors(parent):
        if parent in ancestors_cache:
            return ancestors_cache[parent]
        seen, found, pending = set(), set(), [parent]
        while pending:
            ident = pending.pop()
            if ident in seen:
                continue
            seen.add(ident)
            if ident in ancestors_cache:
                found.update(ancestors_cache[ident])
                continue
            if ident in target_ids:
                found.add(ident)
            pending.extend(parent_links.get(ident, ()))
        ancestors_cache[parent] = frozenset(found)
        return ancestors_cache[parent]

    if target_ids and has_exons:
        for _, fields in gff_rows(path):
            if fields[2].lower() != "exon":
                continue
            _, parents = feature_links(fields[8])
            targets = set()
            for parent in parents:
                targets.update(noncoding_ancestors(parent))
            if not targets:
                continue
            exon = (reference_names[fields[0]], int(fields[3]) - 1, int(fields[4]), fields[6])
            for target in targets:
                exons_by_target[target].append(exon)
    for reference, start, end, strand, category, ident in noncoding:
        if ident in exons_by_target:
            for exon_ref, exon_start, exon_end, exon_strand in exons_by_target[ident]:
                if exon_ref != reference or exon_strand not in {strand, ".", "?"}:
                    continue
                overlap_start, overlap_end = max(start, exon_start), min(end, exon_end)
                if overlap_end > overlap_start:
                    intervals[(reference, strand, category)].append((overlap_start, overlap_end))
        else:
            intervals[(reference, strand, category)].append((start, end))
    return {key: Intervals(values) for key, values in intervals.items()}


def classify_file(sam, gff, rna_reference=False):
    counts = {category: 0 for category in CATEGORIES}
    with pysam.AlignmentFile(str(sam), "r") as alignments:
        lengths = dict(zip(alignments.references, alignments.lengths))
        index = annotation_index(gff, lengths, rna_reference)
        for read in alignments.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.reference_start < 0 or read.reference_end is None or read.reference_end > lengths[read.reference_name]:
                raise ValueError(f"Alignment exceeds reference bounds in {sam}: {read.query_name}")
            strand = "-" if read.is_reverse else "+"
            blocks = read.get_blocks()
            category = "intergenic_unannotated"
            for candidate in PRIORITY:
                regions = index.get((read.reference_name, strand, candidate))
                if regions and any(regions.overlaps(start, end) for start, end in blocks):
                    category = candidate
                    break
            counts[category] += 1
    return counts


def classify(dna_sam, dna_gff, rna_sam=None, rna_gff=None):
    if (rna_sam is None) != (rna_gff is None):
        raise ValueError("--rna-sam and --rna-gff must be supplied together")
    counts = classify_file(dna_sam, dna_gff)
    if rna_sam is not None:
        rna_counts = classify_file(rna_sam, rna_gff, rna_reference=True)
        for category in CATEGORIES:
            counts[category] += rna_counts[category]
    return {"aligned": sum(counts.values()), **counts}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rna-sam", type=Path)
    parser.add_argument("--rna-gff", type=Path)
    parser.add_argument("--dna-sam", required=True, type=Path)
    parser.add_argument("--dna-gff", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if (args.rna_sam is None) != (args.rna_gff is None):
        parser.error("--rna-sam and --rna-gff must be supplied together")
    counts = classify(args.dna_sam, args.dna_gff, args.rna_sam, args.rna_gff)
    with args.output.open("w", newline="") as destination:
        writer = csv.writer(destination, delimiter="\t")
        writer.writerow(("metric", "count"))
        writer.writerows(counts.items())


if __name__ == "__main__":
    main()
