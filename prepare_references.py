#!/usr/bin/env python3
"""Discover reference folders and prepare separate RNA and genomic DNA inputs.

Reference files are read recursively; generated files always live in the run's
output directory. RNA sequence IDs are namespaced because transcript coordinates
must never be confused with genomic coordinates when SAM files are combined.
"""

import argparse
import gzip
import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import quote, unquote


FASTA_EXTENSIONS = (".fa", ".fasta", ".fna")
ANNOTATION_EXTENSIONS = (".gff", ".gff3", ".gbff", ".gb", ".gbk")


def open_text(path):
    return gzip.open(path, "rt") if path.name.lower().endswith(".gz") else path.open()


def extension(path):
    return Path(path.stem if path.suffix.lower() == ".gz" else path.name).suffix.lower()


def stem(path):
    name = path.name[:-3] if path.name.lower().endswith(".gz") else path.name
    return str(Path(name).with_suffix("")).lower()


def molecule(path):
    """NCBI extracted CDS are transcript targets, not whole genomic DNA."""
    name = path.name.lower()
    if "cds_from_genomic" in name or "rna_from_genomic" in name:
        return "rna"
    if re.search(r"(?:^|[._-])rna(?:[._-]|$)", name):
        return "rna"
    if any(re.search(r"(?:^|[._-])rna(?:[._-]|$)", p.lower()) for p in path.parts[:-1]):
        return "rna"
    return "dna"


def discover(value):
    path = Path(value).expanduser().resolve()
    if path.is_dir():
        files = sorted((p for p in path.rglob("*") if p.is_file()), key=lambda p: str(p).casefold())
    else:
        # Older configs accepted a shared filename root. Keep these runs usable
        # while new configs can point at the top-level organism folder.
        files = sorted((p for p in path.parent.glob(path.name + ".*") if p.is_file()),
                       key=lambda p: str(p).casefold())
        if path.is_file():
            files.insert(0, path)
    fasta = [p for p in files if extension(p) in FASTA_EXTENSIONS]
    annotations = [p for p in files if extension(p) in ANNOTATION_EXTENSIONS]
    if not fasta:
        raise ValueError(f"No .fa/.fasta/.fna[.gz] reference files found in {value}")
    return fasta, annotations


def sequence_id(header):
    if not header.strip():
        raise ValueError("Empty reference FASTA header")
    ident = header.split()[0]
    # Older NCBI genomic FASTAs use gi|...|ref|ACCESSION|. Current files use
    # the accession directly; use the same identifier for their annotations.
    match = re.search(r"(?:^|\|)(?:ref|gb|emb|dbj)\|([^|]+)", ident)
    return match.group(1) if match else ident


def attributes(text):
    result = {}
    for item in text.strip().strip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key.strip()] = unquote(value.strip())
        elif " " in item.strip():
            key, value = item.strip().split(None, 1)
            result[key] = value.strip('"')
    return result


def is_rna(feature_type):
    low = feature_type.lower()
    return "rna" in low or low in {"transcript", "primary_transcript", "ribozyme"}


def relevant(feature_type):
    # Exon hierarchy distinguishes mature ncRNA from introns on genomic DNA.
    return feature_type.lower() in {"gene", "cds", "exon"} or is_rna(feature_type)


@dataclass
class Feature:
    ref: str
    type: str
    start: int
    end: int
    strand: str
    attrs: dict


def gff_features(path):
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Invalid GFF row at {path}:{number}")
            if not relevant(fields[2]):
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as exc:
                raise ValueError(f"Invalid GFF coordinates at {path}:{number}") from exc
            yield Feature(fields[0], fields[2], start, end, fields[6], attributes(fields[8]))


def genbank_features(path):
    try:
        from Bio import SeqIO
    except ImportError as exc:
        raise ValueError("Biopython is required to read GBFF/GenBank references; install the pipeline environment") from exc
    with open_text(path) as handle:
        records = 0
        for record in SeqIO.parse(handle, "genbank"):
            records += 1
            if "rna" in record.annotations.get("molecule_type", "").lower():
                # Some RefSeq RNA records have only a source/gene (for example
                # pseudogene transcripts). Their complete sequence still has
                # transcript coordinates and belongs to the other-RNA class.
                yield Feature(record.id, "transcript", 1, len(record.seq), "+", {"ID": record.id})
            for index, feature in enumerate(record.features, 1):
                if not relevant(feature.type) or feature.location is None:
                    continue
                attrs = {key: str(values[0]) for key, values in feature.qualifiers.items() if values}
                if feature.type == "CDS":
                    attrs.pop("translation", None)
                attrs.setdefault("ID", f"{record.id}_{feature.type}_{index}")
                parts = list(feature.location.parts)
                for part in parts:
                    strand = "-" if part.strand == -1 else "+" if part.strand == 1 else "."
                    yield Feature(record.id, feature.type, int(part.start) + 1, int(part.end), strand, dict(attrs))
        if not records:
            raise ValueError(f"No GenBank records found in {path}")


def read_annotation(path):
    return gff_features(path) if extension(path) in {".gff", ".gff3"} else genbank_features(path)


def choose_annotations(files):
    """Choose GFF over redundant GBFF of the same assembly/shard."""
    selected = {}
    for path in files:
        key = (path.parent, stem(path))
        rank = 0 if extension(path) in {".gff", ".gff3"} else 1
        previous = selected.get(key)
        if previous is None or rank < previous[0]:
            selected[key] = (rank, path)
        elif rank == previous[0]:
            raise ValueError(f"Multiple annotations for {path.parent / stem(path)}: {previous[1]} and {path}")
    return [entry[1] for entry in selected.values()]


def write_feature(handle, feature, *, ref=None):
    attrs = dict(feature.attrs)
    gene_id = next((attrs[key] for key in ("gene_id", "locus", "locus_tag", "gene", "ID", "Name", "Parent")
                    if attrs.get(key)), None)
    if not gene_id:
        gene_id = f"{feature.ref}_{feature.type}_{feature.start}_{feature.end}"
    attrs["gene_id"] = gene_id
    encoded = ";".join(f"{quote(str(key), safe='_.:-')}={quote(str(value), safe='_.:-,')}"
                       for key, value in attrs.items())
    handle.write("\t".join((ref or feature.ref, "reference", feature.type, str(feature.start),
                            str(feature.end), ".", feature.strand, "0" if feature.type == "CDS" else ".", encoded)) + "\n")


def combine_fasta(paths, destination, *, rna=False):
    """Stream sequences, recording lengths/headers and rejecting duplicate conflicts.

    Keeping only a line plus a wrapping buffer avoids materializing complete
    chromosomes several times while concatenating large host assemblies.
    """
    lengths, headers, digests = {}, {}, {}
    with destination.open("w") as output:
        for path in paths:
            count, header, pending = 0, None, ""

            def finish_record():
                if header is None:
                    return
                if not length:
                    raise ValueError(f"Empty reference sequence {ident} in {path}")
                digest = hasher.digest()
                if ident in lengths:
                    if digests[ident] != digest:
                        raise ValueError(f"Conflicting reference sequences share ID {ident} in {path}")
                    return
                if pending:
                    output.write(pending + "\n")
                lengths[ident], headers[ident], digests[ident] = length, (header, path), digest

            with open_text(path) as source:
                for line in source:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        finish_record()
                        count += 1
                        header, pending = line[1:], ""
                        ident = sequence_id(header)
                        hasher, length = hashlib.sha256(), 0
                        if ident not in lengths:
                            output.write(f">{'RNA__' if rna else ''}{ident}\n")
                        continue
                    if header is None:
                        raise ValueError(f"Sequence before FASTA header in {path}")
                    hasher.update(line.upper().encode())
                    length += len(line)
                    if ident not in lengths:
                        pending += line
                        wrapped = len(pending) // 80 * 80
                        for offset in range(0, wrapped, 80):
                            output.write(pending[offset:offset + 80] + "\n")
                        pending = pending[wrapped:]
                finish_record()
            if not count:
                raise ValueError(f"No reference sequences found in {path}")
    return lengths, headers


def validate_feature(feature, lengths, path):
    if feature.ref not in lengths:
        raise ValueError(f"Annotation reference {feature.ref} in {path} has no matching FASTA sequence")
    if feature.start < 1 or feature.end < feature.start or feature.end > lengths[feature.ref]:
        raise ValueError(f"Annotation coordinates outside FASTA sequence {feature.ref} in {path}: {feature.start}..{feature.end}")


def header_annotations(header):
    return dict(re.findall(r"\[([^=\]]+)=([^\]]*)\]", header))


def rna_type(header, path):
    attrs = header_annotations(header)
    if "cds_from_genomic" in path.name.lower() or "protein_id" in attrs:
        return "CDS"
    key = attrs.get("gbkey", "")
    product = attrs.get("product", "").lower()
    if key.lower() == "cds":
        return "CDS"
    if key.lower() == "rrna" or "ribosomal rna" in product:
        return "rRNA"
    if key.lower() == "trna" or "transfer rna" in product:
        return "tRNA"
    if key and is_rna(key):
        return key
    return None


def prepare_organism(value, output_dir, label):
    fastas, annotations = discover(value)
    dna_files = [p for p in fastas if molecule(p) == "dna"]
    rna_files = [p for p in fastas if molecule(p) == "rna"]
    if not dna_files:
        raise ValueError(f"No whole genomic DNA FASTA found for {label} in {value}")
    dna_fasta, dna_gff = output_dir / f"{label}.dna.fasta", output_dir / f"{label}.dna.gff3"
    dna_lengths, _ = combine_fasta(dna_files, dna_fasta)
    dna_annotations = choose_annotations([p for p in annotations if molecule(p) == "dna"])
    if not dna_annotations:
        raise ValueError(f"No genomic GFF/GFF3 or GBFF annotation found for {label} in {value}")
    # Extracted bacterial RNA headers can share locus/protein names with their
    # genomic annotation. Retain just these requested loci, not whole genomes.
    rna_headers = {}
    if rna_files:
        rna_fasta, rna_gff = output_dir / f"{label}.rna.fasta", output_dir / f"{label}.rna.gff3"
        rna_lengths, rna_headers = combine_fasta(rna_files, rna_fasta, rna=True)
    wanted = set()
    for header, _ in rna_headers.values():
        a = header_annotations(header)
        wanted.update(a[key] for key in ("locus_tag", "gene", "protein_id") if a.get(key))
    genomic_by_key = {}
    count = 0
    with dna_gff.open("w") as output:
        output.write("##gff-version 3\n")
        for path in dna_annotations:
            for feature in read_annotation(path):
                validate_feature(feature, dna_lengths, path)
                write_feature(output, feature)
                count += 1
                for key in ("locus_tag", "locus", "gene", "protein_id"):
                    name = feature.attrs.get(key)
                    if name in wanted:
                        genomic_by_key.setdefault(name, []).append(feature)
    if not count:
        raise ValueError(f"No CDS or RNA/gene annotation features found for {label} in {value}")
    result = {"dna": {"fasta": str(dna_fasta), "gff": str(dna_gff)}, "rna": None}
    if not rna_files:
        return result

    rna_annotations = choose_annotations([p for p in annotations if molecule(p) == "rna"])
    annotated = set()
    with rna_gff.open("w") as output:
        output.write("##gff-version 3\n")
        for path in rna_annotations:
            for feature in read_annotation(path):
                validate_feature(feature, rna_lengths, path)
                write_feature(output, feature, ref="RNA__" + feature.ref)
                if feature.type.lower() != "gene":
                    annotated.add(feature.ref)
        for ident, length in rna_lengths.items():
            header, path = rna_headers[ident]
            a = header_annotations(header)
            gene = next((a[key] for key in ("locus_tag", "gene", "protein_id") if a.get(key)), ident)
            # Transcript coordinates are always on the reference's forward
            # strand, including transcripts extracted from a minus-strand gene.
            write_feature(output, Feature(ident, "transcript", 1, length, "+", {"gene_id": gene}),
                          ref="RNA__" + ident)
            if ident in annotated:
                continue
            feature_type = rna_type(header, path)
            candidates = []
            for key in ("locus_tag", "gene", "protein_id"):
                candidates.extend(genomic_by_key.get(a.get(key), []))
            if not feature_type:
                types = {f.type for f in candidates if f.type.lower() != "gene"}
                if "CDS" in types:
                    feature_type = "CDS"
                else:
                    feature_type = next((t for t in sorted(types) if is_rna(t)), None)
            if not feature_type and "rna_from_genomic" in path.name.lower():
                # An NCBI RNA extraction is explicitly an RNA reference even
                # when the producer omitted product/gbkey qualifiers.
                feature_type = "misc_RNA"
            if feature_type:
                write_feature(output, Feature(ident, feature_type, 1, length, "+", {"gene_id": gene}),
                              ref="RNA__" + ident)
            else:
                raise ValueError(f"Cannot classify RNA reference {ident} in {path}: supply matching RNA GBFF/GFF or NCBI header qualifiers")
    result["rna"] = {"fasta": str(rna_fasta), "gff": str(rna_gff)}
    return result


def prepare_references(decoy_dir, bacteria_dir, host_dir, output_dir):
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    decoy_files, _ = discover(decoy_dir)
    decoy_fasta = output_dir / "decoy.fasta"
    combine_fasta(decoy_files, decoy_fasta)
    manifest = {
        "decoy": {"fasta": str(decoy_fasta)},
        "bacteria": prepare_organism(bacteria_dir, output_dir, "bacteria"),
        "host": prepare_organism(host_dir, output_dir, "host") if host_dir else None,
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--decoy-dir", required=True)
    parser.add_argument("--bacteria-dir", required=True)
    parser.add_argument("--host-dir", default="")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        prepare_references(args.decoy_dir, args.bacteria_dir, args.host_dir, args.output_dir)
    except (ValueError, OSError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
