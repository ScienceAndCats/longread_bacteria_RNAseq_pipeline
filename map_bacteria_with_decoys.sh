#!/usr/bin/env bash
set -euo pipefail

CONFIG_FILE="${1:-config.env}"
if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "Config file not found: $CONFIG_FILE" >&2
  echo "Usage: bash $0 [config.env]" >&2
  exit 1
fi

# shellcheck source=/dev/null
source "$CONFIG_FILE"

: "${FASTQ_DIR:=.}"
: "${FASTQ_GLOB:=*.fastq}"
: "${OUTPUT_DIR:=.}"
: "${SAMPLE_READS:=false}"
: "${SAMPLE_SIZE:=5000}"
: "${SAMPLE_SEED:=1}"
: "${DECOY_MINIMAP2_REFERENCE:?Set DECOY_MINIMAP2_REFERENCE in $CONFIG_FILE}"
: "${BACTERIA_MINIMAP2_REFERENCE:?Set BACTERIA_MINIMAP2_REFERENCE in $CONFIG_FILE}"
: "${HOST_MINIMAP2_REFERENCE:=}"
: "${MINIMAP2_PRESET:=map-ont}"
: "${ADAPTER_NANOPORE:=TTTCTGTTGGTGCTGATATTGC}"
: "${MIN_READ_LENGTH:=22}"
: "${THREADS:=16}"
: "${GENE_POSITION_BINS:=100}"
: "${METAGENE_MIN_FEATURE_READS:=10}"
: "${CSV_CONVERSION_SCRIPT:=bacteria_with_decoys_csvConversion.py}"

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be a positive integer" >&2; exit 1; }
[[ "$GENE_POSITION_BINS" =~ ^[1-9][0-9]*$ ]] || { echo "GENE_POSITION_BINS must be a positive integer" >&2; exit 1; }
[[ "$METAGENE_MIN_FEATURE_READS" =~ ^[1-9][0-9]*$ ]] || { echo "METAGENE_MIN_FEATURE_READS must be a positive integer" >&2; exit 1; }
case "$MINIMAP2_PRESET" in
  map-ont|map-hifi|map-pb) ;;
  *) echo "MINIMAP2_PRESET must be a long-read preset: map-ont, map-hifi, or map-pb" >&2; exit 1 ;;
esac

shopt -s nullglob

# Locate the single annotation whose filename starts with the reference
# basename and a .gff suffix (for example, reference.gff or reference.gff3).
find_gff_annotation() {
  local reference_basename="$1"
  local reference_label="$2"
  local candidates=("${reference_basename}".gff*)

  if (( ${#candidates[@]} == 0 )); then
    echo "ERROR: $reference_label annotation not found; expected ${reference_basename}.gff*" >&2
    return 1
  fi
  if (( ${#candidates[@]} > 1 )); then
    echo "ERROR: Multiple $reference_label annotations match ${reference_basename}.gff*: ${candidates[*]}" >&2
    return 1
  fi
  printf '%s\n' "${candidates[0]}"
}

# Prefer the locus attribute for featureCounts gene IDs, falling back first to
# locus_tag and then to gene when the preferred attributes do not occur.
find_gene_identifier_attribute() {
  local annotation_file="$1"

  if [[ "$annotation_file" == *.gz ]]; then
    gzip -cd -- "$annotation_file"
  else
    cat -- "$annotation_file"
  fi | awk -F '\t' '
    !/^#/ && NF >= 9 {
      attribute_count = split($9, attributes, ";")
      for (i = 1; i <= attribute_count; i++) {
        sub(/^[[:space:]]+/, "", attributes[i])
        split(attributes[i], parts, /[=[:space:]]/)
        if (parts[1] == "locus") has_locus = 1
        if (parts[1] == "locus_tag") has_locus_tag = 1
        if (parts[1] == "gene") has_gene = 1
      }
    }
    END {
      if (has_locus) print "locus"
      else if (has_locus_tag) print "locus_tag"
      else if (has_gene) print "gene"
      else exit 1
    }
  '
}

BACTERIA_GFF=$(find_gff_annotation "$BACTERIA_MINIMAP2_REFERENCE" "bacterial")
if ! BACTERIA_GENE_ATTRIBUTE=$(find_gene_identifier_attribute "$BACTERIA_GFF"); then
  echo "ERROR: bacterial annotation has none of the locus, locus_tag, or gene attributes: $BACTERIA_GFF" >&2
  exit 1
fi
HOST_GFF=""
HOST_GENE_ATTRIBUTE=""
if [[ -n "$HOST_MINIMAP2_REFERENCE" ]]; then
  HOST_GFF=$(find_gff_annotation "$HOST_MINIMAP2_REFERENCE" "host")
  if ! HOST_GENE_ATTRIBUTE=$(find_gene_identifier_attribute "$HOST_GFF"); then
    echo "ERROR: host annotation has none of the locus, locus_tag, or gene attributes: $HOST_GFF" >&2
    exit 1
  fi
fi

ALIGNMENT_DIR="$OUTPUT_DIR/minimap2_alignments"
DECOY_ALIGNMENT_DIR="$ALIGNMENT_DIR/decoy"
BACTERIA_ALIGNMENT_DIR="$ALIGNMENT_DIR/bacteria"
HOST_ALIGNMENT_DIR="$ALIGNMENT_DIR/host"
mkdir -p "$OUTPUT_DIR" "$DECOY_ALIGNMENT_DIR" "$BACTERIA_ALIGNMENT_DIR"
if [[ -n "$HOST_MINIMAP2_REFERENCE" ]]; then
  mkdir -p "$HOST_ALIGNMENT_DIR"
fi
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "$CSV_CONVERSION_SCRIPT" != /* ]]; then
  CSV_CONVERSION_SCRIPT="$SCRIPT_DIR/$CSV_CONVERSION_SCRIPT"
fi

# Reuse a minimap2 index or build one from a FASTA beside the configured
# basename. The basename convention also keeps annotation discovery unchanged.
prepare_minimap2_index() {
  local index_basename="$1"
  local reference_label="$2"
  local threads="$3"
  local fasta_extension fasta_file
  if [[ -f "${index_basename}.mmi" ]]; then
    echo "using existing $reference_label minimap2 index: ${index_basename}.mmi"
    return 0
  fi

  fasta_file=""
  for fasta_extension in fa fasta fna fa.gz fasta.gz fna.gz; do
    if [[ -f "${index_basename}.${fasta_extension}" ]]; then
      fasta_file="${index_basename}.${fasta_extension}"
      break
    fi
  done

  if [[ -z "$fasta_file" ]]; then
    echo "ERROR: No reference FASTA found; expected ${index_basename}.{fa,fasta,fna}[.gz]." >&2
    return 1
  fi

  mkdir -p "$(dirname "$index_basename")"
  echo "building $reference_label minimap2 index from $fasta_file ..."
  minimap2 -t "$threads" -d "${index_basename}.mmi" "$fasta_file"
}

prepare_minimap2_index "$DECOY_MINIMAP2_REFERENCE" "decoy" "$THREADS"
prepare_minimap2_index "$BACTERIA_MINIMAP2_REFERENCE" "bacterial" "$THREADS"
if [[ -n "$HOST_MINIMAP2_REFERENCE" ]]; then
  prepare_minimap2_index "$HOST_MINIMAP2_REFERENCE" "host" "$THREADS"
fi

matched_fastqs=("$FASTQ_DIR"/$FASTQ_GLOB)
if (( ${#matched_fastqs[@]} == 0 )); then
  echo "No FASTQ files found in $FASTQ_DIR matching $FASTQ_GLOB" >&2
  exit 1
fi

fastq_files=("${matched_fastqs[@]}")

# Randomly subsample complete FASTQ records with reservoir sampling. This keeps
# memory use bounded by SAMPLE_SIZE even for very large input files.
sample_fastq() {
  local input_file="$1"
  local output_file="$2"
  python3 - "$input_file" "$output_file" "$SAMPLE_SIZE" "$SAMPLE_SEED" <<'PY'
import gzip
import hashlib
import random
import sys

input_path, output_path, sample_size, seed = sys.argv[1:]
sample_size = int(sample_size)
rng = random.Random(f"{seed}:{hashlib.sha256(input_path.encode()).hexdigest()}")
opener = gzip.open if input_path.endswith(".gz") else open
reservoir = []

with opener(input_path, "rt") as source:
    read_count = 0
    while True:
        record = [source.readline() for _ in range(4)]
        if not record[0]:
            break
        if any(line == "" for line in record):
            raise SystemExit(f"Incomplete FASTQ record in {input_path}")
        if read_count < sample_size:
            reservoir.append(record)
        else:
            replacement = rng.randrange(read_count + 1)
            if replacement < sample_size:
                reservoir[replacement] = record
        read_count += 1

with open(output_path, "w") as destination:
    for record in reservoir:
        destination.writelines(record)
PY
}

case "${SAMPLE_READS,,}" in
  true)
    [[ "$SAMPLE_SIZE" =~ ^[1-9][0-9]*$ ]] || { echo "SAMPLE_SIZE must be a positive integer" >&2; exit 1; }
    sample_dir="$OUTPUT_DIR/.bacteria_sampled_fastq"
    rm -rf "$sample_dir"
    mkdir -p "$sample_dir"
    sampled_fastqs=()
    for input_fastq in "${fastq_files[@]}"; do
      sampled_fastq="$sample_dir/$(basename "${input_fastq%.gz}")"
      echo "sampling up to $SAMPLE_SIZE reads from $(basename "$input_fastq") ..."
      sample_fastq "$input_fastq" "$sampled_fastq"
      sampled_fastqs+=("$sampled_fastq")
    done
    fastq_files=("${sampled_fastqs[@]}")
    ;;
  false) ;;
  *) echo "SAMPLE_READS must be true or false" >&2; exit 1 ;;
esac

# Trim the configured long-read adapter and discard reads below the minimum length.
trimmed_fastqs=(); sample_names=()
for input_fastq in "${fastq_files[@]}"; do
  filename=$(basename "$input_fastq"); filename=${filename%.gz}; sample=${filename%.fastq}
  trimmed_fastq="$OUTPUT_DIR/${sample}_${MIN_READ_LENGTH}bp.trim.fastq"
  cutadapt -m "$MIN_READ_LENGTH" -j "$THREADS" -a "$ADAPTER_NANOPORE" \
    -o "$trimmed_fastq" "$input_fastq" \
    > "$OUTPUT_DIR/${sample}_${MIN_READ_LENGTH}bp.cutadapt_log.txt"
  trimmed_fastqs+=("$trimmed_fastq"); sample_names+=("$sample")
done

# Count complete FASTQ records (plain or gzip-compressed).
fastq_read_count() {
  if [[ "$1" == *.gz ]]; then gzip -cd -- "$1"; else cat -- "$1"; fi |
    awk 'END { if (NR % 4) exit 1; print NR / 4 }'
}

# Align reads with the nanopore-aware minimap2 preset, retain only primary
# alignments in SAM, and write unmapped reads for the next classification stage.
# The small TSV report is intentionally stable input for both CSV generators.
minimap2_stage() {
  local reference="$1" sam="$2" report="$3" unmapped="$4" input="$5"
  local input_count unmapped_count aligned_count raw_log="${report%.txt}.minimap2.stderr.txt"
  input_count=$(fastq_read_count "$input")
  minimap2 -t "$THREADS" -ax "$MINIMAP2_PRESET" --secondary=no "${reference}.mmi" "$input" \
    2> "$raw_log" | samtools view -h -F 2304 -o "$sam" -
  samtools fastq -@ "$THREADS" -f 4 -0 "$unmapped" -s /dev/null -n "$sam" 2>> "$raw_log"
  unmapped_count=$(fastq_read_count "$unmapped")
  aligned_count=$((input_count - unmapped_count))
  printf 'metric\tcount\ninput\t%s\naligned\t%s\nunmapped\t%s\n' \
    "$input_count" "$aligned_count" "$unmapped_count" > "$report"
}

# Map to bacterial decoys and pass only primary-unmapped reads onward.
decoy_unmapped_fastqs=()
for idx in "${!trimmed_fastqs[@]}"; do
  i="${trimmed_fastqs[$idx]}"; i_basename="${sample_names[$idx]}_${MIN_READ_LENGTH}bp"
  echo "mapping to bacterial decoys with minimap2: $i_basename ..."
  unmapped="$DECOY_ALIGNMENT_DIR/${i_basename}_unmapped_to_other_bugs.fastq.gz"
  minimap2_stage "$DECOY_MINIMAP2_REFERENCE" \
    "$DECOY_ALIGNMENT_DIR/${i_basename}.mapped_to_other_bugs.sam" \
    "$DECOY_ALIGNMENT_DIR/${i_basename}.mapped_to_other_bugs.minimap2.txt" \
    "$unmapped" "$i"
  decoy_unmapped_fastqs+=("$unmapped")
done

# Map decoy-unmapped reads to the target bacterium.
host_input_fastqs=()
for idx in "${!decoy_unmapped_fastqs[@]}"; do
  i="${decoy_unmapped_fastqs[$idx]}"; i_basename="${sample_names[$idx]}_${MIN_READ_LENGTH}bp"
  echo "mapping to the bacterial reference with minimap2: $i_basename ..."
  host_input="$BACTERIA_ALIGNMENT_DIR/${i_basename}_unmapped_to_bacteria.fastq.gz"
  minimap2_stage "$BACTERIA_MINIMAP2_REFERENCE" \
    "$BACTERIA_ALIGNMENT_DIR/BACTERIA_${i_basename}.sam" \
    "$BACTERIA_ALIGNMENT_DIR/BACTERIA_${i_basename}.minimap2.txt" \
    "$host_input" "$i"
  host_input_fastqs+=("$host_input")
done

BACTERIA_sam_filenames=("$BACTERIA_ALIGNMENT_DIR"/BACTERIA*.sam)

# Optionally classify reads that mapped to neither decoys nor bacteria as host.
if [[ -n "$HOST_MINIMAP2_REFERENCE" ]]; then
  for idx in "${!host_input_fastqs[@]}"; do
    i="${host_input_fastqs[$idx]}"; i_basename="${sample_names[$idx]}_${MIN_READ_LENGTH}bp"
    echo "mapping to the host reference with minimap2: $i_basename ..."
    minimap2_stage "$HOST_MINIMAP2_REFERENCE" \
      "$HOST_ALIGNMENT_DIR/HOST_${i_basename}.sam" \
      "$HOST_ALIGNMENT_DIR/HOST_${i_basename}.minimap2.txt" \
      "$HOST_ALIGNMENT_DIR/${i_basename}_unmapped_to_host.fastq.gz" "$i"
  done
fi

# Create featureCounts summary file.
# Assign mapped reads to bacterial genes, retaining the pipeline's stranded and
# overlapping-feature behavior.
featureCounts -T "$THREADS" -a "$BACTERIA_GFF" -O -s 1 -g "$BACTERIA_GENE_ATTRIBUTE" -t CDS \
  -o "$OUTPUT_DIR/featurecounts_BACTERIA_summary.txt" "${BACTERIA_sam_filenames[@]}"
sed 's/\t/,/g' "$OUTPUT_DIR/featurecounts_BACTERIA_summary.txt" > "$OUTPUT_DIR/featurecounts_BACTERIA_summary.csv"
# Count alignments assigned to an rRNA annotation separately.  The generated
# .summary file supplies the non-duplicated Assigned total used by the read
# disposition CSV; all other bacterial alignments are reported as non-rRNA.
featureCounts -T "$THREADS" -a "$BACTERIA_GFF" -s 1 -g "$BACTERIA_GENE_ATTRIBUTE" -t rRNA \
  -o "$OUTPUT_DIR/featurecounts_BACTERIA_rRNA.txt" "${BACTERIA_sam_filenames[@]}"

# When a host reference is configured, independently assign host alignments to
# its CDS features and keep the results alongside the host minimap2 outputs.
if [[ -n "$HOST_GFF" ]]; then
  HOST_sam_filenames=("$HOST_ALIGNMENT_DIR"/HOST_*.sam)
  featureCounts -T "$THREADS" -a "$HOST_GFF" -O -s 1 -g "$HOST_GENE_ATTRIBUTE" -t CDS \
    -o "$HOST_ALIGNMENT_DIR/featurecounts_HOST_summary.txt" "${HOST_sam_filenames[@]}"
  sed 's/\t/,/g' "$HOST_ALIGNMENT_DIR/featurecounts_HOST_summary.txt" \
    > "$HOST_ALIGNMENT_DIR/featurecounts_HOST_summary.csv"
  featureCounts -T "$THREADS" -a "$HOST_GFF" -s 1 -g "$HOST_GENE_ATTRIBUTE" -t rRNA \
    -o "$HOST_ALIGNMENT_DIR/featurecounts_HOST_rRNA.txt" "${HOST_sam_filenames[@]}"
fi

# Calculate coverage with respect to the bacterial reference with samtools.
for i in "${BACTERIA_sam_filenames[@]}"
do
    i_basename=$(basename "$i" .sam)
    echo "analyzing $i_basename ..."
    samtools view -@ "$THREADS" -bS "$i" > "$BACTERIA_ALIGNMENT_DIR/${i_basename}.bam"
    samtools sort -@ "$THREADS" -o "$BACTERIA_ALIGNMENT_DIR/${i_basename}.sorted.bam" "$BACTERIA_ALIGNMENT_DIR/${i_basename}.bam"
    samtools index -@ "$THREADS" "$BACTERIA_ALIGNMENT_DIR/${i_basename}.sorted.bam"
    samtools coverage -mA "$BACTERIA_ALIGNMENT_DIR/${i_basename}.sorted.bam"
    samtools coverage -m -o "$BACTERIA_ALIGNMENT_DIR/${i_basename}_coverage.txt" "$BACTERIA_ALIGNMENT_DIR/${i_basename}.sorted.bam"
    profile_args=(
      --bam "$BACTERIA_ALIGNMENT_DIR/${i_basename}.sorted.bam"
      --gff "$BACTERIA_GFF"
      --gene-id "$BACTERIA_GENE_ATTRIBUTE"
      --bins "$GENE_POSITION_BINS"
      --min-feature-reads "$METAGENE_MIN_FEATURE_READS"
      --output-prefix "$BACTERIA_ALIGNMENT_DIR/${i_basename}"
    )
    python3 "$SCRIPT_DIR/gene_position_profile.py" "${profile_args[@]}"
done

# Run python script to convert results to csv file.
(
  cd "$OUTPUT_DIR"
  python3 "$CSV_CONVERSION_SCRIPT"
  python3 "$SCRIPT_DIR/read_breakdown_csv.py"
)
