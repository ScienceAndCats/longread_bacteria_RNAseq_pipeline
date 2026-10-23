#!/usr/bin/env bash
set -euo pipefail

CONFIG_FILE="${1:-config.env}"
if (( $# )); then shift; fi
PREPARE_ONLY=false
SKIP_SUMMARY=false
SHARED_PREPARED_DIR=""
while (( $# )); do
  case "$1" in
    --prepare-only) PREPARE_ONLY=true; shift ;;
    --skip-summary) SKIP_SUMMARY=true; shift ;;
    --prepared-dir)
      [[ $# -ge 2 && -n "$2" ]] || { echo "--prepared-dir requires a directory" >&2; exit 1; }
      SHARED_PREPARED_DIR="$2"; shift 2 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done
if [[ "$PREPARE_ONLY" == true && -n "$SHARED_PREPARED_DIR" ]]; then
  echo "--prepare-only and --prepared-dir cannot be combined" >&2; exit 1
fi
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
# Legacy basename settings remain usable, but new configurations use folders.
DECOY_REFERENCE_DIR="${DECOY_REFERENCE_DIR-${DECOY_BBSPLIT_REFERENCE:-}}"
BACTERIA_REFERENCE_DIR="${BACTERIA_REFERENCE_DIR-${BACTERIA_BBSPLIT_REFERENCE:-}}"
HOST_REFERENCE_DIR="${HOST_REFERENCE_DIR-${HOST_BBSPLIT_REFERENCE:-}}"
: "${DECOY_REFERENCE_DIR:?Set DECOY_REFERENCE_DIR in $CONFIG_FILE}"
: "${BACTERIA_REFERENCE_DIR:?Set BACTERIA_REFERENCE_DIR in $CONFIG_FILE}"
HOST_ENABLED=false
[[ -z "$HOST_REFERENCE_DIR" ]] || HOST_ENABLED=true
: "${BBSPLIT_K:=11}"
: "${MIN_READ_LENGTH:=20}"
: "${THREADS:=16}"
: "${GENE_POSITION_BINS:=100}"
: "${METAGENE_MIN_FEATURE_READS:=10}"
: "${CSV_CONVERSION_SCRIPT:=bacteria_with_decoys_csvConversion.py}"

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be a positive integer" >&2; exit 1; }
[[ "$GENE_POSITION_BINS" =~ ^[1-9][0-9]*$ ]] || { echo "GENE_POSITION_BINS must be a positive integer" >&2; exit 1; }
[[ "$METAGENE_MIN_FEATURE_READS" =~ ^[1-9][0-9]*$ ]] || { echo "METAGENE_MIN_FEATURE_READS must be a positive integer" >&2; exit 1; }
[[ "$BBSPLIT_K" =~ ^([1-9]|1[0-5])$ ]] || { echo "BBSPLIT_K must be an integer from 1 to 15" >&2; exit 1; }
[[ "$MIN_READ_LENGTH" =~ ^[1-9][0-9]*$ ]] && (( MIN_READ_LENGTH <= 6000 )) || {
  echo "MIN_READ_LENGTH must be a positive integer no greater than 6000" >&2; exit 1;
}
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fastq_files=(); input_sample_names=()
if [[ "$PREPARE_ONLY" != true ]]; then
  # Use the same file selection and collision checks as scheduled runs before
  # preparing references or writing sample outputs. Names come from the selected
  # filename, including when it is a symlink to a differently named FASTQ.
  mkdir -p "$OUTPUT_DIR/.pipeline"
  python3 - "$SCRIPT_DIR" "$FASTQ_DIR" "$FASTQ_GLOB" \
    "$OUTPUT_DIR/.pipeline/local_samples.json" > "$OUTPUT_DIR/.pipeline/local_samples.inputs" <<'SAMPLES'
import json, os, sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
from sample_manifest import discover_samples
try:
    samples = discover_samples(sys.argv[2], sys.argv[3])
except (ValueError, OSError) as exc:
    raise SystemExit(f"ERROR: {exc}")
Path(sys.argv[4]).write_text(json.dumps(samples, indent=2) + "\n")
for sample in samples:
    for value in (sample["reads"][0]["path"], sample["sample"]):
        sys.stdout.buffer.write(os.fsencode(value) + b"\0")
SAMPLES
  while IFS= read -r -d '' input_fastq && IFS= read -r -d '' input_sample_name; do
    fastq_files+=("$input_fastq")
    input_sample_names+=("$input_sample_name")
  done < "$OUTPUT_DIR/.pipeline/local_samples.inputs"
fi

ALIGNMENT_DIR="$OUTPUT_DIR/bbsplit_alignments"
PREPARED_DIR="${SHARED_PREPARED_DIR:-$ALIGNMENT_DIR}"
if [[ -n "$SHARED_PREPARED_DIR" && ! -f "$PREPARED_DIR/.prepared" ]]; then
  echo "Shared references are incomplete: $PREPARED_DIR" >&2; exit 1
fi
DECOY_ALIGNMENT_DIR="$ALIGNMENT_DIR/decoy"
BACTERIA_ALIGNMENT_DIR="$ALIGNMENT_DIR/bacteria"
HOST_ALIGNMENT_DIR="$ALIGNMENT_DIR/host"
LEFTOVER_READS_DIR="$ALIGNMENT_DIR/leftover_reads"
OVERSIZED_READS_DIR="$ALIGNMENT_DIR/oversized_reads"
mkdir -p "$OUTPUT_DIR" "$DECOY_ALIGNMENT_DIR" "$BACTERIA_ALIGNMENT_DIR" "$LEFTOVER_READS_DIR" "$OVERSIZED_READS_DIR"
if [[ "$HOST_ENABLED" == true ]]; then
  mkdir -p "$HOST_ALIGNMENT_DIR"
fi
if [[ "$CSV_CONVERSION_SCRIPT" != /* ]]; then
  CSV_CONVERSION_SCRIPT="$SCRIPT_DIR/$CSV_CONVERSION_SCRIPT"
fi

# Normalize all reference shards and annotation coordinate systems before mapping.
reference_args=(--decoy-dir "$DECOY_REFERENCE_DIR" --bacteria-dir "$BACTERIA_REFERENCE_DIR"
  --output-dir "$PREPARED_DIR/references")
if [[ "$HOST_ENABLED" == true ]]; then
  reference_args+=(--host-dir "$HOST_REFERENCE_DIR")
fi
if [[ -z "$SHARED_PREPARED_DIR" ]]; then
  rm -f -- "$PREPARED_DIR/.prepared"
  python3 "$SCRIPT_DIR/prepare_references.py" "${reference_args[@]}"
fi
# Read only fixed manifest fields; generated paths are never evaluated as shell code.
mapfile -t reference_paths < <(python3 - "$PREPARED_DIR/references/manifest.json" <<'MANIFEST'
import json, sys
from pathlib import Path
manifest = json.loads(Path(sys.argv[1]).read_text())
print(manifest["decoy"]["fasta"])
for organism in ("bacteria", "host"):
    entry = manifest.get(organism) or {}
    for molecule in ("dna", "rna"):
        reference = entry.get(molecule) or {}
        for key in ("fasta", "gff"):
            print(reference.get(key, ""))
MANIFEST
)
DECOY_FASTA="${reference_paths[0]}"
BACTERIA_DNA_FASTA="${reference_paths[1]}"; BACTERIA_DNA_GFF="${reference_paths[2]}"
BACTERIA_RNA_FASTA="${reference_paths[3]}"; BACTERIA_RNA_GFF="${reference_paths[4]}"
HOST_DNA_FASTA="${reference_paths[5]}"; HOST_DNA_GFF="${reference_paths[6]}"
HOST_RNA_FASTA="${reference_paths[7]}"; HOST_RNA_GFF="${reference_paths[8]}"
BACTERIA_GFF="$PREPARED_DIR/references/bacteria.combined.gff3"
HOST_GFF=""
BACTERIA_GENE_ATTRIBUTE=gene_id
HOST_GENE_ATTRIBUTE=gene_id
if [[ "$HOST_ENABLED" == true ]]; then
  HOST_GFF="$PREPARED_DIR/references/host.combined.gff3"
fi
if [[ -z "$SHARED_PREPARED_DIR" ]]; then
  cat -- "$BACTERIA_DNA_GFF" > "$BACTERIA_GFF"
  if [[ -n "$BACTERIA_RNA_GFF" ]]; then cat -- "$BACTERIA_RNA_GFF" >> "$BACTERIA_GFF"; fi
  if [[ "$HOST_ENABLED" == true ]]; then
    cat -- "$HOST_DNA_GFF" > "$HOST_GFF"
    if [[ -n "$HOST_RNA_GFF" ]]; then cat -- "$HOST_RNA_GFF" >> "$HOST_GFF"; fi
  fi
fi

# Each stage, k value, and reference content uses an isolated index. BBSplit's
# own reference-name cache would otherwise reuse stale sequence after an edit.
run_bbsplit() {
  # BBTools' launcher reconstructs its Java command with eval, so retain Bash
  # quoting for paths containing spaces or shell metacharacters.
  local argument escaped
  local quoted_arguments=()
  for argument in "$@"; do
    printf -v escaped '%q' "$argument"
    quoted_arguments+=("$escaped")
  done
  bbsplit.sh "${quoted_arguments[@]}"
}

prepare_bbsplit_index() {
  local fasta_file="$1" stage="$2"
  local fingerprint index_dir
  if [[ -n "$SHARED_PREPARED_DIR" ]]; then
    # The preparation job already checksummed the large reference files.
    # Workers load its frozen index paths instead of re-reading human genomes.
    index_dir=$(python3 - "$PREPARED_DIR/indexes.json" "$stage" "$BBSPLIT_K" <<'INDEX'
import json, sys
from pathlib import Path
indexes = json.loads(Path(sys.argv[1]).read_text())
if indexes["k"] != sys.argv[3]:
    raise SystemExit("Shared indexes use a different BBSPLIT_K")
print(indexes["stages"][sys.argv[2]])
INDEX
) || return 1
    [[ -f "$index_dir/.ready" ]] || { echo "Shared index is incomplete: $index_dir" >&2; return 1; }
    printf '%s\n' "$index_dir"
    return 0
  fi
  fingerprint=$(sha256sum "$fasta_file") || return 1
  fingerprint=${fingerprint%% *}
  index_dir="$PREPARED_DIR/indexes/$stage/k$BBSPLIT_K/$fingerprint"
  mkdir -p "$index_dir" || return 1
  echo "preparing $stage BBSplit index from $fasta_file ..." >&2
  # Explicit set names keep FASTA filenames out of SAM sequence identifiers.
  if ! run_bbsplit "ref_${stage}=$fasta_file" "path=$index_dir" build=1 "k=$BBSPLIT_K" \
    "threads=$THREADS" mapmode=pb msa=MultiStateAligner9PacBio fastareadlen=6000 \
    qtrim=f interleaved=f secondary=f ambiguous=best trimrefdescription=t \
    > "$index_dir/index.stdout.txt" 2> "$index_dir/index.stderr.txt"; then
    cat "$index_dir/index.stderr.txt" >&2
    return 1
  fi
  touch "$index_dir/.ready"
  printf '%s\n' "$index_dir"
}

DECOY_BBSPLIT_INDEX=$(prepare_bbsplit_index "$DECOY_FASTA" "decoy")
BACTERIA_DNA_INDEX=$(prepare_bbsplit_index "$BACTERIA_DNA_FASTA" "bacteria")
BACTERIA_RNA_INDEX=""; HOST_RNA_INDEX=""
if [[ -n "$BACTERIA_RNA_FASTA" ]]; then
  BACTERIA_RNA_INDEX=$(prepare_bbsplit_index "$BACTERIA_RNA_FASTA" "bacteria_rna")
fi
if [[ "$HOST_ENABLED" == true ]]; then
  HOST_DNA_INDEX=$(prepare_bbsplit_index "$HOST_DNA_FASTA" "host")
  if [[ -n "$HOST_RNA_FASTA" ]]; then
    HOST_RNA_INDEX=$(prepare_bbsplit_index "$HOST_RNA_FASTA" "host_rna")
  fi
fi

if [[ -z "$SHARED_PREPARED_DIR" ]]; then
  python3 - "$PREPARED_DIR/indexes.json" "$BBSPLIT_K" \
    "$DECOY_BBSPLIT_INDEX" "$BACTERIA_DNA_INDEX" "$BACTERIA_RNA_INDEX" \
    "${HOST_DNA_INDEX:-}" "$HOST_RNA_INDEX" <<'INDEXES'
import json, sys
from pathlib import Path
names = ("decoy", "bacteria", "bacteria_rna", "host", "host_rna")
Path(sys.argv[1]).write_text(json.dumps({"k": sys.argv[2], "stages": dict(zip(names, sys.argv[3:]))}, indent=2) + "\n")
INDEXES
  touch "$PREPARED_DIR/.prepared"
fi
if [[ "$PREPARE_ONLY" == true ]]; then
  echo "Shared reference preparation complete: $PREPARED_DIR"
  exit 0
fi

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
opener = gzip.open if input_path.lower().endswith(".gz") else open
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
    for idx in "${!fastq_files[@]}"; do
      input_fastq="${fastq_files[$idx]}"
      sampled_fastq="$sample_dir/${input_sample_names[$idx]}.fastq"
      echo "sampling up to $SAMPLE_SIZE reads from $(basename "$input_fastq") ..."
      sample_fastq "$input_fastq" "$sampled_fastq"
      sampled_fastqs+=("$sampled_fastq")
    done
    fastq_files=("${sampled_fastqs[@]}")
    ;;
  false) ;;
  *) echo "SAMPLE_READS must be true or false" >&2; exit 1 ;;
esac

# Complete all adapter trimming and the minimum length filter before routing
# reads above 6000 bp out of the alignment/analysis pipeline.
trimmed_fastqs=(); sample_names=()
for idx in "${!fastq_files[@]}"; do
  input_fastq="${fastq_files[$idx]}"; sample="${input_sample_names[$idx]}"
  trimmed_fastq="$OUTPUT_DIR/${sample}_${MIN_READ_LENGTH}bp.trim.fastq"
  all_trimmed_fastq="$OUTPUT_DIR/.${sample}_${MIN_READ_LENGTH}bp.adapter_trimmed.fastq"
  bash "$SCRIPT_DIR/trim_fastq.sh" "$CONFIG_FILE" "$input_fastq" "$all_trimmed_fastq" \
    "$OUTPUT_DIR/${sample}_${MIN_READ_LENGTH}bp"
  python3 "$SCRIPT_DIR/split_trimmed_fastq.py" "$all_trimmed_fastq" "$trimmed_fastq" \
    "$OVERSIZED_READS_DIR/${sample}_${MIN_READ_LENGTH}bp_over_6000bp.fastq.gz" \
    "$OUTPUT_DIR/${sample}_${MIN_READ_LENGTH}bp.length_filter.tsv"
  rm -- "$all_trimmed_fastq"
  trimmed_fastqs+=("$trimmed_fastq"); sample_names+=("$sample")
done

# Count complete FASTQ records (plain or gzip-compressed).
fastq_read_count() {
  if [[ "$1" == *.gz ]]; then gzip -cd -- "$1"; else cat -- "$1"; fi |
    awk 'END { if (NR % 4) exit 1; print NR / 4 }'
}

# Align reads with BBSplit's long-read mapper, retain only primary
# alignments in SAM, and write unmapped reads for the next classification stage.
# The small TSV report is intentionally stable input for both CSV generators.
bbsplit_stage() {
  local index_dir="$1" sam="$2" report="$3" unmapped="$4" input="$5"
  local input_count unmapped_count aligned_count primary_count raw_log="${report%.txt}.stderr.txt"
  input_count=$(fastq_read_count "$input")
  run_bbsplit "path=$index_dir" build=1 "in=$input" out=stdout.sam overwrite=t \
    "threads=$THREADS" "k=$BBSPLIT_K" mapmode=pb \
    msa=MultiStateAligner9PacBio fastareadlen=6000 \
    qtrim=f interleaved=f secondary=f ambiguous=best trimrefdescription=t \
    2> "$raw_log" | samtools view -h -F 2304 -o "$sam" -
  samtools fastq -@ "$THREADS" -f 4 -0 "$unmapped" -s /dev/null -n "$sam" 2>> "$raw_log"
  unmapped_count=$(fastq_read_count "$unmapped")
  primary_count=$(samtools view -c "$sam")
  aligned_count=$(samtools view -c -F 4 "$sam")
  if (( primary_count != input_count || aligned_count + unmapped_count != input_count )); then
    echo "ERROR: BBSplit read counts do not balance for $input; see $raw_log" >&2
    return 1
  fi
  printf 'metric\tcount\ninput\t%s\naligned\t%s\nunmapped\t%s\n' \
    "$input_count" "$aligned_count" "$unmapped_count" > "$report"
}

# Map to bacterial decoys and pass only primary-unmapped reads onward.
decoy_unmapped_fastqs=()
for idx in "${!trimmed_fastqs[@]}"; do
  i="${trimmed_fastqs[$idx]}"; i_basename="${sample_names[$idx]}_${MIN_READ_LENGTH}bp"
  echo "mapping to bacterial decoys with BBSplit: $i_basename ..."
  unmapped="$DECOY_ALIGNMENT_DIR/${i_basename}_unmapped_to_other_bugs.fastq.gz"
  bbsplit_stage "$DECOY_BBSPLIT_INDEX" \
    "$DECOY_ALIGNMENT_DIR/${i_basename}.mapped_to_other_bugs.sam" \
    "$DECOY_ALIGNMENT_DIR/${i_basename}.mapped_to_other_bugs.bbsplit.txt" \
    "$unmapped" "$i"
  decoy_unmapped_fastqs+=("$unmapped")
done

# Keep all RNA mapping ahead of DNA, including RNA from the other organism.
# Every stage receives only primary-unmapped reads from the previous stage.
BACTERIA_sam_filenames=(); HOST_sam_filenames=()
for idx in "${!decoy_unmapped_fastqs[@]}"; do
  current_input="${decoy_unmapped_fastqs[$idx]}"
  i_basename="${sample_names[$idx]}_${MIN_READ_LENGTH}bp"
  for stage in bacteria_rna host_rna bacteria_dna host_dna; do
    case "$stage" in
      bacteria_rna)
        [[ -n "$BACTERIA_RNA_INDEX" ]] || continue
        index="$BACTERIA_RNA_INDEX"; organism=bacteria; molecule=rna; prefix=BACTERIA ;;
      host_rna)
        [[ -n "$HOST_RNA_INDEX" ]] || continue
        index="$HOST_RNA_INDEX"; organism=host; molecule=rna; prefix=HOST ;;
      bacteria_dna)
        index="$BACTERIA_DNA_INDEX"; organism=bacteria; molecule=dna; prefix=BACTERIA ;;
      host_dna)
        [[ "$HOST_ENABLED" == true ]] || continue
        index="$HOST_DNA_INDEX"; organism=host; molecule=dna; prefix=HOST ;;
    esac
    stage_dir="$ALIGNMENT_DIR/$organism/$molecule"
    mkdir -p "$stage_dir"
    unmapped="$stage_dir/${i_basename}_unmapped.fastq.gz"
    echo "mapping to $organism $molecule references with BBSplit: $i_basename ..."
    bbsplit_stage "$index" "$stage_dir/${prefix}_${i_basename}.sam" \
      "$stage_dir/${prefix}_${i_basename}.bbsplit.txt" "$unmapped" "$current_input"
    current_input="$unmapped"
  done
  cp -- "$current_input" "$LEFTOVER_READS_DIR/${i_basename}_leftover.fastq.gz"
  finalize_args=(--alignment-dir "$ALIGNMENT_DIR" --sample-prefix "$i_basename")
  if [[ "$HOST_ENABLED" == true ]]; then finalize_args+=(--host-enabled); fi
  if [[ -n "$BACTERIA_RNA_INDEX" ]]; then finalize_args+=(--bacteria-rna); fi
  if [[ -n "$HOST_RNA_INDEX" ]]; then finalize_args+=(--host-rna); fi
  python3 "$SCRIPT_DIR/finalize_sample_alignments.py" "${finalize_args[@]}"
  # This compatibility FASTQ contains reads assigned to neither decoy nor bacteria.
  # Include host RNA reads removed before the bacterial DNA stage.
  bacterial_unmapped="$BACTERIA_ALIGNMENT_DIR/${i_basename}_unmapped_to_bacteria.fastq.gz"
  if [[ -n "$HOST_RNA_INDEX" ]]; then
    samtools fastq -@ "$THREADS" -F 2308 -0 "$BACTERIA_ALIGNMENT_DIR/.${i_basename}_host_rna.fastq.gz" \
      -s /dev/null -n "$HOST_ALIGNMENT_DIR/rna/HOST_${i_basename}.sam"
    cat -- "$BACTERIA_ALIGNMENT_DIR/.${i_basename}_host_rna.fastq.gz" \
      "$BACTERIA_ALIGNMENT_DIR/dna/${i_basename}_unmapped.fastq.gz" > "$bacterial_unmapped"
    rm -- "$BACTERIA_ALIGNMENT_DIR/.${i_basename}_host_rna.fastq.gz"
  else
    cp -- "$BACTERIA_ALIGNMENT_DIR/dna/${i_basename}_unmapped.fastq.gz" "$bacterial_unmapped"
  fi
  for organism in bacteria host; do
    if [[ "$organism" == host && "$HOST_ENABLED" != true ]]; then continue; fi
    if [[ "$organism" == bacteria ]]; then
      prefix=BACTERIA; dna_gff="$BACTERIA_DNA_GFF"; rna_gff="$BACTERIA_RNA_GFF"
      BACTERIA_sam_filenames+=("$BACTERIA_ALIGNMENT_DIR/${prefix}_${i_basename}.sam")
    else
      prefix=HOST; dna_gff="$HOST_DNA_GFF"; rna_gff="$HOST_RNA_GFF"
      HOST_sam_filenames+=("$HOST_ALIGNMENT_DIR/${prefix}_${i_basename}.sam")
    fi
    organism_dir="$ALIGNMENT_DIR/$organism"
    category_args=(--dna-sam "$organism_dir/dna/${prefix}_${i_basename}.sam" --dna-gff "$dna_gff"
      --output "$organism_dir/${prefix}_${i_basename}.categories.tsv")
    if [[ -n "$rna_gff" ]]; then
      category_args+=(--rna-sam "$organism_dir/rna/${prefix}_${i_basename}.sam" --rna-gff "$rna_gff")
    fi
    python3 "$SCRIPT_DIR/classify_alignments.py" "${category_args[@]}"
  done
done

# Create featureCounts summary file.
# Assign mapped reads to bacterial genes, retaining the pipeline's stranded and
# overlapping-feature behavior.
python3 "$SCRIPT_DIR/count_features.py" --threads "$THREADS" --annotation "$BACTERIA_GFF" \
  --allow-overlap --gene-id "$BACTERIA_GENE_ATTRIBUTE" --feature-type CDS \
  --output "$OUTPUT_DIR/featurecounts_BACTERIA_summary.txt" "${BACTERIA_sam_filenames[@]}"
sed 's/\t/,/g' "$OUTPUT_DIR/featurecounts_BACTERIA_summary.txt" > "$OUTPUT_DIR/featurecounts_BACTERIA_summary.csv"
# Retain the independent rRNA featureCounts output for per-feature analysis.
# Exclusive read categories for the breakdown come from classify_alignments.py.
python3 "$SCRIPT_DIR/count_features.py" --threads "$THREADS" --annotation "$BACTERIA_GFF" \
  --gene-id "$BACTERIA_GENE_ATTRIBUTE" --feature-type rRNA \
  --output "$OUTPUT_DIR/featurecounts_BACTERIA_rRNA.txt" "${BACTERIA_sam_filenames[@]}"

# When a host reference is configured, independently assign host alignments to
# its CDS features and keep the results alongside the host BBSplit outputs.
if [[ -n "$HOST_GFF" ]]; then
  python3 "$SCRIPT_DIR/count_features.py" --threads "$THREADS" --annotation "$HOST_GFF" \
    --allow-overlap --gene-id "$HOST_GENE_ATTRIBUTE" --feature-type CDS \
    --output "$HOST_ALIGNMENT_DIR/featurecounts_HOST_summary.txt" "${HOST_sam_filenames[@]}"
  sed 's/\t/,/g' "$HOST_ALIGNMENT_DIR/featurecounts_HOST_summary.txt" \
    > "$HOST_ALIGNMENT_DIR/featurecounts_HOST_summary.csv"
  python3 "$SCRIPT_DIR/count_features.py" --threads "$THREADS" --annotation "$HOST_GFF" \
    --gene-id "$HOST_GENE_ATTRIBUTE" --feature-type rRNA \
    --output "$HOST_ALIGNMENT_DIR/featurecounts_HOST_rRNA.txt" "${HOST_sam_filenames[@]}"
fi

# Calculate coverage with respect to the bacterial reference with samtools.
for i in "${BACTERIA_sam_filenames[@]}"
do
    i_basename=$(basename "$i" .sam)
    echo "analyzing $i_basename ..."
    samtools view -@ "$THREADS" -bS "$i" > "$BACTERIA_ALIGNMENT_DIR/${i_basename}.bam"
    samtools sort -@ "$THREADS" -o "$BACTERIA_ALIGNMENT_DIR/${i_basename}.sorted.bam" "$BACTERIA_ALIGNMENT_DIR/${i_basename}.bam"
    samtools index -@ "$THREADS" "$BACTERIA_ALIGNMENT_DIR/${i_basename}.sorted.bam"
    # Genome coverage describes DNA fallback alignments; RNA has transcript coordinates.
    dna_bam="$BACTERIA_ALIGNMENT_DIR/dna/${i_basename}.sorted.bam"
    samtools sort -@ "$THREADS" -o "$dna_bam" "$BACTERIA_ALIGNMENT_DIR/dna/${i_basename}.sam"
    samtools coverage -mA "$dna_bam"
    samtools coverage -m -o "$BACTERIA_ALIGNMENT_DIR/${i_basename}_coverage.txt" "$dna_bam"
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

# Run project summaries after all samples, or defer to the dependent Slurm job.
if [[ "$SKIP_SUMMARY" != true ]]; then
(
  cd "$OUTPUT_DIR"
  python3 "$CSV_CONVERSION_SCRIPT"
  python3 "$SCRIPT_DIR/read_breakdown_csv.py" --host-enabled "$HOST_ENABLED"
)
fi
