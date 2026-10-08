#!/usr/bin/env bash
# Trim a sample independently of reference selection/mapping.
set -euo pipefail
[[ $# == 4 ]] || { echo "Usage: bash $0 config.env input.fastq[.gz] output.fastq[.gz] report_prefix" >&2; exit 1; }
source "$1"
input_fastq=$2; trimmed_fastq=$3; report_prefix=$4
: "${ADAPTER_TRIMMING_MODE:=legacy}"
: "${CUTADAPT_ERROR_RATE:=0.15}"
: "${CUTADAPT_MIN_OVERLAP:=12}"
: "${ADAPTER_MAX_TRIMMING_ROUNDS:=10}"
: "${MIN_READ_LENGTH:=20}"
: "${THREADS:=16}"
: "${TRIM_NANOPORE_ADAPTERS:=true}"
: "${TRIM_TRUSEQ_SMALL_RNA:=false}"
: "${TRIM_ILLUMINA_ADAPTERS:=false}"
: "${ADAPTER_NANOPORE:=TTTCTGTTGGTGCTGATATTGC}"
: "${ADAPTER_NEB_5P:=GTTCAGAGTTCTACAGTCCGACGATC}"
: "${ADAPTER_NEB_3P:=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC}"
: "${ADAPTER_TRUSEQ_SMALL_RNA:=TGGAATTCTCGGGTGCCAAGG}"
: "${ADAPTER_ILLUMINA:=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA}"
: "${ADAPTER_ILLUMINA_SMALL_RNA:=TGGAATTCTCGGGTGCCAAGG}"
for boolean in TRIM_NANOPORE_ADAPTERS TRIM_TRUSEQ_SMALL_RNA TRIM_ILLUMINA_ADAPTERS; do
  case "${!boolean,,}" in
    true|false) ;;
    *) echo "$boolean must be true or false" >&2; exit 1 ;;
  esac
done
[[ "$THREADS" =~ ^[1-9][0-9]*$ && "$CUTADAPT_MIN_OVERLAP" =~ ^[1-9][0-9]*$ && "$ADAPTER_MAX_TRIMMING_ROUNDS" =~ ^[1-9][0-9]*$ && "$MIN_READ_LENGTH" =~ ^[0-9]+$ ]] || {
  echo "Invalid thread, overlap, round limit or minimum length" >&2; exit 1;
}
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "$ADAPTER_TRIMMING_MODE" == neb_e7330 ]]; then
  neb_args=(--input "$input_fastq" --output "$trimmed_fastq" --report-prefix "$report_prefix"
    --five "$ADAPTER_NEB_5P" --three "$ADAPTER_NEB_3P" --nanopore "$ADAPTER_NANOPORE"
    --error-rate "$CUTADAPT_ERROR_RATE" --overlap "$CUTADAPT_MIN_OVERLAP"
    --max-rounds "$ADAPTER_MAX_TRIMMING_ROUNDS" --minimum-length "$MIN_READ_LENGTH")
  if [[ -n "${ADAPTER_NANOPORE_RC:-}" ]]; then neb_args+=(--nanopore-rc "$ADAPTER_NANOPORE_RC"); fi
  if [[ "${TRIM_NANOPORE_ADAPTERS,,}" == false ]]; then neb_args+=(--no-nanopore); fi
  if [[ "${TRIM_TRUSEQ_SMALL_RNA,,}" == true ]]; then neb_args+=(--truseq "$ADAPTER_TRUSEQ_SMALL_RNA"); fi
  python3 "$script_dir/trim_neb_adapters.py" "${neb_args[@]}" > "${report_prefix}.cutadapt_log.txt"
elif [[ "$ADAPTER_TRIMMING_MODE" == legacy ]]; then
  # Keep the original adapter choices, but filter length only after every pass.
  temp_dir=$(mktemp -d "$(dirname "$trimmed_fastq")/.adapter_trim.XXXXXX")
  trap 'rm -rf -- "$temp_dir"' EXIT
  common=(-j "$THREADS" -e "$CUTADAPT_ERROR_RATE" -O "$CUTADAPT_MIN_OVERLAP")
  nanopore_args=()
  if [[ "${TRIM_NANOPORE_ADAPTERS,,}" == true ]]; then nanopore_args=(-a "$ADAPTER_NANOPORE"); fi
  cutadapt "${common[@]}" "${nanopore_args[@]}" -o "$temp_dir/trim.fastq" "$input_fastq" > "$temp_dir/nanopore.log"
  final_input="$temp_dir/trim.fastq"
  if [[ "${TRIM_ILLUMINA_ADAPTERS,,}" == true ]]; then
    cutadapt "${common[@]}" --times 2 -a "$ADAPTER_ILLUMINA" -a "$ADAPTER_ILLUMINA_SMALL_RNA" \
      -o "$temp_dir/illumina.fastq" "$final_input" > "${report_prefix}.illumina_adapter_log.txt"
    final_input="$temp_dir/illumina.fastq"
  fi
  cutadapt "${common[@]}" -m "$MIN_READ_LENGTH" --json "${report_prefix}.cutadapt.json" \
    -o "$trimmed_fastq" "$final_input" > "${report_prefix}.cutadapt_log.txt"
  cat "$temp_dir/nanopore.log" > "${report_prefix}.nanopore_adapter_log.txt"
else
  echo "ADAPTER_TRIMMING_MODE must be legacy or neb_e7330" >&2; exit 1
fi
