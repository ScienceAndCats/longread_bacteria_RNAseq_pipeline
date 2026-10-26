#!/usr/bin/env bash
# Extract a linked insert once, then trim remaining terminal adapter layers.
set -euo pipefail
[[ $# == 4 ]] || { echo "Usage: bash $0 config.env input.fastq[.gz] output.fastq[.gz] report_prefix" >&2; exit 1; }
# shellcheck source=/dev/null
source "$1"
input_fastq=$2; trimmed_fastq=$3; report_prefix=$4
if [[ "$input_fastq" -ef "$trimmed_fastq" ]]; then
  echo "Input and output FASTQ must be different files" >&2
  exit 1
fi
: "${ADAPTER_NEB_5P=GTTCAGAGTTCTACAGTCCGACGATC}"
: "${ADAPTER_NEB_3P=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC}"
: "${ADAPTER_SEQUENCES=TTTCTGTTGGTGCTGATATTGC GTTCAGAGTTCTACAGTCCGACGATC AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT AATGATACGGCGACCACCGAGATCTACAC CAAGCAGAAGACGGCATACGAGAT}"
: "${ADAPTER_REVERSE_COMPLEMENTS=true}"
: "${CUTADAPT_ERROR_RATE=0.15}"
: "${CUTADAPT_MIN_OVERLAP=12}"
: "${ADAPTER_MAX_TRIMMING_ROUNDS=3}"
: "${MIN_READ_LENGTH=20}"
: "${THREADS=16}"

for positive in THREADS CUTADAPT_MIN_OVERLAP ADAPTER_MAX_TRIMMING_ROUNDS; do
  [[ "${!positive}" =~ ^[1-9][0-9]*$ ]] || {
    echo "$positive must be a positive integer" >&2; exit 1;
  }
done
[[ "$MIN_READ_LENGTH" =~ ^[1-9][0-9]{0,3}$ ]] && (( MIN_READ_LENGTH <= 6000 )) || {
  echo "MIN_READ_LENGTH must be a positive integer no greater than 6000" >&2; exit 1;
}
case "${ADAPTER_REVERSE_COMPLEMENTS,,}" in
  true|false) ;;
  *) echo "ADAPTER_REVERSE_COMPLEMENTS must be true or false" >&2; exit 1 ;;
esac

# Validate both linked flanks and the independent terminal adapter list.
validated=$(python3 - "$CUTADAPT_ERROR_RATE" "$ADAPTER_NEB_5P" "$ADAPTER_NEB_3P" \
  "$ADAPTER_SEQUENCES" "$ADAPTER_REVERSE_COMPLEMENTS" <<'PY'
import math
import sys

rate, five, three, text, reverse = sys.argv[1:]
try:
    rate = float(rate)
except ValueError:
    raise SystemExit("CUTADAPT_ERROR_RATE must be a finite number with 0 <= value < 1")
if not math.isfinite(rate) or not 0 <= rate < 1:
    raise SystemExit("CUTADAPT_ERROR_RATE must be a finite number with 0 <= value < 1")
print(rate)
alphabet = set("ACGTRYSWKMBDHVN")
for name, sequence in (("ADAPTER_NEB_5P", five), ("ADAPTER_NEB_3P", three)):
    sequence = sequence.strip().upper()
    if not sequence or set(sequence) - alphabet:
        raise SystemExit(f"{name} must be a nonempty DNA IUPAC sequence")
    print(sequence)
complement = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
adapters = []
for sequence in text.upper().split():
    if set(sequence) - alphabet:
        raise SystemExit("ADAPTER_SEQUENCES must contain whitespace-separated DNA IUPAC sequences")
    for adapter in (sequence, sequence.translate(complement)[::-1]) if reverse.lower() == "true" else (sequence,):
        if adapter not in adapters:
            adapters.append(adapter)
for adapter in adapters:
    print(adapter)
PY
)
mapfile -t validated_lines <<< "$validated"
error_rate=${validated_lines[0]}
linked_five=${validated_lines[1]}
linked_three=${validated_lines[2]}
adapters=("${validated_lines[@]:3}")
five_args=(); three_args=()
for adapter in "${adapters[@]}"; do
  # X is Cutadapt's non-internal adapter sentinel: full and partial terminal
  # matches are allowed, while matches wholly inside an insert are excluded.
  five_args+=(-g "X$adapter")
  three_args+=(-a "${adapter}X")
done

mkdir -p -- "$(dirname "$trimmed_fastq")" "$(dirname "$report_prefix")"
for old_log in "${report_prefix}".adapter_round*.5p.log "${report_prefix}".adapter_round*.3p.log; do
  [[ -f "$old_log" ]] || continue
  round_suffix=${old_log#"$report_prefix".adapter_round}
  if [[ "$round_suffix" =~ ^[1-9][0-9]*\.(5p|3p)\.log$ ]]; then rm -- "$old_log"; fi
done
temp_dir=$(mktemp -d "$(dirname "$trimmed_fastq")/.adapter_trim.XXXXXX")
trap 'rm -rf -- "$temp_dir"' EXIT
current_input=$input_fastq
common=(-j "$THREADS" -e "$error_rate" -O "$CUTADAPT_MIN_OVERLAP" --times 1)
# This is a single extraction pass, independent of the terminal round limit.
# Require both unanchored flanks so outer library sequence can be removed.
# The rightmost 3' match retains intervening adapter copies for terminal cleanup.
# Normalize reverse-oriented linked inserts, preserving their original headers.
cutadapt "${common[@]}" --revcomp --rename '{header}' \
  -g "neb_linked=${linked_five};required...${linked_three};required;rightmost" \
  --json "${report_prefix}.linked_adapter.json" -o "$temp_dir/linked.fastq" "$current_input" \
  > "${report_prefix}.linked_adapter.log"
current_input="$temp_dir/linked.fastq"
if (( ${#adapters[@]} )); then
  for (( round=1; round<=ADAPTER_MAX_TRIMMING_ROUNDS; round++ )); do
    # Run passes sequentially so THREADS bounds concurrent Cutadapt workers.
    # Keep only the current and next temporary FASTQ while each pass runs.
    cutadapt "${common[@]}" "${five_args[@]}" -o "$temp_dir/five.fastq" "$current_input" \
      > "${report_prefix}.adapter_round${round}.5p.log"
    if [[ "$current_input" == "$temp_dir/three.fastq" || "$current_input" == "$temp_dir/linked.fastq" ]]; then
      rm -- "$current_input"
    fi
    cutadapt "${common[@]}" "${three_args[@]}" -o "$temp_dir/three.fastq" "$temp_dir/five.fastq" \
      > "${report_prefix}.adapter_round${round}.3p.log"
    rm -- "$temp_dir/five.fastq"
    current_input="$temp_dir/three.fastq"
  done
fi

# Neither extraction nor terminal trimming filters reads: the final native
# Cutadapt report retains the original count and applies the minimum last.
final_output="$temp_dir/final.fastq"
if [[ "$trimmed_fastq" == *.gz ]]; then final_output+=".gz"; fi
cutadapt -j "$THREADS" -m "$MIN_READ_LENGTH" --json "${report_prefix}.cutadapt.json" \
  -o "$final_output" "$current_input" > "${report_prefix}.cutadapt_log.txt"
mv -- "$final_output" "$trimmed_fastq"
