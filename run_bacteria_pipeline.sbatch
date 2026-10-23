#!/usr/bin/env bash
#SBATCH --job-name=bacteria_pipeline
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=bacteria_pipeline-%A_%a.out

set -euo pipefail

# Submit with: bash submit_bacteria_pipeline.sh config.env
# The helper supplies the array range (up to 40 samples, 8 concurrent by
# default) and the preparation/aggregation dependencies automatically.
ROLE="${1:?Use submit_bacteria_pipeline.sh to submit this template}"
MANIFEST="${2:?Missing run manifest}"
PIPELINE_DIR="${3:?Missing pipeline directory}"

# Account/email are passed to sbatch by the helper from the submitted config.
# Read that run's snapshot before activation, using only Bash to bootstrap.
RUN_CONFIG="$(dirname "$MANIFEST")/config.env"
[[ -f "$RUN_CONFIG" ]] || { echo "Run configuration not found: $RUN_CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$RUN_CONFIG"
if [[ "${PIPELINE_SKIP_CONDA:-false}" != true ]]; then
  [[ -f "${CONDA_SH:-}" ]] || { echo "Set CONDA_SH to your conda.sh path in the config file" >&2; exit 1; }
  [[ -n "${CONDA_ENV:-}" ]] || { echo "Set CONDA_ENV in the config file" >&2; exit 1; }
  source "$CONDA_SH"
  conda activate "$CONDA_ENV"
fi

cd "$PIPELINE_DIR"
echo "Job started: $(date); host: $(hostname); role: $ROLE"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-not set}; sample index: ${SLURM_ARRAY_TASK_ID:-not an array task}"
python3 "$PIPELINE_DIR/slurm_pipeline.py" "$ROLE" "$MANIFEST"
echo "Job finished: $(date)"
