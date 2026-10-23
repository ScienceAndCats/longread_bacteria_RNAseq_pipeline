#!/usr/bin/env bash
#SBATCH --job-name=bacteria_pipeline
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=bacteria_pipeline-%A_%a.out
#SBATCH --mail-type=END,FAIL
# Optional: uncomment and fill in your address.
##SBATCH --mail-user=YOUR_USERNAME@example.edu

set -euo pipefail

# Submit with: bash submit_bacteria_pipeline.sh config.env
# The helper supplies the array range (up to 40 samples, 8 concurrent by
# default) and the preparation/aggregation dependencies automatically.
ROLE="${1:?Use submit_bacteria_pipeline.sh to submit this template}"
MANIFEST="${2:?Missing run manifest}"
PIPELINE_DIR="${3:?Missing pipeline directory}"

# Fill in your conda installation, or export CONDA_SH before submission.
# For an environment already activated by your site's module/job setup,
# export PIPELINE_SKIP_CONDA=true before submission instead.
if [[ "${PIPELINE_SKIP_CONDA:-false}" != true ]]; then
  CONDA_SH="${CONDA_SH:-/path/to/YOUR_USERNAME/miniforge3/etc/profile.d/conda.sh}"
  [[ -f "$CONDA_SH" ]] || { echo "Set CONDA_SH to your conda.sh path" >&2; exit 1; }
  source "$CONDA_SH"
  conda activate bacteria-RNAseq-pipeline
fi

cd "$PIPELINE_DIR"
echo "Job started: $(date); host: $(hostname); role: $ROLE"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-not set}; sample index: ${SLURM_ARRAY_TASK_ID:-not an array task}"
python3 "$PIPELINE_DIR/slurm_pipeline.py" "$ROLE" "$MANIFEST"
echo "Job finished: $(date)"
