#!/usr/bin/env bash
# Submit all samples plus an automatic dependent aggregation job.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec python3 "$SCRIPT_DIR/slurm_pipeline.py" submit --config "${1:-config.env}"
