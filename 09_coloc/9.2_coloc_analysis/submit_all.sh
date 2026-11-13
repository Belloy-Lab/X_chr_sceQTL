#!/bin/bash
set -euo pipefail
mkdir -p logs
export LSF_DOCKER_VOLUMES="/path/to/project:/path/to/project $HOME:$HOME"
RSCRIPT="/path/to/09_coloc/9.2_coloc_analysis/coloc.R"

# master
if [[ -z "${JOB_INDEX:-}" ]]; then

  N_JOBS=$(($(wc -l < job_table.txt) - 1))
  echo "Submitting ${N_JOBS} jobs..."

  for ((IDX=1; IDX<=N_JOBS; IDX++)); do
    bsub \
      -G <compute-group> \
      -g <job-group> \
      -q <queue> \
      -sp 100 \
      -R "rusage[mem=16GB] span[hosts=1]" \
      -a 'docker(<docker-image>)' \
      -env "all,JOB_INDEX=${IDX}" \
      -J "coloc_${IDX}" \
      -o logs/%J.out \
      -e logs/%J.err \
      "bash $0"
  done
  exit 0
fi

# worker
IDX=${JOB_INDEX}
echo "Running job index: $IDX"
Rscript "$RSCRIPT" "$IDX"