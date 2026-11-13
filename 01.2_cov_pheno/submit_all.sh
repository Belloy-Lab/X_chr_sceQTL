#!/bin/bash
set -euo pipefail
mkdir -p logs

# CONFIG (LSF)
DOCKER_VOLUMES="/path/to/project:/path/to/project $HOME:$HOME"
RSCRIPT="/path/to/01.2_cov_pheno/R_code_job_submit.R"
COMPUTE_GROUP="<compute-group>"
JOB_GROUP="<job-group>"
QUEUE="<queue>"
DOCKER_IMAGE="<docker-image>"
N_JOBS=21   # one per (sex x cell type)

export LSF_DOCKER_VOLUMES="$DOCKER_VOLUMES"

# master
if [[ -z "${JOB_INDEX:-}" ]]; then
  echo "Submitting ${N_JOBS} jobs..."

  for IDX in $(seq 1 "${N_JOBS}"); do
    bsub \
      -G "${COMPUTE_GROUP}" \
      -g "${JOB_GROUP}" \
      -q "${QUEUE}" \
      -R "rusage[mem=100GB]" \
      -a "docker(${DOCKER_IMAGE})" \
      -env "all,JOB_INDEX=${IDX}" \
      -J "eQTL_${IDX}" \
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