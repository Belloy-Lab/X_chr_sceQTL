#!/bin/bash
sex_groups=("both" "female" "male")
cell_types=("Ast" "End" "Exc" "Inh" "Mic" "Oli" "OPC")
mkdir -p logs

R_SCRIPT=/path/to/07_twas/7.5_generate_weights/run_weights.R
SW_DIR=/path/to/07_twas/software

# chmod shared executables once
chmod +x \
  "$SW_DIR/plink/plink" \
  "$SW_DIR/gcta_nr_robust" \
  "$SW_DIR/gemma-0.98.5-linux-static-AMD64"

for sex in "${sex_groups[@]}"; do
  for ct in "${cell_types[@]}"; do
    bsub \
      -G <compute-group> \
      -g <job-group> \
      -q <queue> \
      -sp 100 \
      -R "rusage[mem=32GB] span[hosts=1]" \
      -a 'docker(<docker-image>)' \
      -env "LSF_DOCKER_VOLUMES=/path/to/project:/path/to/project $HOME:$HOME" \
      -o logs/weight_${sex}_${ct}.out \
      -e logs/weight_${sex}_${ct}.err \
      -J weight_${sex}_${ct} \
      Rscript "$R_SCRIPT" "$sex" "$ct"
  done
done