#!/bin/bash
set -euo pipefail

JOB_TABLE="TWAS_job_reference_LBD.txt"
FUSION_SCRIPT="../software/fusion_twas-master/FUSION.assoc_test.R"

TOTAL=$(($(wc -l < "$JOB_TABLE") - 1))
COUNT=0
START_TIME=$(date +%s)

echo "Total jobs: $TOTAL"
echo "====================================="

tail -n +2 "$JOB_TABLE" | while IFS=$'\t' read -r \
  PAR_type trait sex cell_type gwas_model \
  sumstats pos_path weights_dir REF_LD_CHR chr \
  out_dir out_prefix out_file log_file
do
  COUNT=$((COUNT+1))

  NOW=$(date +%s)
  ELAPSED=$((NOW - START_TIME))
  AVG=$((ELAPSED / COUNT))
  REMAIN=$((AVG * (TOTAL - COUNT)))
  PERCENT=$((COUNT * 100 / TOTAL))

  echo "[$COUNT/$TOTAL | ${PERCENT}%] $trait $sex $cell_type $gwas_model chr$chr | ETA: ${REMAIN}s"

  mkdir -p "$out_dir"
  mkdir -p "$(dirname "$log_file")"

  if [[ -f "$out_file" ]]; then
    echo "[SKIP] Output exists: $out_file"
    continue
  fi

  Rscript "$FUSION_SCRIPT" \
    --sumstats "$sumstats" \
    --weights "$pos_path" \
    --weights_dir "$weights_dir" \
    --ref_ld_chr "$REF_LD_CHR" \
    --chr "$chr" \
    --out "$out_file" \
    > "$log_file" 2>&1

done

echo "====================================="
echo "All jobs completed."