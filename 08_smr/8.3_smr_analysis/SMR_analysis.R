# SMR: each XWAS .ma against the cell-type eQTL BESD
library(data.table)
# CONFIG
smr_bin      <- "smr"   # SMR binary (v1.3.1) placed in this folder; download from https://yanglab.westlake.edu.cn/software/smr/
plink2_bin   <- "plink2"
src_bfile    <- "../../00.1_X_chr_bfile_QC/chrX_QC_EU.0.G"
bfile_23     <- "chrX_QC_EU.23.G"
beqtl_summary <- "../8.1_data_clean_eQTL/eQTL_besd"
gwas_dir     <- "../8.2_xwas_clean"
peqtl_smr    <- 1   # test all probes (SMR default 5e-8)
thread_num   <- 10
# 1. bfile with numeric chromosome (X -> 23)
system(paste("chmod +x", smr_bin))
system(paste(
  plink2_bin,
  "--bfile", src_bfile,
  "--output-chr 26",
  "--make-bed",
  "--out", bfile_23
))
# 2. One job per .ma
gwas_files <- dir(gwas_dir, pattern = "\\.ma$", full.names = TRUE)
analysis_ref <- data.table(gwas_path = gwas_files)
analysis_ref[, trait      := gsub("\\.ma$", "", basename(gwas_path))]
analysis_ref[, out_prefix := paste0("smr_res_", trait)]
print(analysis_ref)
# 3. Run SMR
for (i in seq_len(nrow(analysis_ref))) {
  gwas <- analysis_ref$gwas_path[i]
  out  <- analysis_ref$out_prefix[i]
  
  cmd <- paste(
    paste0("./", smr_bin),
    "--bfile", bfile_23,
    "--gwas-summary", gwas,
    "--beqtl-summary", beqtl_summary,
    "--peqtl-smr", peqtl_smr,
    "--out", out,
    "--thread-num", thread_num
  )
  
  cat("\n=============================\n")
  cat("Running:", out, "\n")
  cat("=============================\n")
  system(cmd)
}