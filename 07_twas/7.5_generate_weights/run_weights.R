# FUSION weights for one sex x cell type (non-PAR genes); usage: Rscript run_weights.R <sex> <cell_type>
args <- commandArgs(trailingOnly = TRUE)
sex_input      <- args[1]
celltype_input <- args[2]
cat("Running:", sex_input, celltype_input, "\n")
library(data.table)
# CONFIG
bfile_base <- "../7.2_extract_bfile/"
pheno_base <- "../7.4_extract_pheno/"
covar_base <- "../7.3_extract_covariate/"
sw_dir        <- "../software/"
plink_bin     <- paste0(sw_dir, "plink/plink")
gcta_bin      <- paste0(sw_dir, "gcta_nr_robust")
gemma_bin     <- paste0(sw_dir, "gemma-0.98.5-linux-static-AMD64")
fusion_script <- paste0(sw_dir, "fusion_twas-master/FUSION.compute_weights_XCHR.R")
r_libs_user   <- "path/to/R_library"
# 1. Reference table for this sex x cell type
combo_dir <- file.path(bfile_base, sex_input, celltype_input)
bed_files <- list.files(combo_dir, pattern = "\\.bed$")
genes <- sort(unique(sub(paste0("_", sex_input, "_", celltype_input, "\\.bed$"), "", bed_files)))

if (length(genes) == 0) {
  stop("No .bed files found in ", combo_dir)
}

df_reference <- data.table(
  sex       = sex_input,
  cell_type = celltype_input,
  gene      = genes,
  PAR_type  = "non_PAR"
)

file_tag <- paste0(df_reference$gene, "_", df_reference$sex, "_", df_reference$cell_type)
df_reference$bfile <- paste0(bfile_base, df_reference$sex, "/", df_reference$cell_type, "/", file_tag)
df_reference$pheno <- paste0(pheno_base, df_reference$sex, "/", df_reference$cell_type, "/", file_tag, ".pheno")
df_reference$covar <- paste0(covar_base, df_reference$sex, "/", df_reference$cell_type, "/",
                             df_reference$sex, "_", df_reference$cell_type, ".covar")
# 2. Output paths
dir.create("tmp", showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("weights", sex_input, celltype_input), showWarnings = FALSE, recursive = TRUE)

df_reference$tmp <- paste0("tmp/", df_reference$sex, "_", df_reference$cell_type, "_", df_reference$gene)
df_reference$out <- paste0("weights/", df_reference$sex, "/", df_reference$cell_type, "/", df_reference$gene)

write.csv(df_reference,
          paste0("job_reference_table_", sex_input, "_", celltype_input, ".csv"),
          row.names = FALSE)
# 3. Ensure executables are runnable (also done once in submit_all.sh)
for (bin in c(plink_bin, gcta_bin, gemma_bin)) {
  if (file.access(bin, mode = 1) != 0) system(paste("chmod +x", shQuote(bin)))
}
# 4. Compute weights per gene
for (i in seq_len(nrow(df_reference))) {
  cmd <- paste(
    paste0("R_LIBS_USER=", r_libs_user, " Rscript"),
    fusion_script,
    "--save_hsq TRUE",
    "--bfile", df_reference$bfile[i],
    "--pheno", df_reference$pheno[i],
    "--covar", df_reference$covar[i],
    "--tmp",   df_reference$tmp[i],
    "--out",   df_reference$out[i],
    "--hsq_p 1 --verbose 1",
    "--noclean TRUE",
    "--PATH_plink", plink_bin,
    "--PATH_gemma", gemma_bin,
    "--PATH_gcta",  gcta_bin
  )
  cat("\n=== Running index:", i, "===\n")
  cat(cmd, "\n\n")
  system(cmd)
}