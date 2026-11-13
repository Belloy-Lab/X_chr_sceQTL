# Per cell type x sex stratum: sample IDs present in the snRNA covariate file and passing chrX QC
rm(list = ls())
library(data.table)
library(dplyr)
# CONFIG
pseudobulk_dir <- "path/to/pseudobulk/"
fam_pass_qc    <- "../00.1_X_chr_bfile_QC/chrX_QC_EU.0.G.fam"   # V5: 1 = male, 2 = female
sex_strata <- list(both = NA, male = 1, female = 2)
cov_files  <- dir(pseudobulk_dir, pattern = ".covariates.tsv")
cell_types <- sub("\\.covariates\\.tsv$", "", cov_files)
cov_files  <- paste0(pseudobulk_dir, cov_files)
df_fam_pass_qc <- fread(fam_pass_qc)
for (sex in names(sex_strata)) {
  sex_code <- sex_strata[[sex]]
  keep_ids <- if (is.na(sex_code)) {
    df_fam_pass_qc$V2
  } else {
    df_fam_pass_qc$V2[df_fam_pass_qc$V5 == sex_code]
  }
  
  for (i in seq_along(cell_types)) {
    ct <- cell_types[i]
    cat("Extract valid samples (", sex, ") from:", ct, "\n")
    df_cov <- fread(cov_files[i])
    df_cov$id <- NULL
    df_sample_plink_format <- data.frame(FID = colnames(df_cov),
                                         IID = colnames(df_cov))
    df_sample_plink_format <- filter(df_sample_plink_format, IID %in% keep_ids)
    cat(nrow(df_sample_plink_format), sex,
        "samples were kept (overlapped with passed QC samples)", "\n")
    fwrite(df_sample_plink_format, paste0("valid_sample_id_", sex, "_", ct, ".txt"),
           sep = "\t", quote = FALSE, row.names = FALSE)
  }
}