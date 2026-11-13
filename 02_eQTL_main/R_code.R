# eQTL mapping with plink2 --glm for 7 settings (PAR both/female/male; X eXCI/rXCI/female/male) x 7 cell types
rm(list = ls())
library(data.table)
# CONFIG
cell_types <- c("Ast", "End", "Exc", "Inh", "Mic", "Oli", "OPC")
cov_pheno_dir <- "../01.2_cov_pheno/"
valid_id_dir  <- "../01.1_valid_sample_id/"
bfile_dir_X   <- "../00.1_X_chr_bfile_QC/chrX_QC_EU.0.G"
bfile_dir_PAR <- "../00.2_PAR_chr_bfile_QC/chrX_QC_EU.0.G"
# fixed covariate list for PAR; X uses all non-sex columns of .cov.txt
par_covar_name <- "age_death,pmi,ROS_study,Total_Genes_Detected,gPC1,gPC2,gPC3,ePC1,ePC2,ePC3,ePC4,ePC5,ePC6,ePC7,ePC8,ePC9,ePC10,ePC11,ePC12,ePC13,ePC14,ePC15,ePC16,ePC17,ePC18,ePC19,ePC20,ePC21,ePC22,ePC23,ePC24,ePC25,ePC26,ePC27,ePC28,ePC29,ePC30"
# model_type (--xchr-model): 1 = eXCI, 2 = rXCI; NA for PAR
analysis_grid <- data.frame(
  folder     = c("2_eQTL_PAR_chr_Both", "2_eQTL_PAR_chr_Female", "2_eQTL_PAR_chr_Male",
                 "2_eQTL_X_chr_eXCI",   "2_eQTL_X_chr_Female",   "2_eQTL_X_chr_Male",
                 "2_eQTL_X_chr_rXCI"),
  chr_type   = c("PAR", "PAR", "PAR", "X", "X", "X", "X"),
  sex        = c("both", "female", "male", "both", "female", "male", "both"),
  model_type = c(NA, NA, NA, 1, 2, 1, 2),
  stringsAsFactors = FALSE
)
for (a in seq_len(nrow(analysis_grid))) {
  
  folder     <- analysis_grid$folder[a]
  chr_type   <- analysis_grid$chr_type[a]
  sex        <- analysis_grid$sex[a]
  model_type <- analysis_grid$model_type[a]
  
  cat("\n==============================================================\n")
  cat("Analysis", a, "/", nrow(analysis_grid), ":", folder,
      "| chr =", chr_type, "| sex =", sex,
      if (!is.na(model_type)) paste("| model =", model_type) else "", "\n")
  cat("==============================================================\n")
  
  bfile_dir <- if (chr_type == "X") bfile_dir_X else bfile_dir_PAR
  
  eqtl_parameter_df <- data.frame(cell_type = cell_types, stringsAsFactors = FALSE)
  eqtl_parameter_df$cov_dir <- paste0(cov_pheno_dir, sex, "_",
                                      eqtl_parameter_df$cell_type, ".cov.txt")
  eqtl_parameter_df$pheno_dir <- paste0(cov_pheno_dir, sex, "_",
                                        eqtl_parameter_df$cell_type, ".pheno_chrX.tsv")
  eqtl_parameter_df$valid_sample_id <- paste0(valid_id_dir, "valid_sample_id_",
                                              sex, "_", eqtl_parameter_df$cell_type, ".txt")
  eqtl_parameter_df$bfile_dir <- bfile_dir
  eqtl_parameter_df$out_dir <- paste0(getwd(), "/", folder, "/", eqtl_parameter_df$cell_type)
  
  dir.create(folder, showWarnings = FALSE)
  for (d in eqtl_parameter_df$cell_type) {
    dir.create(file.path(folder, d), showWarnings = FALSE, recursive = TRUE)
  }
  
  glm_part <- if (sex == "both") {
    "--glm sex hide-covar omit-ref log10"
  } else {
    "--glm hide-covar omit-ref log10"
  }
  
  for (i in seq_len(nrow(eqtl_parameter_df))) {
    out_prefix <- file.path(
      eqtl_parameter_df$out_dir[i],
      basename(eqtl_parameter_df$out_dir[i])
    )
    
    if (chr_type == "PAR") {
      covar_name <- par_covar_name
    } else {
      cov_df <- fread(eqtl_parameter_df$cov_dir[i])
      covar_cols <- colnames(cov_df)
      covar_cols <- covar_cols[!covar_cols %in% c("msex", "FID", "IID")]
      covar_name <- paste(covar_cols, collapse = ",")
    }
    
    if (chr_type == "X") {
      cmd <- paste(
        "plink2",
        paste0("--bfile ", eqtl_parameter_df$bfile_dir[i]),
        paste0("--keep ", eqtl_parameter_df$valid_sample_id[i]),
        paste0("--covar ", eqtl_parameter_df$cov_dir[i]),
        paste0("--covar-name ", covar_name),
        "--covar-variance-standardize",
        paste0("--chr ", chr_type),
        paste0("--xchr-model ", model_type),
        glm_part,
        paste0("--pheno ", eqtl_parameter_df$pheno_dir[i]),
        paste0("--out ", out_prefix)
      )
    }
    
    if (chr_type == "PAR") {
      cmd <- paste(
        "plink2",
        paste0("--bfile ", eqtl_parameter_df$bfile_dir[i]),
        paste0("--keep ", eqtl_parameter_df$valid_sample_id[i]),
        paste0("--covar ", eqtl_parameter_df$cov_dir[i]),
        paste0("--covar-name ", covar_name),
        "--covar-variance-standardize",
        "--chr PAR1,PAR2",
        glm_part,
        paste0("--pheno ", eqtl_parameter_df$pheno_dir[i]),
        paste0("--out ", out_prefix)
      )
    }
    
    system(cmd)
  }
}