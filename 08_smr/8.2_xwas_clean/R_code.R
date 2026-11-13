# XWAS summary stats -> SMR .ma files (one per trait x model), variant IDs aligned to the chrX bfile
rm(list = ls())
library(dplyr)
library(data.table)
# CONFIG
bim_file  <- "path/to/chrX_reference.bim"
gwas_root <- "path/to/XWAS"
# 1. bfile variant IDs
bfile_bim <- fread(bim_file)
bfile_id  <- data.frame(ID = bfile_bim$V2)
rm(bfile_bim)
# 2. Column schemas and alignment function
schema_cleaned <- list(EA = "ALLELE1", NEA = "ALLELE0", BP = "BP",
                       BETA = "BETA", SE = "SE", P = "P", N = "N_incl",
                       FREQ = "A1FREQ")
schema_ebi     <- list(EA = "effect_allele", NEA = "other_allele", BP = "base_pair_location",
                       BETA = "beta", SE = "standard_error", P = "p_value", N = "N_incl",
                       FREQ = "effect_allele_frequency")
schemas <- list(cleaned = schema_cleaned, ebi = schema_ebi)
# allele-reversed matches: swap A1/A2, negate beta, freq = 1 - freq
prepare_smr_sumstats <- function(df_GWAS, bfile_id,
                                 EA_col, NEA_col, BP_col,
                                 BETA_col, SE_col, P_col, N_col, FREQ_col,
                                 chr_prefix = "chrX") {
  df <- copy(df_GWAS)
  df$id1 <- paste0(chr_prefix, ":", df[[BP_col]], ":", df[[NEA_col]], ":", df[[EA_col]])
  df$id2 <- paste0(chr_prefix, ":", df[[BP_col]], ":", df[[EA_col]], ":", df[[NEA_col]])
  df$match <- ifelse(df$id1 %in% bfile_id$ID, "id1",
                     ifelse(df$id2 %in% bfile_id$ID, "id2", NA))
  df <- df[!is.na(df$match), ]
  df$SNP <- ifelse(df$match == "id1", df$id1, df$id2)
  
  df_out <- df %>%
    transmute(
      SNP  = SNP,
      A1   = ifelse(match == "id1", .data[[EA_col]],  .data[[NEA_col]]),
      A2   = ifelse(match == "id1", .data[[NEA_col]], .data[[EA_col]]),
      freq = ifelse(match == "id1", .data[[FREQ_col]], 1 - .data[[FREQ_col]]),
      b    = ifelse(match == "id1", .data[[BETA_col]], -.data[[BETA_col]]),
      se   = .data[[SE_col]],
      p    = .data[[P_col]],
      n    = .data[[N_col]]
    ) %>%
    na.omit()
  
  df_out <- dplyr::select(df_out, SNP, A1, A2, freq, b, se, p, n)
  df_out
}
# 3. Job table; n_manual = sample size when the file has no N_incl column
jobs <- data.frame(
  name = c("AD_rXCI", "AD_eXCI", "AD_male", "AD_female",
           "MSA_rXCI", "MSA_eXCI", "MSA_male", "MSA_female",
           "PD_rXCI", "PD_eXCI", "PD_male", "PD_female",
           "LBD_rXCI", "LBD_eXCI", "LBD_male", "LBD_female"),
  rel_path = c(
    "AD/XWAS_rXCI.txt",
    "AD/XWAS_eXCI.txt",
    "AD/XWAS_male.txt",
    "AD/XWAS_female.txt",
    "MSA/XWAS_rXCI.txt",
    "MSA/XWAS_eXCI.txt",
    "MSA/XWAS_male.txt",
    "MSA/XWAS_female.txt",
    "PD/XWAS_rXCI.txt",
    "PD/XWAS_eXCI.txt",
    "PD/XWAS_male.txt",
    "PD/XWAS_female.txt",
    "LBD/XWAS_rXCI.txt",
    "LBD/XWAS_eXCI.txt",
    "LBD/XWAS_male.txt",
    "LBD/XWAS_female.txt"
  ),
  schema = c("cleaned", "cleaned", "cleaned", "cleaned",
             "ebi",     "cleaned", "ebi",     "ebi",
             "cleaned", "cleaned", "cleaned", "cleaned",
             "ebi",     "cleaned", "ebi",     "ebi"),
  n_manual = c(NA, NA, NA, NA,
               8016, NA, 3859, 4157,
               NA, NA, NA, NA,
               6614, NA, 3606, 3008),
  stringsAsFactors = FALSE
)
# 4. Build .ma files
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, ]
  df_GWAS <- fread(file.path(gwas_root, job$rel_path))
  if (!is.na(job$n_manual)) df_GWAS$N_incl <- job$n_manual
  
  sc <- schemas[[job$schema]]
  df_SMR <- prepare_smr_sumstats(
    df_GWAS  = df_GWAS,
    bfile_id = bfile_id,
    EA_col   = sc$EA,   NEA_col = sc$NEA, BP_col   = sc$BP,
    BETA_col = sc$BETA, SE_col  = sc$SE,  P_col    = sc$P,
    N_col    = sc$N,    FREQ_col = sc$FREQ
  )
  
  fwrite(df_SMR, paste0(job$name, ".ma"), sep = "\t", quote = FALSE)
  cat(sprintf("[%d/%d] %s -> %s.ma  (%d variants)\n",
              i, nrow(jobs), job$name, job$name, nrow(df_SMR)))
}
