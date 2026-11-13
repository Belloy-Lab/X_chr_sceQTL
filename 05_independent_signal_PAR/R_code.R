# Independent PAR signals: combine clumped significant eQTLs across settings (wide), sex-heterogeneity test, second LD clumping across settings
rm(list = ls())
library(data.table)
library(dplyr)
library(tidyr)
library(ieugwasr)
# CONFIG
cis_all_rds <- "../03_merge_eQTL_res_PAR/df_eqtl_cis_all.rds"
rsid_annotation_rds <- "../00.3_rsid_annotation/rsid_gene_annotation.rds"
afreq_prefix      <- "../00.2_PAR_chr_bfile_QC/chrX_QC_EU.0.G"
afreq_both_path   <- paste0(afreq_prefix, ".afreq")
afreq_female_path <- paste0(afreq_prefix, ".female.afreq")
afreq_male_path   <- paste0(afreq_prefix, ".male.afreq")
clump_dir        <- "../04_significant_eQTL_PAR/"
clump_both_rds   <- paste0(clump_dir, "eqtl_sig_clumped_all_PAR_chr_2stepFDR_after_prune_both_par.rds")
clump_female_rds <- paste0(clump_dir, "eqtl_sig_clumped_all_PAR_chr_2stepFDR_after_prune_female_par.rds")
clump_male_rds   <- paste0(clump_dir, "eqtl_sig_clumped_all_PAR_chr_2stepFDR_after_prune_male_par.rds")
ld_ref_bfile <- "../00.2_PAR_chr_bfile_QC/chrX_QC_EU.0.G"
# 1. Annotate cis-eQTL results with rsID, gene, allele frequency
df_eqtl_cis_all <- readRDS(cis_all_rds)
df_eqtl_cis_all$fdr <- NULL
df_rsid_gene_annotation <- readRDS(rsid_annotation_rds)
df_eqtl_cis_all <- left_join(df_eqtl_cis_all,df_rsid_gene_annotation,by='ID')
afreq_both <- fread(afreq_both_path)
afreq_both$OBS_CT <- NULL
colnames(afreq_both)[5] <- "EAF_both"
afreq_female <- fread(afreq_female_path)
afreq_female$OBS_CT <- NULL
colnames(afreq_female)[5] <- "EAF_female"
afreq_male <- fread(afreq_male_path)
afreq_male$OBS_CT <- NULL
colnames(afreq_male)[5] <- "EAF_male"
df_eqtl_cis_all <- left_join(df_eqtl_cis_all,afreq_both,by=c("#CHROM","ID", "REF", "ALT"))
df_eqtl_cis_all <- left_join(df_eqtl_cis_all,afreq_female,by=c("#CHROM","ID", "REF", "ALT"))
df_eqtl_cis_all <- left_join(df_eqtl_cis_all,afreq_male,by=c("#CHROM","ID", "REF", "ALT"))
# 2. Keep cell-type-gene-variant pairs significant in >= 1 setting; pivot wide across settings
df_clump_sig_both <- readRDS(clump_both_rds)
df_clump_sig_female_model_2 <- readRDS(clump_female_rds)
df_clump_sig_male_model_2 <- readRDS(clump_male_rds)
df_clump_sig_all_model <- rbind(df_clump_sig_both,
                                df_clump_sig_female_model_2,
                                df_clump_sig_male_model_2)
df_eqtl_cis_clump_sig_all_model <- semi_join(df_eqtl_cis_all,df_clump_sig_all_model,
                                             by=c("ID","cell_type","gene"))
df_eqtl_cis_clump_sig_all_model <- df_eqtl_cis_clump_sig_all_model %>%
  left_join(
    df_clump_sig_all_model %>%
      select(ID, cell_type, gene, setting, significant_by_2step_FDR),
    by = c("ID", "cell_type", "gene", "setting")
  )
df_eqtl_cis_clump_sig_all_model$TEST <- NULL
df_eqtl_cis_clump_sig_all_model$ERRCODE <- NULL
# A1 should equal ALT
identical(df_eqtl_cis_clump_sig_all_model$A1,df_eqtl_cis_clump_sig_all_model$ALT)
df_eqtl_cis_clump_sig_all_model$A1 <- NULL
id_cols <- c(
  "#CHROM", "POS", "ID", "rsid","REF", "ALT",
  "cell_type", "gene", "par_type",
  "EAF_both", "EAF_female", "EAF_male"
)
value_cols <- c(
  "OBS_CT", "BETA", "SE", "T_STAT", "LOG10_P", "pvalue",
  "significant_by_2step_FDR"
)
df_eqtl_wide <- df_eqtl_cis_clump_sig_all_model %>%
  pivot_wider(
    id_cols   = all_of(id_cols),
    names_from  = setting,
    values_from = all_of(value_cols),
    names_glue  = "{.value}_{setting}"
  )
# 3. Min p-value across settings, male/female beta ratio, sex-heterogeneity test
p_cols <- grep("^pvalue_", colnames(df_eqtl_wide), value = TRUE)
df_eqtl_wide$min_pvalue <- apply(
  df_eqtl_wide[, p_cols],
  1,
  function(x) min(x, na.rm = TRUE)
)
df_eqtl_wide$min_pvalue_model <- apply(
  df_eqtl_wide[, p_cols],
  1,
  function(x) {
    p_cols[which.min(x)]
  }
)
df_eqtl_wide$min_pvalue_model <- sub("pvalue_", "", df_eqtl_wide$min_pvalue_model)
df_eqtl_wide$male_female_ratio <- (df_eqtl_wide$BETA_male_par)/(df_eqtl_wide$BETA_female_par)
df_eqtl_wide$abs_male_female_ratio <- abs(df_eqtl_wide$male_female_ratio)
eps <- 1e-6
df_eqtl_wide <- df_eqtl_wide %>%
  mutate(
    Z_het = (BETA_male_par - BETA_female_par) /
      sqrt(SE_male_par^2 + SE_female_par^2),
    P_het = 2 * pnorm(-abs(Z_het)),
    # BH-FDR is computed in step 4b on the final independent lead variants
    P_het_FDR = NA_real_,
    
    male_sign   = sign(ifelse(abs(BETA_male_par) < eps, 0, BETA_male_par)),
    female_sign = sign(ifelse(abs(BETA_female_par) < eps, 0, BETA_female_par)),
    opposite_dir = (male_sign * female_sign) == -1,
    
    delta_abs = abs(BETA_male_par) -
      abs(BETA_female_par)
  ) %>%
  mutate(
    sex_het_label_nominal = case_when(
      P_het >= 0.05 ~ "No sex heterogeneity",
      
      !opposite_dir & delta_abs > 0 ~ "Male-biased (concordant)",
      !opposite_dir & delta_abs < 0 ~ "Female-biased (concordant)",
      
      opposite_dir & delta_abs > 0 ~ "Male-dominant (opposite)",
      opposite_dir & delta_abs < 0 ~ "Female-dominant (opposite)"
    ))
saveRDS(df_eqtl_wide,"df_eqtl_wide_before_second_clump.rds")
# 4. Second LD clumping across settings (ld_clump_local patched with --allow-extra-chr for PAR)
ld_clump_local_new <- function(dat, clump_kb, clump_r2, clump_p, bfile, plink_bin)
{
  shell <- ifelse(Sys.info()["sysname"] == "Windows", "cmd", "sh")
  fn <- tempfile()
  
  write.table(
    data.frame(SNP = dat$rsid, P = dat$pval),
    file = fn,
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  
  cmd <- paste0(
    shQuote(plink_bin, type = shell),
    " --bfile ", shQuote(bfile, type = shell),
    " --allow-extra-chr",
    " --clump ", shQuote(fn, type = shell),
    " --clump-p1 ", clump_p,
    " --clump-r2 ", clump_r2,
    " --clump-kb ", clump_kb,
    " --out ", shQuote(fn, type = shell)
  )
  
  system(cmd)
  
  res <- read.table(paste0(fn, ".clumped"), header = TRUE)
  unlink(paste0(fn, "*"))
  
  dat[dat$rsid %in% res$SNP, ]
}
ns <- asNamespace("ieugwasr")
unlockBinding("ld_clump_local", ns)
assign("ld_clump_local", ld_clump_local_new, envir = ns)
lockBinding("ld_clump_local", ns)
column_order <- colnames(df_eqtl_wide)
df_eqtl_wide$snp <- df_eqtl_wide$rsid
df_eqtl_wide$rsid <- df_eqtl_wide$ID
df_eqtl_wide$pval <- df_eqtl_wide$min_pvalue
df_eqtl_wide$id <- paste0(df_eqtl_wide$gene,"_",df_eqtl_wide$cell_type)
df_eqtl_wide <- ld_clump(
  df_eqtl_wide,
  clump_kb  = 1000,
  clump_r2  = 0.1,
  clump_p   = 1,
  bfile     = ld_ref_bfile,
  plink_bin = plinkbinr::get_plink_exe()
)
df_eqtl_wide$rsid <- df_eqtl_wide$snp
df_eqtl_wide$snp <- NULL
df_eqtl_wide$pval <- NULL
df_eqtl_wide$id <- NULL
df_eqtl_wide <- select(df_eqtl_wide,column_order)
# 4b. Sex-heterogeneity BH-FDR on the final independent lead variants (PAR tested separately from non-PAR)
het_idx <- !is.na(df_eqtl_wide$P_het)
df_eqtl_wide$P_het_FDR[het_idx] <- p.adjust(df_eqtl_wide$P_het[het_idx], method = "BH")
saveRDS(df_eqtl_wide,"df_eqtl_wide_after_second_clump.rds")
