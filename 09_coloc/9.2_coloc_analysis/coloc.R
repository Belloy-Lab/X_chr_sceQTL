# coloc.abf for one job (XWAS x cell type x setting), over all genes; usage: Rscript coloc.R <job_id>
library(data.table)
library(dplyr)
library(coloc)
# CONFIG
eqtl_rds       <- "../../03_merge_eQTL_res_X/df_eqtl_cis_all.rds"
afreq_file     <- "../../00.1_X_chr_bfile_QC/chrX_QC_EU.0.G.afreq"
gwas_dir       <- "../9.1_xwas_clean/"
job_table_file <- "job_table.txt"
# 1. cis-eQTL with EAF; XWAS list (order must match job_table)
df_eQTL_all      <- readRDS(eqtl_rds)
df_eQTL_distinct <- distinct(df_eQTL_all, cell_type, gene, setting)
cell_types       <- unique(df_eQTL_distinct$cell_type)
df_afreq <- fread(afreq_file)
df_afreq$OBS_CT <- NULL
df_eQTL_all <- left_join(df_eQTL_all, df_afreq)
GWAS_list <- dir(gwas_dir, pattern = ".ma")
GWAS_list <- data.frame(trait     = substr(GWAS_list, 1, nchar(GWAS_list) - 3),
                        GWAS_path = paste0(gwas_dir, GWAS_list))
GWAS_list$GWAS_type <- sub(".*_", "", GWAS_list$trait)
# 2. This job
args   <- commandArgs(trailingOnly = TRUE)
job_id <- as.numeric(args[1])
job_table <- fread(job_table_file)
job <- job_table[job_id]
i           <- job$gwas_i
cell_type_i <- job$cell_type
setting_i   <- job$setting
cat("JOB:", job_id, "GWAS:", GWAS_list$trait[i],
    "cell:", cell_type_i, "setting:", setting_i, "\n")
# 3. XWAS, duplicated SNPs dropped
df_GWAS_all <- fread(GWAS_list$GWAS_path[i]) %>%
  group_by(SNP) %>%
  filter(n() == 1) %>%
  ungroup()
# 4. coloc per gene
all_coloc_results_i   <- list()
all_variant_results_i <- list()
k <- 1
v <- 1
df_eqtl_sub <- df_eQTL_all %>%
  filter(cell_type == cell_type_i,
         setting   == setting_i)
genes <- unique(df_eqtl_sub$gene)
for (g in genes) {
  
  df_eqtl <- df_eqtl_sub %>% filter(gene == g)
  
  df_gwas <- semi_join(df_GWAS_all, df_eqtl, by = c("SNP" = "ID"))
  df_eqtl <- semi_join(df_eqtl, df_gwas, by = c("ID" = "SNP"))
  
  if (nrow(df_gwas) == 0) next
  
  df_gwas <- arrange(df_gwas, SNP)
  df_eqtl <- arrange(df_eqtl, ID)
  
  # same variants, order and A1 required
  if (!identical(df_gwas$SNP, df_eqtl$ID)) next
  if (!all(df_gwas$A1 == df_eqtl$A1))      next
  
  df_gwas$BP <- as.integer(sub(".*:(\\d+):.*", "\\1", df_gwas$SNP))
  
  D1 <- list(
    pvalues  = df_gwas$p,
    type     = "cc",
    snp      = df_gwas$SNP,
    position = df_gwas$BP,
    N        = max(df_gwas$n),
    MAF      = pmin(df_gwas$freq, 1 - df_gwas$freq),
    s        = df_gwas$s[1]
  )
  D2 <- list(
    pvalues  = df_eqtl$pvalue,
    type     = "quant",
    snp      = df_eqtl$ID,
    position = df_eqtl$POS,
    N        = df_eqtl$OBS_CT,
    MAF      = pmin(df_eqtl$ALT_FREQS, 1 - df_eqtl$ALT_FREQS)
  )
  
  coloc_results <- coloc.abf(dataset1 = D1, dataset2 = D2)
  
  all_coloc_results_i[[k]] <- data.frame(
    job_id    = job_id,
    gene      = g,
    cell_type = cell_type_i,
    eQTL_type = setting_i,
    trait     = GWAS_list$trait[i],
    GWAS_type = GWAS_list$GWAS_type[i],
    PP0    = coloc_results$summary["PP.H0.abf"],
    PP1    = coloc_results$summary["PP.H1.abf"],
    PP2    = coloc_results$summary["PP.H2.abf"],
    PP3    = coloc_results$summary["PP.H3.abf"],
    PP4    = coloc_results$summary["PP.H4.abf"],
    n_snps = nrow(df_gwas)
  )
  k <- k + 1
  
  coloc_variant_res <- coloc_results[["results"]]
  coloc_variant_res$job_id    <- job_id
  coloc_variant_res$gene      <- g
  coloc_variant_res$cell_type <- cell_type_i
  coloc_variant_res$eQTL_type <- setting_i
  coloc_variant_res$trait     <- GWAS_list$trait[i]
  coloc_variant_res$GWAS_type <- GWAS_list$GWAS_type[i]
  all_variant_results_i[[v]] <- coloc_variant_res
  v <- v + 1
}
# 5. Save
summary_df <- bind_rows(all_coloc_results_i)
variant_df <- bind_rows(all_variant_results_i)
fwrite(summary_df, paste0("coloc_summary_job_", job_id, ".txt"), sep = "\t")
fwrite(variant_df, paste0("coloc_variant_job_", job_id, ".txt"), sep = "\t")
cat("Saved JOB:", job_id, "\n")