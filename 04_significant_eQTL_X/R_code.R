# Significant non-PAR X cis-eQTLs: two-step FDR, then LD clumping (r2 < 0.1) for independent top signals
rm(list = ls())
library(data.table)
library(dplyr)
library(ieugwasr)
# CONFIG
cis_all_rds <- "../03_merge_eQTL_res_X/df_eqtl_cis_all.rds"
rsid_annotation_rds <- "../00.3_rsid_annotation/rsid_gene_annotation.rds"
# QC bfile (also LD reference) and allele-frequency files
bfile_prefix      <- "../00.1_X_chr_bfile_QC/chrX_QC_EU.0.G"
afreq_both_path   <- paste0(bfile_prefix, ".afreq")
afreq_female_path <- paste0(bfile_prefix, ".female.afreq")
afreq_male_path   <- paste0(bfile_prefix, ".male.afreq")
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
# 2. Two-step FDR (gene-wise Bonferroni, then BH across genes) per cell type x setting
cell_type_setting_combination <- distinct(df_eqtl_cis_all, cell_type, setting)
df_eqtl_cis_sort_list <- vector(
  mode = "list",
  length = nrow(cell_type_setting_combination)
)
cell_type_setting_combination$n_gene_included <- NA
cell_type_setting_combination$n_gene_with_cis_eqtl <- NA
for (i in 1:nrow(cell_type_setting_combination)) {
  ct <- cell_type_setting_combination$cell_type[i]
  st <- cell_type_setting_combination$setting[i]
  cat("Identifying significant eQTL singals | cell type:",ct,"| setting:",st,"\n")
  
  df_eqtl_cis <- filter(df_eqtl_cis_all, cell_type == ct & setting == st)
  
  df_eqtl_cis <- df_eqtl_cis %>%
    group_by(gene) %>%
    mutate(p.adjusted_per_gene = p.adjust(pvalue, method = "bonferroni")) %>%
    ungroup()
  
  fdr <- df_eqtl_cis %>%
    group_by(gene) %>%
    summarize(p.gene_level = min(p.adjusted_per_gene)) %>%
    ungroup() %>%
    mutate(
      BH_fdr = p.adjust(p.gene_level, method = "BH")
    )
  
  thr.BH_fdr <- fdr %>% filter(BH_fdr < 0.05) %>% .$p.gene_level %>% max
  
  df_eqtl_cis <- df_eqtl_cis %>%
    mutate(
      significant_by_2step_FDR = ifelse(
        p.adjusted_per_gene <= thr.BH_fdr,
        "Yes",
        "No"
      )
    )
  
  setDT(df_eqtl_cis)
  
  n_unique_gene <- uniqueN(
    df_eqtl_cis[significant_by_2step_FDR == "Yes", gene]
  )
  
  cell_type_setting_combination$n_gene_included[i] <- length(unique(df_eqtl_cis$gene))
  cell_type_setting_combination$n_gene_with_cis_eqtl[i] <- n_unique_gene
  
  cat("N of genes with cis variants =",n_unique_gene,"\n")
  
  save_name <- paste0("eqtl_res_X_chr_",ct,"_",st,".rds")
  saveRDS(df_eqtl_cis,save_name)
  
  df_eqtl_cis_sort <- arrange(df_eqtl_cis,desc(significant_by_2step_FDR))
  df_eqtl_cis_sort <- select(df_eqtl_cis_sort,cell_type,setting,gene,significant_by_2step_FDR)
  colnames(df_eqtl_cis_sort)[4] <- "Has_cis_eQTL"
  df_eqtl_cis_sort <- distinct(df_eqtl_cis_sort,gene,cell_type,setting,.keep_all = T)
  
  df_eqtl_cis_sort_list[[i]] <- as.data.table(df_eqtl_cis_sort)
}
writexl::write_xlsx(cell_type_setting_combination,"N_of_cis_eQTL_X.xlsx")
df_eqtl_cis_gene_level <- rbindlist(
  df_eqtl_cis_sort_list,
  use.names = TRUE,
  fill = TRUE
)
writexl::write_xlsx(df_eqtl_cis_gene_level,"Details_of_cis_eQTL_X.xlsx")
# 3. LD clumping per setting
setting_all <- unique(cell_type_setting_combination$setting)
for (i in 1:length(setting_all)) {
  analysis_list <- filter(cell_type_setting_combination,setting==setting_all[i])
  
  df_eqtl_sig_list <- list()
  df_eqtl_sig_clumped_list <- list()
  
  for (j in 1:nrow(analysis_list)) {
    cat("Processing setting:",analysis_list$setting[j],"| cell type:",analysis_list$cell_type[j],"\n")
    read_name <- paste0("eqtl_res_X_chr_",analysis_list$cell_type[j],"_",analysis_list$setting[j],".rds")
    df <- readRDS(read_name)
    df_eqtl_sig <- filter(df,significant_by_2step_FDR=='Yes')
    df_eqtl_sig_list[[j]] <- df_eqtl_sig
    
    column_order <- colnames(df_eqtl_sig)
    if (nrow(df_eqtl_sig) == 0) {
      cat("  -> No significant eQTLs, skip clumping\n")
      df_eqtl_sig_clumped_list[[j]] <- df_eqtl_sig
      next
    }
    
    df_eqtl_sig_clumped <- df_eqtl_sig
    df_eqtl_sig_clumped$snp <- df_eqtl_sig_clumped$rsid
    df_eqtl_sig_clumped$rsid <- df_eqtl_sig_clumped$ID
    df_eqtl_sig_clumped$pval <- df_eqtl_sig_clumped$pvalue
    df_eqtl_sig_clumped$id <- df_eqtl_sig_clumped$gene
    
    df_eqtl_sig_clumped <- ld_clump(
      df_eqtl_sig_clumped,
      clump_kb  = 1000,
      clump_r2  = 0.1,
      clump_p   = 1,
      bfile     = bfile_prefix,
      plink_bin = plinkbinr::get_plink_exe()
    )
    
    df_eqtl_sig_clumped$rsid <- df_eqtl_sig_clumped$snp
    df_eqtl_sig_clumped$snp <- NULL
    df_eqtl_sig_clumped$pval <- NULL
    df_eqtl_sig_clumped$id <- NULL
    df_eqtl_sig_clumped <- select(df_eqtl_sig_clumped,column_order)
    
    df_eqtl_sig_clumped_list[[j]] <- df_eqtl_sig_clumped
  }
  df_eqtl_sig_all <- bind_rows(df_eqtl_sig_list)
  df_eqtl_sig_clumped_all <- bind_rows(df_eqtl_sig_clumped_list)
  saveRDS(df_eqtl_sig_all, paste0("eqtl_sig_all_X_chr_2stepFDR_before_prune_",setting_all[i],".rds"))
  saveRDS(df_eqtl_sig_clumped_all, paste0("eqtl_sig_clumped_all_X_chr_2stepFDR_after_prune_",setting_all[i],".rds"))
  writexl::write_xlsx(df_eqtl_sig_clumped_all, paste0("eqtl_sig_clumped_all_X_chr_2stepFDR_after_prune_",setting_all[i],".xlsx"))
}