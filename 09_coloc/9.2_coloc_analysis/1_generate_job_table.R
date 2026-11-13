# coloc job table: XWAS .ma x (cell type, setting)
library(data.table)
library(dplyr)
# CONFIG
eqtl_rds   <- "../../03_merge_eQTL_res_X/df_eqtl_cis_all.rds"
afreq_file <- "../../00.1_X_chr_bfile_QC/chrX_QC_EU.0.G.afreq"
gwas_dir   <- "../9.1_xwas_clean/"
# 1. cis-eQTL with EAF
df_eQTL_all      <- readRDS(eqtl_rds)
df_eQTL_distinct <- distinct(df_eQTL_all, cell_type, gene, setting)
cell_types       <- unique(df_eQTL_distinct$cell_type)
df_afreq <- fread(afreq_file)
df_afreq$OBS_CT <- NULL
df_eQTL_all <- left_join(df_eQTL_all, df_afreq)
# 2. XWAS list
GWAS_list <- dir(gwas_dir, pattern = ".ma")
GWAS_list <- data.frame(trait     = substr(GWAS_list, 1, nchar(GWAS_list) - 3),
                        GWAS_path = paste0(gwas_dir, GWAS_list))
GWAS_list$GWAS_type <- sub(".*_", "", GWAS_list$trait)
# 3. Job table
group_df <- df_eQTL_distinct %>%
  distinct(cell_type, setting)
job_table <- expand.grid(
  gwas_i  = seq_len(nrow(GWAS_list)),
  group_i = seq_len(nrow(group_df))
) %>%
  as_tibble() %>%
  left_join(group_df %>% mutate(group_i = row_number()),
            by = "group_i") %>%
  mutate(job_id = row_number())
fwrite(job_table, "job_table.txt", sep = "\t")