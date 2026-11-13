# cis-eQTL results -> SMR input: flist, one .esd per probe (gene_cell type_setting), dense BESD
library(data.table)
library(dplyr)
# CONFIG
eqtl_rds  <- "../../03_merge_eQTL_res_X/df_eqtl_cis_all.rds"
afreq_txt <- "../../00.1_X_chr_bfile_QC/chrX_QC_EU.0.G.afreq"
tss_rds   <- "../../03_merge_eQTL_res_X/TSS_all_chr.rds"
esd_dir   <- "esd_files"
flist_out <- "my.flist"
besd_out  <- "eQTL_besd"
smr_bin   <- "smr"   # SMR binary (v1.3.1) placed in this folder; download from https://yanglab.westlake.edu.cn/software/smr/
# 1. Load cis-eQTL; male model 1 (0/2) rescaled to the model 2 (0/1) scale
df_merged <- readRDS(eqtl_rds)
df_merged$gene_id <- paste0(df_merged$gene, "_", df_merged$cell_type, "_", df_merged$setting)
df_merged <- df_merged %>%
  mutate(
    BETA    = ifelse(setting == "male_xchr_model_1", BETA / 2, BETA),
    SE      = ifelse(setting == "male_xchr_model_1", SE   / 2, SE),
    setting = ifelse(setting == "male_xchr_model_1", "male_xchr_model_2", setting)
  )
# 2. Annotate EAF, p-value, gene start
df_eaf <- fread(afreq_txt)
df_eaf$OBS_CT <- NULL
df_merged <- left_join(df_merged, df_eaf, by = c("#CHROM", "ID", "REF", "ALT"))
df_merged$pvalue <- 10^(-df_merged$LOG10_P)
tss_df <- readRDS(tss_rds)
tss_df <- dplyr::select(tss_df, gene = gene_symbol, gene_start)
tss_df <- dplyr::distinct(tss_df, gene, .keep_all = TRUE)
df_merged <- left_join(df_merged, tss_df, by = "gene")
# 3. flist
my_flist <- distinct(df_merged, gene_id, .keep_all = TRUE)
my_flist <- my_flist[, c("#CHROM", "gene_id", "gene_start", "gene")]
my_flist$GeneticDistance <- 0
my_flist$Orientation     <- NA
my_flist$PathOfEsd       <- paste0(esd_dir, "/", my_flist$gene_id, ".esd")
my_flist <- dplyr::select(my_flist,
                          "#CHROM", gene_id, GeneticDistance, "gene_start",
                          "gene", Orientation, PathOfEsd)
colnames(my_flist) <- c("Chr", "ProbeID", "GeneticDistance", "ProbeBp",
                        "Gene", "Orientation", "PathOfEsd")
write.table(my_flist, flist_out,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE,
            fileEncoding = "UTF-8", eol = "\n")
# 4. .esd per probe
dir.create(esd_dir, showWarnings = FALSE)
df_dt    <- as.data.table(df_merged)
df_split <- split(df_dt, df_dt$gene_id)
n_total  <- length(df_split)
for (i in seq_along(df_split)) {
  probe_ID <- names(df_split)[i]
  esd <- df_split[[i]][, .(`#CHROM`, ID, POS, ALT, REF, ALT_FREQS, BETA, SE, pvalue)]
  setnames(
    esd,
    old = c("#CHROM", "ID", "POS", "ALT", "REF", "ALT_FREQS", "BETA", "SE", "pvalue"),
    new = c("Chr", "SNP", "Bp", "A1", "A2", "Freq", "Beta", "se", "p")
  )
  fwrite(esd, file = file.path(esd_dir, paste0(probe_ID, ".esd")), sep = "\t", quote = FALSE)
  
  if (i %% 100 == 0 || i == n_total) {
    cat(sprintf("[%s] %d/%d (%.1f%%)\n",
                format(Sys.time(), "%H:%M:%S"), i, n_total, 100 * i / n_total))
  }
}
# 5. Dense BESD
system(paste("chmod +x", smr_bin))
system(paste0("./", smr_bin, " --eqtl-flist ", flist_out, " --make-besd-dense --out ", besd_out))