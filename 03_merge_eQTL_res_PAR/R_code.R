# Merge PAR eQTL results: gene TSS from GTF, cis window (TSS +/- 1 Mb), gene-wise cis-FDR
rm(list = ls())
library(data.table)
library(dplyr)
# CONFIG
eqtl_main_dir  <- "../02_eQTL_main/"
male_par_dir   <- paste0(eqtl_main_dir, "2_eQTL_PAR_chr_Male/")
female_par_dir <- paste0(eqtl_main_dir, "2_eQTL_PAR_chr_Female/")
both_par_dir   <- paste0(eqtl_main_dir, "2_eQTL_PAR_chr_Both/")
gtf_path <- "path/to/gencode.basic.annotation.gtf.gz"
cis_window <- 1000000
# PAR coordinates (hg38)
PAR1_start <- 10001
PAR1_end   <- 2781479
PAR2_start <- 155701383
PAR2_end   <- 156030895
make_glm_df <- function(base_dir) {
  cell_types <- c("Ast","End","Exc","Inh","Mic","Oli","OPC")
  out <- data.frame(path=character(), cell_type=character(), gene=character(),
                    stringsAsFactors = FALSE)
  
  for (ct in cell_types) {
    files <- list.files(file.path(base_dir, ct),
                        pattern = "glm\\.linear$",
                        full.names = TRUE)
    if (!length(files)) next
    genes <- sapply(strsplit(basename(files), "\\."), `[`, 2)
    out <- rbind(out, data.frame(path=files, cell_type=ct, gene=genes,
                                 stringsAsFactors = FALSE))
  }
  out
}
# 1. Collect result paths
df_male_par <- make_glm_df(male_par_dir)
df_male_par$setting <- "male_par"
df_female_par <- make_glm_df(female_par_dir)
df_female_par$setting <- "female_par"
df_both_par <- make_glm_df(both_par_dir)
df_both_par$setting <- "both_par"
df_all <- rbind(df_male_par,
                df_female_par,
                df_both_par)
dim(df_all)
# 2. Gene-level TSS from GTF
gtf <- fread(gtf_path)
if (ncol(gtf) < 9) stop("GTF has fewer than 9 columns")
setnames(
  gtf,
  names(gtf)[1:9],
  c("seqname","source","feature","start","end","score","strand","frame","attribute")
)
gene_anno <- gtf[feature == "gene"]
gene_anno[, gene_id := sub('.*gene_id "([^"]+)".*', "\\1", attribute)]
gene_anno[, gene_symbol := sub('.*gene_name "([^"]+)".*', "\\1", attribute)]
gene_anno[, gene_type := sub('.*gene_type "([^"]+)".*', "\\1", attribute)]
# unmatched regex -> NA
gene_anno[gene_id == attribute, gene_id := NA_character_]
gene_anno[gene_symbol == attribute, gene_symbol := NA_character_]
gene_anno[gene_type == attribute, gene_type := NA_character_]
gene_anno$gene_id <- sub("\\.[0-9]+$", "", gene_anno$gene_id)
gene_map <- unique(gene_anno[, .(gene_id, gene_symbol, gene_type)])
tx <- gtf[feature == "transcript"]
tx[, gene_id := sub('.*gene_id "([^"]+)".*', "\\1", attribute)]
tx[gene_id == attribute, gene_id := NA_character_]
tx$gene_id <- sub("\\.[0-9]+$", "", tx$gene_id)
tx[, tss := fifelse(strand == "+", start, end)]
tss_tx <- tx[, .(
  chr    = as.character(seqname),
  tss    = as.integer(tss),
  strand = as.character(strand),
  gene_id = gene_id
)]
# gene TSS: '+' smallest, '-' largest transcript TSS
tss_gene_one_df <- tss_tx[
  !is.na(gene_id),
  {stopifnot(length(unique(strand)) == 1)
    if (strand[1] == "+") .SD[which.min(tss)] else .SD[which.max(tss)]},
  by = gene_id
]
tss_gene_one_df <- merge(
  tss_gene_one_df,
  gene_anno[, .(
    gene_id,
    gene_start = as.integer(start),
    gene_end   = as.integer(end)
  )],
  by = "gene_id",
  all.x = TRUE
)
tss_gene_one_df <- merge(
  tss_gene_one_df,
  gene_map,
  by = "gene_id",
  all.x = TRUE
)
saveRDS(tss_gene_one_df,"TSS_all_chr.rds")
# 3. Tested chrX genes: PAR membership and cis window
tss_gene_X_chr <- filter(tss_gene_one_df,tss_gene_one_df$chr=='chrX')
unique_eQTL_genes <- unique(df_all$gene)
tss_gene_X_chr <- filter(tss_gene_X_chr,tss_gene_X_chr$gene_symbol %in% unique_eQTL_genes)
nrow(tss_gene_X_chr)==length(unique_eQTL_genes)   # check
setnames(tss_gene_X_chr, old = "gene_symbol", new = "gene")
tss_gene_X_chr$tss_in_PAR <- with(
  tss_gene_X_chr,
  (tss >= PAR1_start & tss <= PAR1_end) |
    (tss >= PAR2_start & tss <= PAR2_end)
)
# PAR_type: whole gene body within PAR1/PAR2
tss_gene_X_chr$PAR_type <- with(
  tss_gene_X_chr,
  ifelse(
    gene_start >= PAR1_start & gene_end <= PAR1_end, "PAR1",
    ifelse(
      gene_start >= PAR2_start & gene_end <= PAR2_end, "PAR2",
      "non_PAR"
    )
  )
)
tss_gene_X_chr$cis_low <- pmax(tss_gene_X_chr$tss - cis_window, 1)
tss_gene_X_chr$cis_high <- tss_gene_X_chr$tss + cis_window
writexl::write_xlsx(tss_gene_X_chr,"tss_gene_X_chr.xlsx")
df_all <- left_join(df_all,tss_gene_X_chr,by='gene')
which(is.na(df_all$cis_low))
df_all <- filter(df_all,PAR_type %in% c("PAR1","PAR2"))
# 4. cis variants per gene; gene-wise cis-FDR
cis_list   <- list()
process_list <- list()
for (i in 1:nrow(df_all)) {
  df_eqtl <- fread(df_all$path[i])
  df_eqtl$pvalue <- 10^(-df_eqtl$LOG10_P)
  
  gene <- df_all$gene[i]
  gene_id <- df_all$gene_id[i]
  cell_type <- df_all$cell_type[i]
  cis_high <- df_all$cis_high[i]
  cis_low <- df_all$cis_low[i]
  par_type <- df_all$PAR_type[i]
  setting <- df_all$setting[i]
  
  df_eqtl_cis <- filter(df_eqtl, POS >= cis_low & POS <= cis_high)
  df_eqtl_cis$cell_type <- cell_type
  df_eqtl_cis$gene <- gene
  df_eqtl_cis$par_type <- par_type
  df_eqtl_cis$setting <- setting
  df_eqtl_cis <- filter(df_eqtl_cis,!is.na(df_eqtl_cis$LOG10_P))
  df_eqtl_cis$fdr <- p.adjust(df_eqtl_cis$pvalue,method = 'fdr')
  
  cat("process gene:",gene,"| cell type:",cell_type,"| setting:",setting,"|",nrow(df_eqtl_cis),"cis variants extracted \n")
  
  df_process <- df_all[i,]
  df_process$path <- NULL
  df_process$n_cis_variants_all <- nrow(df_eqtl_cis)
  
  cis_list[[i]]   <- df_eqtl_cis
  process_list[[i]] <- df_process
}

# 5. Combine and save
df_eqtl_cis_all   <- rbindlist(cis_list,   use.names = TRUE, fill = TRUE)
df_process_all <- rbindlist(process_list, use.names = TRUE, fill = TRUE)
saveRDS(df_eqtl_cis_all,   file = "df_eqtl_cis_all.rds")
writexl::write_xlsx(df_process_all,"df_process_all.xlsx")