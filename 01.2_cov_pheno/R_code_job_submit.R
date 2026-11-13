# Per (sex x cell type) job: chrX phenotype table (TMM/voom, ComBat, inverse-normal) and covariate table (covariates + ePCs + gPCs); job index selects the row of analysis_ref_df
rm(list = ls())
library(tidyverse)
library(data.table)
library(edgeR)
library(sva)
library(limma)
library(HGNChelper)
library(plinkbinr)
# CONFIG
pseudobulk_dir <- "path/to/pseudobulk/"   # *.pseudobulk.rds and *.covariates.tsv per cell type
valid_id_dir   <- "../01.1_valid_sample_id/"
gencode_gtf    <- "path/to/gencode.basic.annotation.gtf.gz"
metadata_xlsx  <- "path/to/metadata.xlsx"   # specimenID / batch / msex
bfile_dir      <- "path/to/wgs_genotypes"   # autosomal genotypes for gPCs
# 0. Analysis reference table (one row per sex x cell type)
gene_exp_rds_path <- dir(pseudobulk_dir, pattern = "rds")
gene_exp_rds_path_df <- data.frame(
  cell_type     = sub("\\.pseudobulk\\.rds$", "", gene_exp_rds_path),
  exp_file_path = paste0(pseudobulk_dir, gene_exp_rds_path),
  stringsAsFactors = FALSE
)
cov_path <- dir(pseudobulk_dir, pattern = "covariates.tsv")
cov_path_df <- data.frame(
  cell_type     = sub(".covariates.tsv", "", cov_path),
  cov_file_path = paste0(pseudobulk_dir, cov_path),
  stringsAsFactors = FALSE
)
analysis_ref_df <- expand.grid(
  sex       = c("both", "female", "male"),
  cell_type = c("Ast", "End", "Exc", "Inh", "Mic", "Oli", "OPC"),
  stringsAsFactors = FALSE
)
analysis_ref_df$valid_id <- paste0(valid_id_dir, "valid_sample_id_",
                                   analysis_ref_df$sex, "_", analysis_ref_df$cell_type, ".txt")
analysis_ref_df$valid_id_bothsex_for_exp_norm <- paste0(valid_id_dir, "valid_sample_id_",
                                                        "both_", analysis_ref_df$cell_type, ".txt")
analysis_ref_df <- left_join(analysis_ref_df, gene_exp_rds_path_df)
analysis_ref_df <- left_join(analysis_ref_df, cov_path_df)
analysis_ref_df$exp_output_path <- paste0(analysis_ref_df$sex, "_",
                                          analysis_ref_df$cell_type, ".pheno_chrX.tsv")
analysis_ref_df$cov_output_path <- paste0(analysis_ref_df$sex, "_",
                                          analysis_ref_df$cell_type, ".cov.txt")
analysis_ref_df$bfile_for_gPC <- paste0("bfile_for_gPC/", analysis_ref_df$sex, "_", analysis_ref_df$cell_type)
analysis_ref_df$prune_path    <- paste0("prune/", analysis_ref_df$sex, "_", analysis_ref_df$cell_type)
analysis_ref_df$genetic_pc    <- paste0("genetic_PC/genetic_PC_", analysis_ref_df$sex, "_", analysis_ref_df$cell_type, ".txt")
analysis_ref_df$gPC_path      <- paste0("gPC/gPC_", analysis_ref_df$sex, "_", analysis_ref_df$cell_type)
dir.create("bfile_for_gPC", showWarnings = FALSE)
dir.create("genetic_PC",    showWarnings = FALSE)
dir.create("prune",         showWarnings = FALSE)
dir.create("gPC",           showWarnings = FALSE)
# 1. Gene-symbol exclusion list: non-approved symbols whose HGNC-suggested symbol collides with another gene
expr_list <- lapply(gene_exp_rds_path_df$exp_file_path, readRDS)
names(expr_list) <- basename(gene_exp_rds_path_df$exp_file_path)
ref <- rownames(expr_list[[1]])
identical_flag <- sapply(expr_list, function(x) identical(ref, rownames(x)))
print(identical_flag)
if (all(identical_flag)) {
  cat("\nAll rownames are IDENTICAL across pseudobulk matrices\n")
  cat("Generating gene symbol exclusion list \n")
  
  exp_df <- expr_list[[1]]
  exp_df_check_gene_symbol <- checkGeneSymbols(rownames(exp_df))
  
  exp_df_check_gene_symbol_num <- as.data.frame(table(exp_df_check_gene_symbol$Suggested.Symbol))
  exp_df_check_gene_symbol_num <- filter(exp_df_check_gene_symbol_num, Freq > 1)
  
  exp_df_check_gene_symbol_duplicate <- filter(exp_df_check_gene_symbol,
                                               Suggested.Symbol %in% exp_df_check_gene_symbol_num$Var1)
  exp_df_check_gene_symbol_duplicate_remove <- filter(exp_df_check_gene_symbol_duplicate, Approved == FALSE)
  genes_to_remove <- exp_df_check_gene_symbol_duplicate_remove$x
  
} else {
  stop(
    "\nERROR: rownames are NOT identical across pseudobulk matrices.\n",
    "Non-identical matrices:\n  ",
    paste(names(identical_flag)[!identical_flag], collapse = "\n  "),
    call. = FALSE
  )
}
# 2. chrX gene set (GENCODE)
gencode_df <- fread(gencode_gtf)
# GTF columns: V1 = seqname, V3 = feature, V9 = attributes
gencode_df_x_chr <- gencode_df[gencode_df$V1 == "chrX", ]
gencode_df_x_chr <- gencode_df_x_chr[gencode_df_x_chr$V3 == "gene", ]
gencode_df_x_chr$gene_id   <- sub('.*gene_id "([^"]+)".*', "\\1", gencode_df_x_chr$V9)
gencode_df_x_chr$gene_id   <- sub("\\..*$", "", gencode_df_x_chr$gene_id)
gencode_df_x_chr$gene_name <- sub('.*gene_name "([^"]+)".*', "\\1", gencode_df_x_chr$V9)
gencode_df_x_chr$gene_type <- sub('.*gene_type "([^"]+)".*', "\\1", gencode_df_x_chr$V9)
gencode_df_x_chr_hgnc_check <- checkGeneSymbols(gencode_df_x_chr$gene_name)
gencode_df_x_chr_hgnc_check <- filter(gencode_df_x_chr_hgnc_check, gencode_df_x_chr_hgnc_check$Approved == TRUE)
# 3. Per-job processing
args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])
cat("Running job index:", i, "\n")
# 3.1 Expression normalization (on the both-sex sample set)
cat("Normalization of gene expression matrix from the pseudobulk matrix for:",
    analysis_ref_df$sex[i], analysis_ref_df$cell_type[i], "\n")
valid_id          <- fread(analysis_ref_df$valid_id[i])
valid_id_for_norm <- fread(analysis_ref_df$valid_id_bothsex_for_exp_norm[i])
expression_pseudobulk_df <- readRDS(analysis_ref_df$exp_file_path[i])
expression_pseudobulk_df <- as.data.frame(expression_pseudobulk_df)
expression_pseudobulk_df <- select(expression_pseudobulk_df, all_of(valid_id_for_norm$IID))
dge  <- DGEList(counts = expression_pseudobulk_df)
keep <- filterByExpr(dge)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
v <- voom(dge, plot = TRUE)
logcpm <- v$E
mean_logcpm <- apply(logcpm, 1, mean)
logcpm <- logcpm[mean_logcpm > 2.0, ]
# ComBat batch correction, preserving sex
metadata <- readxl::read_excel(metadata_xlsx)
metadata <- select(metadata, specimenID, batch, msex)
metadata$msex <- factor(metadata$msex)
metadata <- metadata[match(colnames(logcpm), metadata$specimenID), ]
stopifnot(all(metadata$specimenID == colnames(logcpm)))
batch_n <- table(metadata$batch)
print(batch_n)
if (any(batch_n < 2)) {
  warning("Some batches have <2 samples; ComBat can be unstable.")
}
batches <- metadata$batch
mod     <- model.matrix(~ msex, data = metadata)
logcpm  <- ComBat(dat = logcpm, batch = batches, mod = mod)
logcpm_for_PC <- logcpm
# rank-based inverse normal transformation
logcpm <- t(apply(logcpm, 1, rank, ties.method = "average"))
logcpm <- qnorm(logcpm / (ncol(logcpm) + 1))
# chrX phenotype table
pheno_chrX <- as.data.frame(t(logcpm))
n_remove <- sum(colnames(pheno_chrX) %in% genes_to_remove)
cat("Removing", n_remove,
    "duplicate/alias genes before gene symbol updates in all chromosomes in",
    analysis_ref_df$sex[i], analysis_ref_df$cell_type[i], "\n")
pheno_chrX <- pheno_chrX[, !colnames(pheno_chrX) %in% genes_to_remove, drop = FALSE]
pheno_chrX_hgnc_gene_check <- checkGeneSymbols(colnames(pheno_chrX))
pheno_chrX_hgnc_gene_check <- filter(pheno_chrX_hgnc_gene_check,
                                     Suggested.Symbol %in% gencode_df_x_chr_hgnc_check$Suggested.Symbol)
pheno_chrX <- select(pheno_chrX, pheno_chrX_hgnc_gene_check$x)
n_alias <- sum(pheno_chrX_hgnc_gene_check$Approved == FALSE)
cat("Updated", n_alias, "non-approved (alias) gene symbols from",
    analysis_ref_df$sex[i], analysis_ref_df$cell_type[i], "\n")
colnames(pheno_chrX) <- pheno_chrX_hgnc_gene_check$Suggested.Symbol
pheno_chrX$FID <- rownames(pheno_chrX)
pheno_chrX$IID <- pheno_chrX$FID
pheno_chrX <- select(pheno_chrX, FID, IID, everything())
pheno_chrX <- filter(pheno_chrX, pheno_chrX$IID %in% valid_id$IID)
pheno_chrX <- pheno_chrX[match(valid_id$IID, pheno_chrX$IID), ]
fwrite(pheno_chrX, analysis_ref_df$exp_output_path[i], sep = "\t", quote = FALSE, row.names = FALSE)
# 3.2 Covariate table
cat("Formatting covariate table for:",
    analysis_ref_df$sex[i], analysis_ref_df$cell_type[i], "\n")
cov_df <- fread(analysis_ref_df$cov_file_path[i])
cov_df <- column_to_rownames(cov_df, var = "id")
cov_df <- as.data.frame(t(cov_df))
cov_df <- cov_df[valid_id$IID, ]
cov_df$FID <- rownames(cov_df)
cov_df$IID <- cov_df$FID
cov_df <- select(cov_df, FID, IID, everything())
# PMI: impute NA with the median
cov_df$pmi <- as.numeric(cov_df$pmi)
pmi_median <- median(cov_df$pmi, na.rm = TRUE)
cov_df$pmi[is.na(cov_df$pmi)] <- pmi_median
# ePCs: top 1000 variable non-chrX genes, first 30 PCs
df_exp_PC <- select(as.data.frame(logcpm_for_PC), all_of(valid_id_for_norm$IID))
df_exp_PC <- df_exp_PC[!rownames(df_exp_PC) %in% gencode_df_x_chr$gene_name, , drop = FALSE]
gene_sd   <- apply(df_exp_PC, 1, sd, na.rm = TRUE)
top_genes <- names(sort(gene_sd, decreasing = TRUE))[1:1000]
df_exp_PC <- df_exp_PC[top_genes, ]
pca     <- prcomp(t(df_exp_PC), scale. = TRUE)
expr_pc <- as.data.frame(pca$x[, 1:30])
colnames(expr_pc) <- gsub("PC", "ePC", colnames(expr_pc))
expr_pc$IID <- rownames(expr_pc)
cov_df <- cov_df %>% select(-starts_with("ePC"))
cov_df <- left_join(cov_df, expr_pc, by = "IID")
# gPCs: autosomes, LD-pruned, 3 PCs
cmd <- paste0(
  get_plink_exe(), " ",
  "--bfile ", bfile_dir, " ",
  "--chr 1-22 ",
  "--keep ", analysis_ref_df$valid_id_bothsex_for_exp_norm[i], " ",
  "--snps-only just-acgt ",
  "--make-bed ",
  "--out ", analysis_ref_df$bfile_for_gPC[i]
)
cat(cmd, "\n")
system(cmd)
cmd <- paste0(
  get_plink_exe(), " ",
  "--bfile ", analysis_ref_df$bfile_for_gPC[i], " ",
  "--keep ", analysis_ref_df$valid_id_bothsex_for_exp_norm[i], " ",
  "--geno 0.001 ",
  "--maf 0.01 ",
  "--indep-pairwise 1500 150 0.01 ",
  "--out ", analysis_ref_df$prune_path[i]
)
cat(cmd, "\n")
system(cmd)
cmd <- paste0(
  get_plink_exe(), " ",
  "--bfile ", analysis_ref_df$bfile_for_gPC[i], " ",
  "--keep ", analysis_ref_df$valid_id_bothsex_for_exp_norm[i], " ",
  "--extract ", analysis_ref_df$prune_path[i], ".prune.in ",
  "--pca 3 ",
  "--out ", analysis_ref_df$gPC_path[i]
)
cat(cmd, "\n")
system(cmd)
gPC_res <- fread(paste0(analysis_ref_df$gPC_path[i], ".eigenvec"), header = FALSE)
colnames(gPC_res) <- c("FID", "IID", paste0("gPC", 1:(ncol(gPC_res) - 2)))
cov_df <- cov_df %>% select(-starts_with("gPC"))
cov_df <- left_join(cov_df, gPC_res, by = c("FID", "IID"))
fwrite(cov_df, analysis_ref_df$cov_output_path[i], sep = "\t", quote = FALSE, row.names = FALSE)