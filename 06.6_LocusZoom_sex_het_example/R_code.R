# LocusZoom (locuszoomr) female-vs-male eQTL panels for sex-heterogeneous lead variants (P_het FDR < 0.1), +/-150 kb around the index variant
rm(list = ls())
library(tidyverse)
library(locuszoomr)
library(data.table)
library(EnsDb.Hsapiens.v86)   # BiocManager::install("EnsDb.Hsapiens.v86") if missing
# CONFIG
pl <- "plink1.9"
bfile_dir  <- "../00.1_X_chr_bfile_QC/"
bfile_name <- "chrX_QC_EU.0.G"
eqtl_rds   <- "../03_merge_eQTL_res_X/df_eqtl_cis_all.rds"
# independent lead variants with sex-het test (step 05)
indep_x_rds <- "../05_independent_signal_X/df_eqtl_wide_after_second_clump.rds"
sex_het_fdr <- 0.1        # lead variants plotted: P_het_FDR < sex_het_fdr
flank_bp    <- 150000     # half-window around the index variant
out_dir     <- "LCzoom_sex_het_FDR_lt0.1"
ens_db      <- "EnsDb.Hsapiens.v86"
highlight_col  <- "red"   # target gene colour in the gene track
female_setting <- "female_xchr_model_2"
male_setting   <- "male_xchr_model_1"
tmp_dir <- "temp"         # PLINK LD files
# page / panel layout
pdf_width  <- 3
pdf_height <- 4.875
panel_h    <- c(1.8, 1.8, 0.6)   # female scatter, male scatter, gene track
cex_gene   <- 0.65
mar_scatter <- c(0.6, 3.5, 2.4, 3.5)   # left/right must match across panels
mar_genes   <- c(3.2, 3.5, 0.3, 3.5)
headroom    <- 0.1                     # y-axis headroom above the highest -log10(P)
# 1. cis-eQTL results (female / male settings only)
df_eqtl_all <- readRDS(eqtl_rds)
df_eqtl_all <- dplyr::filter(df_eqtl_all, setting %in% c(female_setting, male_setting))
# 2. Lead variants to plot
cell_order <- c("Exc", "Inh", "Ast", "Oli", "OPC", "Mic", "End")
analysis_list <- readRDS(indep_x_rds) %>%
  dplyr::filter(!is.na(P_het_FDR), P_het_FDR < sex_het_fdr) %>%
  dplyr::select(gene, cell_type, ID, rsid, POS, EAF_female, EAF_male, P_het, P_het_FDR) %>%
  dplyr::mutate(cell_type = factor(cell_type, levels = cell_order)) %>%
  dplyr::arrange(cell_type, P_het) %>%
  dplyr::mutate(cell_type = as.character(cell_type))
print(analysis_list[, c("cell_type", "gene", "rsid", "P_het", "P_het_FDR")])
dir.create(tmp_dir, showWarnings = FALSE)
dir.create(out_dir, showWarnings = FALSE)
# 3. Helpers
# locuszoomr column format
prep_eqtl <- function(df) {
  df %>%
    dplyr::transmute(
      chrom         = "X",
      pos           = POS,
      rsid          = ID,
      other_allele  = REF,
      effect_allele = ALT,
      p             = pvalue,
      beta          = BETA,
      se            = SE
    ) %>%
    dplyr::filter(!is.na(p)) %>%
    dplyr::arrange(p)
}
# LD (r2) with the index variant via PLINK
get_ld <- function(snps, index_snp, tag) {
  snp_file <- file.path(tmp_dir, paste0("snp_list_", tag, ".txt"))
  out_pref <- file.path(tmp_dir, paste0("ld_res_", tag))
  write.table(unique(snps), snp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  command <- paste(
    pl,
    "--bfile", paste0(bfile_dir, bfile_name),
    "--allow-no-sex",
    "--r2",
    "--ld-snp", index_snp,
    "--ld-window 99999",
    "--ld-window-kb 99999",
    "--ld-window-r2 0",
    "--extract", snp_file,
    "--keep-allele-order",
    "--out", out_pref
  )
  system(command, ignore.stdout = TRUE)
  ld_file <- paste0(out_pref, ".ld")
  if (!file.exists(ld_file)) {
    warning("PLINK LD failed for ", tag, "; plotting without LD colours")
    return(data.frame(rsid = index_snp, R2 = 1))
  }
  ld_df <- fread(ld_file)
  dplyr::select(ld_df, rsid = SNP_B, R2)
}
# locus object; protein-coding genes plus the target gene
make_locus <- function(df, index_snp, xr, target_gene) {
  loc <- locus(
    df,
    gene      = NULL,
    seqname   = "X",
    xrange    = xr,
    index_snp = index_snp,
    LD        = "R2",
    ens_db    = ens_db
  )
  # recombination track needs internet access (UCSC); skipped if unavailable
  loc <- tryCatch(link_recomb(loc, genome = "hg38"), error = function(e) loc)
  loc$TX <- dplyr::filter(loc$TX, gene_biotype == "protein_coding" | gene_name == target_gene)
  loc
}
fmt_p <- function(x) formatC(x, format = "e", digits = 1)
add_panel_title <- function(main, sub) {
  mtext(main, side = 3, line = 1.25, font = 2, cex = 0.7)
  mtext(sub,  side = 3, line = 0.3,  font = 1, cex = 0.6)
}
# 4. One PDF per lead variant: female panel, male panel, gene track
for (i in seq_len(nrow(analysis_list))) {
  g       <- analysis_list$gene[i]
  ct      <- analysis_list$cell_type[i]
  top_snp <- analysis_list$ID[i]
  tag     <- paste0(g, "_", ct)
  cat("\n[", i, "/", nrow(analysis_list), "] ", tag, " index ", top_snp, "\n", sep = "")

  xrange_min <- analysis_list$POS[i] - flank_bp
  xrange_max <- analysis_list$POS[i] + flank_bp

  eqtl_sub <- dplyr::filter(df_eqtl_all, gene == g, cell_type == ct,
                            POS >= xrange_min, POS <= xrange_max)
  eqtl_female <- prep_eqtl(dplyr::filter(eqtl_sub, setting == female_setting))
  eqtl_male   <- prep_eqtl(dplyr::filter(eqtl_sub, setting == male_setting))

  # variants tested in both sexes
  eqtl_female <- dplyr::semi_join(eqtl_female, eqtl_male, by = "rsid")
  eqtl_male   <- dplyr::semi_join(eqtl_male, eqtl_female, by = "rsid")

  if (!top_snp %in% eqtl_female$rsid) {
    warning("Index variant ", top_snp, " not tested in both sexes for ", tag, "; skipped")
    next
  }

  ld_df <- get_ld(eqtl_female$rsid, top_snp, tag)
  eqtl_female <- dplyr::left_join(eqtl_female, ld_df, by = "rsid")
  eqtl_male   <- dplyr::left_join(eqtl_male,   ld_df, by = "rsid")

  loc_female <- make_locus(eqtl_female, top_snp, c(xrange_min, xrange_max), g)
  loc_male   <- make_locus(eqtl_male,   top_snp, c(xrange_min, xrange_max), g)

  # shared y axis across the two panels
  ymax <- max(-log10(c(eqtl_female$p, eqtl_male$p)), na.rm = TRUE)
  ylim_eqtl <- c(0, ymax + headroom)

  sub_title <- paste0(analysis_list$rsid[i],
                      "  |  P-het ", fmt_p(analysis_list$P_het[i]),
                      ", FDR ", signif(analysis_list$P_het_FDR[i], 2))
  sub_f <- paste0(sub_title, "  |  EAF ", round(analysis_list$EAF_female[i], 3))
  sub_m <- paste0(sub_title, "  |  EAF ", round(analysis_list$EAF_male[i], 3))

  pdf(file.path(out_dir, paste0("LCzoom_", g, "_", ct, "_female_male.pdf")),
      width = pdf_width, height = pdf_height)
  layout(matrix(1:3, ncol = 1), heights = panel_h)
  op <- par(mgp = c(2, 0.6, 0), tcl = -0.3)

  par(mar = mar_scatter)
  scatter_plot(loc_female, ylim = ylim_eqtl, labels = NULL,
               legend_pos = "topleft", xticks = FALSE, align = FALSE,
               cex = 0.8, cex.axis = 0.8, cex.lab = 0.9)
  add_panel_title(paste0("Female eQTL: ", g, " (", ct, ")"), sub_f)

  par(mar = mar_scatter)
  scatter_plot(loc_male, ylim = ylim_eqtl, labels = NULL,
               legend_pos = "topleft", xticks = FALSE, align = FALSE,
               cex = 0.8, cex.axis = 0.8, cex.lab = 0.9)
  add_panel_title(paste0("Male eQTL: ", g, " (", ct, ")"), sub_m)

  par(mar = mar_genes)
  genetracks(loc_female, highlight = g, highlight_col = highlight_col,
             align = FALSE, cex.axis = 0.8, cex.lab = 0.9, cex.text = cex_gene)

  par(op)
  dev.off()
}
