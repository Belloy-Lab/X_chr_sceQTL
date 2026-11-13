# Merge TWAS + SMR + COLOC; support flags; stacked AD/PD/MSA heatmap (Fig 5B) and LocusZoom plots of the prioritized genes
rm(list = ls())
library(data.table)
library(tidyverse)          # dplyr / tidyr / tibble / stringr / readr / ggplot2
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(locuszoomr)
library(EnsDb.Hsapiens.v86)   # masks dplyr::select / filter, so dplyr:: is used explicitly in section 5
# ==============================================================================
# CONFIG
# ==============================================================================
# --- input result tables ---
twas_rds  <- "../07_twas/7.8_combine_twas_results/twas_res_distinct.rds"
smr_rds   <- "../08_smr/8.4_combine_results/smr_res_distinct.rds"
coloc_rds <- "../09_coloc/9.3_combine_results/coloc_res_distinct.rds"
deg_xlsx  <- "../06.3_sex_DEG/DEG_female_vs_male_X_chr_gene_non_adj_disease.xlsx"
eqtl_rds  <- "../03_merge_eQTL_res_X/df_eqtl_cis_all.rds"
# --- thresholds ---
fdr_smr_cut    <- 0.05    # SMR support: FDR_SMR < this
heidi_cut      <- 0.05    # SMR support: p_HEIDI NA or > this
coloc_pp4_cut  <- 0.6     # COLOC support: PP4 > this
pp4_label_cut  <- 0.5     # in-cell PP4 label shown when PP4 > this
deg_p_cut      <- 0.05    # DEG sex-bias nominal p
# --- combined stacked figure ---
target_traits  <- c("AD", "PD", "MSA")
draw_black_box       <- FALSE   # black boxes are added manually on the final figure, not drawn here
dagger_symbol        <- "\u2020"  # dagger: cTWAS-significant cell whose sex-DEG is concordant with the locus XWAS model
# --- right-side annotation rule (per gene), among cells with SMR or COLOC support:
#     (1) one XWAS model per locus = model of the strongest cTWAS cell in the locus;
#     (2) prefer the cell whose sex-DEG result is concordant with that model
#         (eXCI/female -> female-up; male -> male-up; rXCI -> no difference),
#         ties broken by the strongest cTWAS. DEG / SMR / COLOC columns are read at that cell ---
# genes grouped into their own row slice (locus grouping) per trait, in the stacked figure;
# also defines the loci used for the "one XWAS model per locus" annotation
box_groups <- list(
  AD  = list(c("KRBOX4", "CHST7", "SLC9A7", "JADE3")),
  PD  = list(c("RAB40A", "MORF4L2")),
  MSA = list(c("SLC6A8", "PLXNA3", "FAM3A", "IKBKG"))
)
# --- cell types ---
cell_types <- c("Ast", "End", "Exc", "Inh", "Mic", "OPC", "Oli")
# --- gwas models + palette ---
models <- c("female", "male", "eXCI", "rXCI")
model_colors <- c(
  female = "#e31a1c", male = "#1f78b4", eXCI = "#ff7f00", rXCI = "#33a02c"
)
na_color <- "grey85"
# --- TWAS Z color scale: negative = deep purple, 0 = white, positive = deep gold ---
twas_low_col  <- "#4B0082"
twas_high_col <- "#B8860B"
# --- DEG sex-bias levels + palette ---
deg_levels <- c("Upregulated in female", "Upregulated in male", "No significant difference")
deg_colors <- c(`Upregulated in female`     = "#e31a1c",
                `Upregulated in male`       = "#1f78b4",
                `No significant difference` = na_color)
support_colors <- c(yes = "#333333", no = na_color)   # yes = dark, no = grey
# ==============================================================================
# 1. Merge TWAS + SMR + COLOC
# ==============================================================================
twas_res_distinct  <- readRDS(twas_rds)
smr_res_distinct   <- readRDS(smr_rds)
coloc_res_distinct <- readRDS(coloc_rds)
twas_smr_coloc_merged <- twas_res_distinct %>%
  left_join(smr_res_distinct,
            by = c("Gene", "gwas_model", "cell_type", "trait")) %>%
  left_join(coloc_res_distinct,
            by = c("Gene", "gwas_model", "cell_type", "trait"))
# ==============================================================================
# 2. FDR + support flags (per trait + cell_type + gwas_model)
# ==============================================================================
twas_smr_coloc_merged <- twas_smr_coloc_merged %>%
  group_by(trait, cell_type, gwas_model) %>%
  mutate(
    FDR_SMR       = p.adjust(p_SMR, method = "fdr"),
    TWAS_SUPPORT  = !is.na(eqtl_support),
    SMR_SUPPORT   = FDR_SMR < fdr_smr_cut & (is.na(p_HEIDI) | p_HEIDI > heidi_cut),
    COLOC_SUPPORT = PP4 > coloc_pp4_cut
  ) %>%
  ungroup()
# order by CHR, BP
twas_smr_coloc_merged <- arrange(twas_smr_coloc_merged, CHR, gene_start)
# final support: TWAS + (SMR OR COLOC)
twas_smr_coloc_merged$final_support <-
  twas_smr_coloc_merged$TWAS_SUPPORT &
  (twas_smr_coloc_merged$SMR_SUPPORT | twas_smr_coloc_merged$COLOC_SUPPORT)
traits <- unique(twas_smr_coloc_merged$trait)
# ==============================================================================
# 3. Helper: pad matrices to all 7 cell types (missing columns -> NA / grey)
# ==============================================================================
# pad to all 7 cell types
force_all_celltypes <- function(m, all_ct) {
  missing_ct <- setdiff(all_ct, colnames(m))
  if (length(missing_ct) > 0) {
    add_mat <- matrix(NA_real_, nrow = nrow(m), ncol = length(missing_ct),
                      dimnames = list(rownames(m), missing_ct))
    m <- cbind(m, add_mat)
  }
  m[, all_ct, drop = FALSE]
}
# ==============================================================================
# 4. Sex-DEG lookup per gene x cell type (nominal P < deg_p_cut; logFC > 0 = higher in female)
# ==============================================================================
df_deg_extra <- readxl::read_excel(deg_xlsx) %>%
  mutate(
    DEG_sex_bias = case_when(
      logFC > 0 & P.Value < deg_p_cut ~ "Upregulated in female",
      logFC < 0 & P.Value < deg_p_cut ~ "Upregulated in male",
      TRUE ~ "No significant difference"
    )
  ) %>%
  select(cell_type, Gene = updated_Gene_symbol, DEG_sex_bias) %>%
  distinct(cell_type, Gene, .keep_all = TRUE)
# ==============================================================================
# 5. Fig 5B: AD / PD / MSA stacked heatmap (final_support genes only; shared cell-type columns
#    and legends; right columns: locus XWAS model / DEG sex bias / SMR / COLOC)
# ==============================================================================
target_traits_present <- target_traits[target_traits %in% traits]  # keep only those present in data
# ---- cell_fun factory: capture this trait's matrix and df ----
make_cf <- function(mat_ref, df_ref) {
  function(j, i, x, y, width, height, fill) {
    gene <- rownames(mat_ref)[i]; xax <- colnames(mat_ref)[j]
    key  <- paste(gene, xax, sep = "__")
    idx  <- which(df_ref$key == key)
    if (length(idx) == 0) return()
    row <- df_ref[idx[1], ]
    # black boxes (off by default; added manually on the final figure)
    if (draw_black_box && isTRUE(row$box_cell)) {
      grid.rect(x = x, y = y, width = width * 0.95, height = height * 0.95,
                gp = gpar(col = "black", lwd = 1.5, fill = NA))
    }
    # gwas_model label (rxci/exci/...) intentionally not drawn in the stacked figure
    # SMR (#) and COLOC (PP4) support are only shown for cells passing cTWAS FDR (*)
    twas_sig  <- isTRUE(row$TWAS_SUPPORT)
    # dagger: sex-DEG (nominal P < deg_p_cut) in this cell concordant with the locus XWAS model
    label_top <- paste0(ifelse(twas_sig, "*", ""),
                        ifelse(twas_sig && isTRUE(row$SMR_SUPPORT), "#", ""),
                        ifelse(twas_sig && isTRUE(row$deg_dagger), dagger_symbol, ""))
    if (label_top != "") grid.text(label_top, x, y + unit(1, "pt"),
                                   gp = gpar(fontsize = 7, col = "black"))
    if (twas_sig && isTRUE(row$COLOC_SUPPORT) && !is.na(row$PP4) && row$PP4 > pp4_label_cut) {
      grid.text(sprintf("%.2f", row$PP4), x, y - unit(6, "pt"),
                gp = gpar(fontsize = 5.5, col = "black"))
    }
  }
}
# ---- build per-trait: raw matrix / df / right-annotation source ----
build_trait_block <- function(trait_name) {
  d <- filter(twas_smr_coloc_merged, trait == trait_name)
  sig_gene <- unique(filter(d, TWAS_SUPPORT == TRUE)$Gene)
  d_for_plot <- filter(d, Gene %in% sig_gene) %>%
    group_by(cell_type, Gene) %>%
    slice_max(abs(TWAS.Z_best), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(CHR, gene_start)
  if (nrow(d_for_plot) == 0) return(NULL)
  
  df_local <- d_for_plot %>%
    mutate(x_axis = cell_type) %>%
    select(Gene, x_axis, gwas_model, TWAS.Z_best, TWAS_SUPPORT,
           SMR_SUPPORT, COLOC_SUPPORT, PP4) %>%
    distinct() %>%
    complete(Gene, x_axis,
             fill = list(TWAS.Z_best = NA_real_, TWAS_SUPPORT = FALSE,
                         SMR_SUPPORT = FALSE, COLOC_SUPPORT = FALSE,
                         PP4 = NA_real_, gwas_model = ""))
  df_local$key <- paste(df_local$Gene, df_local$x_axis, sep = "__")
  
  m <- df_local %>%
    select(Gene, x_axis, TWAS.Z_best) %>%
    pivot_wider(names_from = x_axis, values_from = TWAS.Z_best) %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  m <- force_all_celltypes(m, cell_types)
  
  gene_order <- unique(d_for_plot$Gene)
  gene_order <- gene_order[gene_order %in% rownames(m)]
  m <- m[gene_order, , drop = FALSE]
  
  # keep only final_support genes
  sig_final <- unique(filter(d_for_plot, final_support == TRUE)$Gene)
  m <- m[rownames(m) %in% sig_final, , drop = FALSE]
  if (nrow(m) == 0) return(NULL)
  m <- force_all_celltypes(m, cell_types)
  
  # ---- supported cells (final_support) with their sex-DEG label ----
  sup <- d_for_plot %>%
    filter(Gene %in% rownames(m), final_support == TRUE) %>%
    left_join(df_deg_extra, by = c("Gene", "cell_type")) %>%
    mutate(DEG_sex_bias = ifelse(is.na(DEG_sex_bias), "No significant difference", DEG_sex_bias))
  
  # ---- step 1: one XWAS model per locus ----
  # loci = gene groups in box_groups[[trait]] (other genes = own locus); the locus model is
  # the model of the strongest cTWAS cell (max |Z|) across all supported cells of the locus
  grp <- box_groups[[trait_name]]
  sup$locus <- sup$Gene
  for (gi in seq_along(grp)) sup$locus[sup$Gene %in% grp[[gi]]] <- paste0("locus_", gi)
  sup <- sup %>%
    group_by(locus) %>%
    mutate(locus_model = gwas_model[which.max(abs(TWAS.Z_best))]) %>%
    ungroup()
  
  # ---- step 2: best cell per gene ----
  # prefer the supported cell whose sex-DEG result is concordant with the locus model
  # (eXCI/female model -> female-upregulated; male model -> male-upregulated;
  #  rXCI model -> no sex difference); ties broken by the strongest cTWAS (max |Z|)
  anno_src <- sup %>%
    mutate(deg_concordant = case_when(
      locus_model %in% c("eXCI", "female") & DEG_sex_bias == "Upregulated in female"  ~ TRUE,
      locus_model == "male"                & DEG_sex_bias == "Upregulated in male"    ~ TRUE,
      locus_model == "rXCI"                & DEG_sex_bias == "No significant difference" ~ TRUE,
      TRUE ~ FALSE)) %>%
    group_by(Gene) %>%
    arrange(desc(deg_concordant), desc(abs(TWAS.Z_best)), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(Gene, locus, locus_model, best_cell_type = cell_type, best_model = gwas_model,
           best_Z = TWAS.Z_best, SMR_SUPPORT, COLOC_SUPPORT)
  
  # best cell per gene (used for the right-side annotation and the black box)
  df_local$box_cell <- paste(df_local$Gene, df_local$x_axis, sep = "__") %in%
    paste(anno_src$Gene, anno_src$best_cell_type, sep = "__")
  
  # dagger flag per cell: only for prioritized cells (cTWAS FDR<0.05 AND SMR or COLOC support,
  # i.e. final_support), flagged when the sex-DEG result (nominal P < deg_p_cut) is concordant
  # with the locus XWAS model, same rule as the best-cell selection above
  # (eXCI/female -> female-up; male -> male-up; rXCI -> no sex difference)
  df_local <- df_local %>%
    left_join(anno_src %>% select(Gene, locus_model), by = "Gene") %>%
    left_join(df_deg_extra %>% select(Gene, cell_type, cell_DEG = DEG_sex_bias),
              by = c("Gene", "x_axis" = "cell_type")) %>%
    mutate(cell_prioritized = TWAS_SUPPORT %in% TRUE & (SMR_SUPPORT %in% TRUE | COLOC_SUPPORT %in% TRUE),
           deg_dagger = cell_prioritized & case_when(
             locus_model %in% c("eXCI", "female") & cell_DEG == "Upregulated in female"     ~ TRUE,
             locus_model == "male"                & cell_DEG == "Upregulated in male"       ~ TRUE,
             locus_model == "rXCI"                & cell_DEG == "No significant difference" ~ TRUE,
             TRUE ~ FALSE))
  
  # DEG sex bias taken at the annotated (best) cell type, so the right-side column matches that cell
  anno_src <- anno_src %>%
    left_join(df_deg_extra, by = c("Gene", "best_cell_type" = "cell_type")) %>%
    mutate(DEG_sex_bias = ifelse(is.na(DEG_sex_bias), "No significant difference", DEG_sex_bias))
  
  anno_src <- anno_src[match(rownames(m), anno_src$Gene), ]
  
  list(mat = m, df = df_local, anno = anno_src)
}
# ---- build one block per trait, drop empty ones ----
blocks <- lapply(target_traits_present, build_trait_block)
names(blocks) <- target_traits_present
present <- target_traits_present[!sapply(blocks, is.null)]
# ====== use row_split to isolate specified gene groups into their own slices ======
# each gene group (box_groups) per trait -> its own slice, with a gap from other genes
make_box_split <- function(gene_names, groups) {
  # groups: list of character vectors (gene groups to isolate for this trait); may be NULL
  n <- length(gene_names)
  if (n == 0) return(factor(integer(0)))
  if (is.null(groups)) groups <- list()
  # which box group each gene belongs to (NA if none)
  membership <- rep(NA_integer_, n)
  for (gi in seq_along(groups)) {
    membership[gene_names %in% groups[[gi]]] <- gi
  }
  # walk row order: adjacent genes "same group (incl. both non-group)" stay together; else new slice
  grp <- integer(n)
  g <- 1L; grp[1] <- 1L
  if (n > 1) {
    for (idx in 2:n) {
      prev <- membership[idx - 1]; cur <- membership[idx]
      same <- (is.na(prev) && is.na(cur)) ||
        (!is.na(prev) && !is.na(cur) && prev == cur)
      if (!same) g <- g + 1L
      grp[idx] <- g
    }
  }
  factor(grp, levels = unique(grp))
}
if (length(present) > 0) {
  # symmetric range: max |Z| across the three diseases (no clipping, keep strongest signal); 0-centered, symmetric legend
  all_abs <- unlist(lapply(present, function(t) abs(as.vector(blocks[[t]]$mat))))
  zmax <- max(all_abs, na.rm = TRUE)            # for outlier robustness use quantile(all_abs, 0.99, na.rm = TRUE)
  if (is.na(zmax) || zmax == 0) zmax <- 1
  col_fun_combined <- colorRamp2(c(-zmax, 0, zmax), c(twas_low_col, "white", twas_high_col))
  
  ht_list <- NULL
  for (k in seq_along(present)) {
    t  <- present[k]
    bk <- blocks[[t]]
    m_plot <- bk$mat          # no manual clipping; values beyond +/-zmax saturate via colorRamp2
    
    # this block's row_split: isolate box_groups[[t]] gene groups into their own slices
    split_fac <- make_box_split(rownames(m_plot), box_groups[[t]])
    nsl <- nlevels(split_fac)
    rt  <- rep("", nsl)
    rt[ceiling(nsl / 2)] <- t   # disease name only on the middle slice, avoid repeating per segment
    
    is_first <- (k == 1)
    is_last  <- (k == length(present))
    
    ha_right <- rowAnnotation(
      `Locus XWAS model`    = factor(bk$anno$locus_model,    levels = models),
      `DEG sex bias`        = factor(bk$anno$DEG_sex_bias,   levels = deg_levels),
      `SMR`                 = factor(ifelse(bk$anno$SMR_SUPPORT   %in% TRUE, "yes", "no"), levels = c("yes","no")),
      `COLOC`               = factor(ifelse(bk$anno$COLOC_SUPPORT %in% TRUE, "yes", "no"), levels = c("yes","no")),
      col = list(`Locus XWAS model`       = model_colors,
                 `DEG sex bias`           = deg_colors,
                 `SMR`                    = support_colors,
                 `COLOC`                  = support_colors),
      gp = gpar(col = "black", lwd = 0.8),
      annotation_name_gp = gpar(fontsize = 7), annotation_name_rot = 45,
      simple_anno_size = unit(8, "mm"), gap = unit(1, "mm"),
      show_annotation_name = is_last,
      show_legend = is_first,
      annotation_legend_param = list(
        `Locus XWAS model`       = list(title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
        `DEG sex bias`           = list(title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
        `SMR`                    = list(title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
        `COLOC`                  = list(title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 8))
      )
    )
    
    ht <- Heatmap(
      m_plot,
      name = if (is_first) "TWAS Z" else paste0("TWAS Z_", t),  # unique name per block, avoid duplicate-name warning
      col = col_fun_combined, na_col = na_color,
      cluster_rows = FALSE, cluster_columns = FALSE,
      row_split = split_fac,            # isolate box_groups gene groups into their own slices
      row_gap   = unit(0.8, "mm"),      # gap width, adjust as needed
      show_row_names = TRUE, row_names_gp = gpar(fontsize = 6), row_names_side = "left",
      show_column_names = is_last, column_names_gp = gpar(fontsize = 8),
      column_names_rot = 45,
      row_title = rt, row_title_gp = gpar(fontsize = 11, fontface = "bold"),
      right_annotation = ha_right,
      cell_fun = make_cf(m_plot, bk$df),
      show_heatmap_legend = is_first,
      heatmap_legend_param = list(
        direction = "horizontal",                 # horizontal at top
        title_gp  = gpar(fontsize = 8, fontface = "bold"),
        labels_gp = gpar(fontsize = 8)
      )
    )
    
    ht_list <- if (is.null(ht_list)) ht else ht_list %v% ht
  }
  
  # height auto-scales with total gene count, plus a bit per gap (bump manually if tight)
  total_rows <- sum(sapply(present, function(t) nrow(blocks[[t]]$mat)))
  n_gaps_total <- sum(sapply(present, function(t)
    nlevels(make_box_split(rownames(blocks[[t]]$mat), box_groups[[t]])) - 1L))
  # cairo_pdf (not pdf) so the Unicode dagger symbol is embedded correctly
  cairo_pdf("AD_PD_MSA_twas_smr_coloc_stacked.pdf",
            width = 11, height = max(6, 1.5 + total_rows * 0.3 + n_gaps_total * 0.08))
  draw(ht_list,
       heatmap_legend_side    = "top",     # TWAS Z legend on top
       annotation_legend_side = "right",   # annotation legend on the right
       column_title = NULL,                # set a custom title here if wanted
       column_title_gp = gpar(fontsize = 13, fontface = "bold"))
  dev.off()
}
# ==============================================================================
# 6. LocusZoom plots for the prioritized genes (AD / PD / MSA; Supplementary Figures)
#    -> locuszoom_plots/LCzoom_<gene>_<trait>_xwas_eqtl_eqtl_top_shared.pdf
#    One plot per trait+Gene (strongest |TWAS Z| among final_support rows).
#    Top panel: XWAS; following panels: cis-eQTL per cell type (top cell type first);
#    bottom: gene track (target gene highlighted). LD (r2) from the chrX QC bfile.
#    Index variant = lead eQTL variant of the top cell type; common y axis across eQTL panels.
# ==============================================================================

# ---- CONFIG (section 6) ----
lz_out_dir     <- "locuszoom_plots"
lz_xwas_dir    <- "../08_smr/8.2_xwas_clean/"   # folder holding <trait>_<model>.ma XWAS files
lz_bfile       <- "../00.1_X_chr_bfile_QC/chrX_QC_EU.0.G"   # LD reference
lz_plink       <- "plink1.9"
LZ_INDEX_SOURCES <- "eqtl_top"   # index variant: lead eQTL variant of the top cell type ("xwas_top" = XWAS lead variant)
LZ_Y_MODES       <- "shared"     # common y axis across eQTL panels ("free" = autoscale each panel)
SHOW_ONLY_PROTEIN_CODING <- TRUE          # FALSE = show all biotypes
KEEP_BIOTYPES <- c("protein_coding")      # e.g. c("protein_coding", "lncRNA")
# layout / size
PANEL_UNIT   <- 1.25   # height of each scatter panel
GENE_UNIT    <- 1.2    # height of the gene track
IN_PER_UNIT  <- 1      # inches per unit
PDF_WIDTH    <- 7      # page width
GENE_CEXTEXT <- 0.7    # gene label size
GENE_MAXROWS <- 7      # max stacked gene rows
GENE_HL_COL  <- "red"  # highlight colour for the target gene
TITLE_CEX    <- 0.72   # panel title size
LAB_CEX      <- 0.6    # y-axis title (-log10P) size
AXIS_CEX     <- 0.9    # y-axis tick label size
# panel margins (smaller = more compact)
MAR_TOP_XWAS <- 2      # XWAS panel top margin (two-line title)
MAR_TOP_EQTL <- 2      # eQTL panel top margin (one-line title)
MAR_BOTTOM   <- 0.5    # scatter panel bottom margin (no x axis)
MAR_GENE_TOP <- 0.8    # gene track top margin
MAR_LR       <- 4      # left/right margins (same for all panels so x axes align)

if (!dir.exists(lz_xwas_dir)) stop("XWAS .ma folder not found: ", lz_xwas_dir, " (edit lz_xwas_dir)")
dir.create(lz_out_dir, showWarnings = FALSE)
lz_temp_dir <- file.path(lz_out_dir, "temp")
dir.create(lz_temp_dir, showWarnings = FALSE)
lz_temp_abs <- normalizePath(lz_temp_dir)

pretty_model <- function(x) {
  map <- c(both_model_1 = "eXCI",
           both_model_2 = "rXCI",
           female_xchr_model_2 = "Female",
           male_xchr_model_1   = "Male")
  out <- unname(map[as.character(x)])
  ifelse(is.na(out), as.character(x), out)
}

safe_neglog10_max <- function(p) {
  p <- p[is.finite(p) & p > 0]
  if (length(p) == 0) return(1)
  -log10(min(p))
}

# genes to plot: final_support, AD/PD/MSA, one row per trait + Gene (max |TWAS Z|)
df_lz <- twas_smr_coloc_merged %>%
  dplyr::filter(final_support == TRUE, trait %in% target_traits) %>%
  dplyr::mutate(abs_best_TWAS_Z = abs(TWAS.Z_best)) %>%
  dplyr::arrange(dplyr::desc(abs_best_TWAS_Z)) %>%
  dplyr::distinct(trait, Gene, .keep_all = TRUE) %>%
  dplyr::select(cell_type, Gene, trait, gwas_model, eqtl_setting_coloc)
cat("LocusZoom: plotting", nrow(df_lz), "trait-gene pairs\n")

# cis-eQTL results
df_eqtl_all <- readRDS(eqtl_rds)

# rsID lookup, used only to display the index SNP in the panel title
# (variants are matched / LD-computed by chr:pos:ref:alt ID; falls back to the ID if no rsID)
rsid_map <- readRDS("../00.3_rsid_annotation/rsid_gene_annotation.rds") %>%
  dplyr::select(ID, rsid) %>%
  dplyr::distinct(ID, .keep_all = TRUE)

# ---- plotting function: one multi-panel LocusZoom PDF for a given index source + y mode ----
make_lczoom <- function(index_source, y_mode,
                        gene_i, trait_i, eqtl_model, xwas_model, top_ct,
                        df_eqtl, df_xwas, lz_cell_types,
                        seqname_val, xrange_min, xrange_max,
                        eqtl_plots, i) {
  
  ## ---- index SNP ----
  if (index_source == "eqtl_top") {
    df_top <- dplyr::filter(df_eqtl, cell_type == top_ct)
    if (nrow(df_top) == 0) df_top <- df_eqtl
    index_snp <- df_top$ID[which.min(df_top$pvalue)]
    index_lab <- paste0("eQTL top (", top_ct, ")")
  } else if (index_source == "xwas_top") {
    index_snp <- df_xwas$SNP[which.min(df_xwas$p)]
    index_lab <- paste0(trait_i, " XWAS top")
  } else {
    stop("index_source must be 'eqtl_top' or 'xwas_top'")
  }
  # index SNP label for the title: rsID if available, otherwise the chr:pos:ref:alt ID
  index_show <- rsid_map$rsid[match(index_snp, rsid_map$ID)]
  index_show <- ifelse(is.na(index_show) | index_show == "", index_snp, index_show)
  
  ## ---- plink LD (r2 vs index SNP) ----
  ld_tag <- paste0(i, "_", index_source)
  snp_list_file <- file.path(lz_temp_abs, paste0("snp_list", ld_tag, ".txt"))
  write.table(unique(df_eqtl$ID), snp_list_file,
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  command <- paste(lz_plink,
                   "--bfile", lz_bfile,
                   "--allow-no-sex",
                   "--r2",
                   "--ld-snp", index_snp,
                   "--ld-window 99999",
                   "--ld-window-kb 99999",
                   "--ld-window-r2 0",
                   "--extract", snp_list_file,
                   "--keep-allele-order",
                   "--out", file.path(lz_temp_abs, paste0("ld_res", ld_tag)),
                   sep = " ")
  system(command)
  
  ld_file <- file.path(lz_temp_abs, paste0("ld_res", ld_tag, ".ld"))
  if (!file.exists(ld_file)) { cat("  !! LD file missing, skipped\n"); return(invisible(NULL)) }
  ld_df <- fread(ld_file)
  ld_df <- dplyr::select(ld_df, rsid = SNP_B, R2)
  
  ## ---- XWAS plotting table ----
  xwas_plot <- df_xwas %>%
    dplyr::transmute(chrom = seqname_val,
                     pos   = BP,
                     rsid  = SNP,
                     other_allele  = A2,
                     effect_allele = A1,
                     p     = p,
                     beta  = b,
                     se    = se) %>%
    dplyr::left_join(ld_df, by = "rsid")
  
  ## ---- eQTL plotting tables (attach LD) ----
  eqtl_plots_ld <- lapply(eqtl_plots, function(d)
    dplyr::left_join(dplyr::select(d, -any_of("R2")), ld_df, by = "rsid"))
  
  ## ---- biotype filter ----
  filter_TX <- function(loc) {
    if (SHOW_ONLY_PROTEIN_CODING) loc$TX <- dplyr::filter(loc$TX, gene_biotype %in% KEEP_BIOTYPES)
    loc
  }
  
  ## ---- locus: XWAS ----
  loc_xwas <- locus(xwas_plot, gene = NULL,
                    seqname = seqname_val,
                    xrange  = c(xrange_min, xrange_max),
                    index_snp = index_snp, LD = "R2",
                    ens_db = "EnsDb.Hsapiens.v86")
  loc_xwas <- tryCatch(link_recomb(loc_xwas, genome = "hg38"),
                       error = function(e) { cat("  (recombination track unavailable, skipped)\n"); loc_xwas })
  loc_xwas <- filter_TX(loc_xwas)
  
  ## ---- locus: eQTL per cell type ----
  loc_eqtl_list <- lapply(lz_cell_types, function(ct) {
    dat <- eqtl_plots_ld[[ct]]
    idx_use <- if (index_snp %in% dat$rsid) index_snp else NULL
    lo <- locus(dat, gene = NULL,
                seqname = seqname_val,
                xrange  = c(xrange_min, xrange_max),
                index_snp = idx_use, LD = "R2",
                ens_db = "EnsDb.Hsapiens.v86")
    lo$recomb <- loc_xwas$recomb
    lo
  })
  names(loc_eqtl_list) <- lz_cell_types
  
  ## ---- y limits ----
  ylim_xwas <- c(0, safe_neglog10_max(xwas_plot$p) + 0.5)
  ymax_each <- sapply(lz_cell_types, function(ct)
    safe_neglog10_max(eqtl_plots_ld[[ct]]$p))
  ylim_eqtl_shared <- c(0, max(ymax_each) + 0.5)
  
  ## ---- draw ----
  n_ct      <- length(lz_cell_types)
  n_panels  <- 1 + n_ct
  row_units <- c(rep(PANEL_UNIT, n_panels), GENE_UNIT)
  pdf_h     <- sum(row_units) * IN_PER_UNIT
  
  out_pdf <- file.path(lz_out_dir, paste0("LCzoom_", gene_i, "_", trait_i,
                                          "_xwas_eqtl_", index_source, "_", y_mode, ".pdf"))
  pdf(out_pdf, width = PDF_WIDTH, height = pdf_h)
  layout(matrix(1:(n_panels + 1), ncol = 1), heights = row_units)
  
  # XWAS on top
  par(mar = c(MAR_BOTTOM, MAR_LR, MAR_TOP_XWAS, MAR_LR))
  scatter_plot(loc_xwas, ylim = ylim_xwas, labels = NULL,
               legend_pos = "topleft", xaxt = "n", xlab = "",
               cex.lab = LAB_CEX, cex.axis = AXIS_CEX)
  title(main = paste0(trait_i, " XWAS (", pretty_model(xwas_model), ")  -  ", gene_i,
                      "\nindex [", index_lab, "]: ", index_show),
        line = 0.6, cex.main = TITLE_CEX)
  
  # eQTL per cell type
  for (k in seq_len(n_ct)) {
    ct <- lz_cell_types[k]
    ylim_ct <- if (y_mode == "shared") ylim_eqtl_shared else c(0, ymax_each[[ct]] + 0.5)
    
    par(mar = c(MAR_BOTTOM, MAR_LR, MAR_TOP_EQTL, MAR_LR))
    scatter_plot(loc_eqtl_list[[ct]], ylim = ylim_ct, labels = NULL,
                 legend_pos = if (k == 1) "topleft" else NULL,
                 xaxt = "n", xlab = "",
                 cex.lab = LAB_CEX, cex.axis = AXIS_CEX)
    tag <- if (ct == top_ct) "  [top]" else ""
    title(main = paste0(gene_i, " eQTL - ", ct, " (", pretty_model(eqtl_model), ")", tag),
          line = 0.6, cex.main = TITLE_CEX)
  }
  
  # gene track: target gene highlighted
  par(mar = c(3, MAR_LR, MAR_GENE_TOP, MAR_LR))
  genetracks(loc_xwas,
             cex.text      = GENE_CEXTEXT,
             maxrows       = GENE_MAXROWS,
             highlight     = gene_i,
             highlight_col = GENE_HL_COL)
  dev.off()
  cat("  written: ", out_pdf, "\n", sep = "")
}

# ---- main loop ----
for (i in seq_len(nrow(df_lz))) {
  
  gene_i      <- df_lz$Gene[i]
  trait_i     <- df_lz$trait[i]
  eqtl_model  <- df_lz$eqtl_setting_coloc[i]
  xwas_model  <- df_lz$gwas_model[i]
  top_ct      <- df_lz$cell_type[i]
  
  cat("\n==== [", i, "/", nrow(df_lz), "] ", gene_i, " | ", trait_i,
      " | eqtl:", eqtl_model, " | gwas:", xwas_model, " ====\n", sep = "")
  
  ## ---- cis-eQTL for this gene / eQTL model ----
  df_eqtl <- dplyr::filter(df_eqtl_all, gene == gene_i, setting == eqtl_model)
  
  ## ---- XWAS summary statistics ----
  xwas_file <- file.path(lz_xwas_dir, paste0(trait_i, "_", xwas_model, ".ma"))
  if (!file.exists(xwas_file)) { cat("  !! XWAS file missing:", xwas_file, ", skipped\n"); next }
  df_xwas <- fread(xwas_file)
  df_xwas$BP <- as.integer(sub("^[^:]*:([^:]*):.*", "\\1", df_xwas$SNP))
  df_xwas <- dplyr::filter(df_xwas, BP >= min(df_eqtl$POS), BP <= max(df_eqtl$POS))
  
  ## ---- shared SNPs ----
  df_xwas <- dplyr::semi_join(df_xwas, df_eqtl, by = c("SNP" = "ID"))
  df_eqtl <- dplyr::semi_join(df_eqtl, df_xwas, by = c("ID"  = "SNP"))
  if (nrow(df_eqtl) == 0 || nrow(df_xwas) == 0) {
    cat("  !! no shared SNPs, skipped\n"); next
  }
  
  ## ---- cell-type order: top cell type first, others by min p ----
  ct_rank <- df_eqtl %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarise(minp = min(pvalue), .groups = "drop")
  lz_cell_types <- c(top_ct, setdiff(ct_rank$cell_type[order(ct_rank$minp)], top_ct))
  lz_cell_types <- lz_cell_types[lz_cell_types %in% df_eqtl$cell_type]
  cat("  # of cell types:", length(lz_cell_types),
      "->", paste(lz_cell_types, collapse = ", "), "\n")
  
  seqname_val <- as.character(df_eqtl[["#CHROM"]][1])
  xrange_min  <- min(df_eqtl$POS, na.rm = TRUE)
  xrange_max  <- max(df_eqtl$POS, na.rm = TRUE)
  
  ## ---- per-cell-type eQTL plotting tables (LD attached inside the function) ----
  eqtl_plots <- lapply(lz_cell_types, function(ct) {
    dplyr::filter(df_eqtl, cell_type == ct) %>%
      dplyr::transmute(chrom = `#CHROM`,
                       pos   = POS,
                       rsid  = ID,
                       other_allele  = REF,
                       effect_allele = ALT,
                       p     = pvalue,
                       beta  = BETA,
                       se    = SE) %>%
      dplyr::arrange(p)
  })
  names(eqtl_plots) <- lz_cell_types
  
  ## ---- index-variant / y-axis versions (see CONFIG) ----
  for (src in LZ_INDEX_SOURCES) {
    for (ym in LZ_Y_MODES) {
      make_lczoom(index_source = src, y_mode = ym,
                  gene_i = gene_i, trait_i = trait_i,
                  eqtl_model = eqtl_model, xwas_model = xwas_model, top_ct = top_ct,
                  df_eqtl = df_eqtl, df_xwas = df_xwas,
                  lz_cell_types = lz_cell_types,
                  seqname_val = seqname_val,
                  xrange_min = xrange_min, xrange_max = xrange_max,
                  eqtl_plots = eqtl_plots, i = i)
    }
  }
}
message("Done: LocusZoom plots in ", lz_out_dir)
