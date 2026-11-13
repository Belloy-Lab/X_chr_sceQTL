# Sex DEG (female vs male) on X-chromosome eQTL input genes per cell type, plus figures; the case/control-adjusted model is computed in-memory only for the logFC correlation
rm(list = ls())
library(tidyverse)
library(edgeR)
library(limma)
library(data.table)
library(HGNChelper)
library(EnhancedVolcano)
library(patchwork)
library(ggrepel)
library(scales)
library(ggsignif)
library(ComplexHeatmap)
library(ggtext)
library(cowplot)
# ==============================================================================
# CONFIG: inputs and parameters (edit here only)
# ==============================================================================
metadata_xlsx <- "path/to/metadata.xlsx"
pseudobulk_dir <- "path/to/pseudobulk/"
cov_pheno_dir  <- "../01.2_cov_pheno/"
# cell types in the order of the pseudobulk rds files (alphabetical, from dir())
cell_types_analysis <- c("Ast", "End", "Exc", "Inh", "Mic", "Oli", "OPC")
# cell-type display order for all figures
cell_types <- c("Exc", "Inh", "Ast", "Oli", "OPC", "Mic", "End")
# main DEG results file
deg_xlsx <- "DEG_female_vs_male_X_chr_gene_non_adj_disease.xlsx"
# ==============================================================================
# 1. Metadata and analysis table
# ==============================================================================
meta_data <- readxl::read_excel(metadata_xlsx)
meta_data$msex <- factor(meta_data$msex, levels = c("1", "0"), labels = c("Male", "Female"))
meta_data$case_control <- ifelse(meta_data$cogdx == 1, "control",
                                 ifelse(meta_data$cogdx %in% 2:6, "case", NA))  # cogdx: 1=control, 2-6=case
meta_data <- select(meta_data, specimenID, batch, Study, age_death, msex, pmi, case_control)

count_files   <- dir(pseudobulk_dir, pattern = "rds")
analysis_list <- tibble(
  cell_types = cell_types_analysis,
  count_path = paste0(pseudobulk_dir, count_files)
)
analysis_list$cov_path <- paste0(cov_pheno_dir, "both_", analysis_list$cell_types, ".cov.txt")
# ==============================================================================
# 2. DEG function (sex effect); optionally adjust for case/control
# ==============================================================================
run_deg <- function(adjust_case_control = FALSE) {
  res_list <- vector("list", nrow(analysis_list))
  for (i in seq_len(nrow(analysis_list))) {
    cell_type <- analysis_list$cell_types[i]
    cat("Run DEG analysis (sex) in", cell_type,
        if (adjust_case_control) "[case/control-adjusted]" else "", "\n")
    df <- as.data.frame(read_rds(analysis_list$count_path[i]))
    
    # drop samples with cell number < 10 (NA in Fujita cov files)
    valid_samples <- fread(analysis_list$cov_path[i])$IID
    meta <- meta_data %>%
      filter(specimenID %in% valid_samples,
             !is.na(Study), !is.na(msex), !is.na(age_death), !is.na(batch))
    meta$pmi[is.na(meta$pmi)] <- median(meta$pmi, na.rm = TRUE)
    df <- select(df, meta$specimenID)
    
    # design (optionally adjust for case/control)
    design <- if (adjust_case_control) {
      model.matrix(~ msex + Study + pmi + age_death + batch + case_control, data = meta)
    } else {
      model.matrix(~ msex + Study + pmi + age_death + batch, data = meta)
    }
    
    # low-expression filter + TMM + voom
    dge <- DGEList(counts = df)
    dge <- dge[filterByExpr(dge), , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge, method = "TMM")
    v <- voom(dge, design)
    
    # drop genes with mean log2CPM < 2.0, then re-voom
    logcpm <- v$E[apply(v$E, 1, mean) > 2.0, ]
    dge$counts <- dge$counts[rownames(logcpm), , drop = FALSE]
    v <- voom(dge, design)
    
    # DE test (coef 2 = msex)
    fit <- eBayes(lmFit(v, design))
    res <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
    res <- rownames_to_column(res, var = "Gene")
    res$gene_cell <- paste0(res$Gene, "_", cell_type)
    tbl <- summary(meta$msex)
    res$compare_summary <- paste(names(tbl), tbl, sep = ": ", collapse = ", ")
    
    # update gene symbols; for duplicates keep the HGNC-approved one
    gsu <- checkGeneSymbols(res$Gene)
    res$updated_Gene_symbol <- gsu$Suggested.Symbol
    res$Approved <- gsu$Approved
    if (any(duplicated(res$updated_Gene_symbol))) {
      res <- res %>%
        group_by(updated_Gene_symbol) %>%
        arrange(desc(Approved == TRUE)) %>%
        slice(1) %>%
        ungroup()
    }
    res <- res[!duplicated(res$updated_Gene_symbol), ]
    res$gene_cell <- paste0(res$updated_Gene_symbol, "_", cell_type)
    res$cell_type <- cell_type
    colnames(res)[1] <- "Old_gene_symbol"
    
    # keep only X-chr genes tested in eQTL for this cell type
    pheno <- fread(paste0(cov_pheno_dir, "both_", cell_type, ".pheno_chrX.tsv"))
    pheno$FID <- NULL; pheno$IID <- NULL
    res <- filter(res, updated_Gene_symbol %in% colnames(pheno))
    
    # FDR restricted to the input genes
    res$adj.P.Val <- p.adjust(res$P.Value, method = "BH")
    res_list[[i]] <- res
  }
  rbindlist(res_list, fill = TRUE)
}
# ==============================================================================
# 3. Non-adjusted DEG -> results table  [OUTPUT: DEG_female_vs_male_X_chr_gene_non_adj_disease.xlsx]
# ==============================================================================
DEG_all <- run_deg(adjust_case_control = FALSE)
writexl::write_xlsx(DEG_all, deg_xlsx)
DEG_all <- readxl::read_excel(deg_xlsx)   # re-read to normalize column types for the plots below

# ==============================================================================
# 4. Volcano plot (version 2: same y axis, compressed x/y axes)
#    [OUTPUT: volcanoplot_..._fdr_significant_version_2.pdf]
# ==============================================================================
# ============================== Tunable parameters ==============================
# ---- x compression (non-End): x <= brk linear, x > brk compressed ----
brk      <- 1                              # x breakpoint
cmp      <- 0.05                           # x compression factor: smaller = stronger
x_breaks <- c(-2, -1, 0, 1, 2, 3, 4, 5, 6) # x ticks for non-End (real logFC values)
# ---- y compression (non-End): y <= y_brk linear, y > y_brk compressed ----
y_brk    <- 100                            # y breakpoint
y_cmp    <- 0.05                           # y compression factor: smaller = stronger
y_breaks <- c(0, 50, 100, 200, 300)        # y ticks for non-End
y_max    <- 300                            # shared y upper limit for the six non-End panels
# ---- End only: real coordinates, no compression ----
end_x_breaks <- c(0.5)                     # x ticks for End (only 0.5)
end_x_right  <- 0.6                         # fixed right end of End x axis
# ---- panel spacing (patchwork): t r b l, in pt ----
panel_margin_t <- 2.5
panel_margin_r <- 3
panel_margin_b <- 2.5
panel_margin_l <- 2.5
title_margin_b <- 1.5      # bottom margin between each panel title and its plot
# ---- axis-title to axis distance ----
xtitle_margin_t <-  0      # x-title to x-axis (smaller = closer)
ytitle_margin_r <- -8      # y-title to y-axis (smaller = closer)
# ---- axis tick label size ----
axis_text_size <- 15       # x/y tick label size (passed to EnhancedVolcano axisLabSize)
# ===============================================================================
# x compression transform
squish_x <- scales::trans_new(
  name      = "squish_x",
  transform = function(x) ifelse(x <= brk, x, brk + (x - brk) * cmp),
  inverse   = function(x) ifelse(x <= brk, x, brk + (x - brk) / cmp)
)
# y compression transform
squish_y <- scales::trans_new(
  name      = "squish_y",
  transform = function(y) ifelse(y <= y_brk, y, y_brk + (y - y_brk) * y_cmp),
  inverse   = function(y) ifelse(y <= y_brk, y, y_brk + (y - y_brk) / y_cmp)
)
plot_list <- list()
for (ct in cell_types) {
  cat("generating volcano plot for:", ct, "\n")
  
  df_plot <- dplyr::filter(DEG_all, cell_type == ct)
  
  df_plot$color_group <- "Not sig"
  df_plot$color_group[df_plot$adj.P.Val < 0.05 & df_plot$logFC > 0] <- "Female high expression"
  df_plot$color_group[df_plot$adj.P.Val < 0.05 & df_plot$logFC < 0] <- "Male high expression"
  
  keyvals <- ifelse(df_plot$adj.P.Val >= 0.05, "grey80",
                    ifelse(df_plot$logFC > 0, "red", "blue"))
  names(keyvals)[keyvals == "red"]    <- "Female-upregulated DEGs"
  names(keyvals)[keyvals == "blue"]   <- "Male-upregulated DEGs"
  names(keyvals)[keyvals == "grey80"] <- "No differential expression"
  
  p <- EnhancedVolcano(
    df_plot,
    lab = df_plot$updated_Gene_symbol,
    x = "logFC",
    y = "adj.P.Val",
    pCutoff = 0.05,
    FCcutoff = 0,
    colCustom = keyvals,
    pointSize = 1.5,
    labSize = 3,
    axisLabSize = axis_text_size,        # axis tick label size (native param)
    title = ct,
    subtitle = NULL,
    drawConnectors = TRUE,
    caption = NULL,
    max.overlaps = 20,
    ylab = expression(-Log[10]("adj.P.Val"))
  )
  
  if (ct != "End") {
    # non-End: compress x and y, unify y to y_max
    p <- p + ggplot2::scale_x_continuous(
      trans  = squish_x,
      breaks = x_breaks,
      expand = ggplot2::expansion(mult = 0.05)
    )
    p <- p + ggplot2::scale_y_continuous(
      trans  = squish_y,
      breaks = y_breaks,
      limits = c(0, y_max),
      expand = ggplot2::expansion(mult = 0.02)
    )
  } else {
    # End: real x/y coordinates, no compression; fixed right end, single x tick
    end_p   <- -log10(df_plot$adj.P.Val); end_p <- end_p[is.finite(end_p)]
    x_left  <- min(df_plot$logFC, na.rm = TRUE) - 0.1   # left end auto from data
    p <- p +
      ggplot2::scale_x_continuous(
        breaks = end_x_breaks,
        expand = ggplot2::expansion(mult = 0.02)
      ) +
      ggplot2::coord_cartesian(
        xlim = c(x_left, end_x_right),
        ylim = c(0, max(end_p, na.rm = TRUE) + 5)
      )
  }
  
  plot_list[[ct]] <- p
}
combined_plot <- wrap_plots(plot_list, ncol = 3) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(margin = ggplot2::margin(t = xtitle_margin_t)),
    axis.title.y = element_text(margin = ggplot2::margin(r = ytitle_margin_r)),
    plot.margin  = ggplot2::margin(panel_margin_t, panel_margin_r,
                                   panel_margin_b, panel_margin_l),
    plot.title   = element_text(margin = ggplot2::margin(b = title_margin_b))
  )
ggsave("volcanoplot_DEG_female_vs_male_X_chr_gene_adj_disease_fdr_significant_version_2.pdf",
       combined_plot, width = 17, height = 15.8, dpi = 300)

# ==============================================================================
# 5. Stacked bar: DEG split by direction (Female-up / Male-up / Non-DEG) (top)
#    and Fujita 2024 cell-type proportion (bottom)
#    [OUTPUT: barplot_DEG_across_model_split_direction.pdf]
# ==============================================================================
# ============================== Tunable parameters ==============================
# cell-type proportion (%) of all nuclei, from Fujita et al. 2024 (mean proportion; "Others" 1.1% not shown)
cell_type_prop   <- c(Exc = 41.1, Inh = 15.2, Ast = 13.8, Oli = 19.6, OPC = 3.6, Mic = 5.0, End = 0.6)
show_prop_label  <- FALSE   # TRUE: print the % value under each proportion bar
count_text_size  <- 4       # size of the count labels in the top panel
count_pad_bottom <- 0.04    # extra space below 0 in the top panel (fraction of y range),
# so labels of thin bottom segments are not cut off
# ===============================================================================
# total input genes per cell type
n_gene_total <- as.data.frame(table(DEG_all$cell_type))
colnames(n_gene_total) <- c("cell_type", "n_total")
# significant DEGs split by direction
DEG_sig_dir <- DEG_all %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::mutate(
    type = dplyr::case_when(
      logFC > 0 ~ "Female-upregulated",
      logFC < 0 ~ "Male-upregulated",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(type)) %>%
  dplyr::group_by(cell_type, type) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop")
# non-significant (Non-DEG) per cell type = total - sig
n_sig_per_ct <- DEG_sig_dir %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(n_sig = sum(count), .groups = "drop")
df_nonsig <- n_gene_total %>%
  dplyr::left_join(n_sig_per_ct, by = "cell_type") %>%
  dplyr::mutate(
    n_sig = ifelse(is.na(n_sig), 0, n_sig),
    count = n_total - n_sig,
    type = "Non-DEG"
  ) %>%
  dplyr::select(cell_type, type, count)
df_long2 <- dplyr::bind_rows(DEG_sig_dir, df_nonsig)
df_long2$cell_type <- factor(df_long2$cell_type, levels = cell_types)
# stacking order: Female above Male
df_long2$type <- factor(
  df_long2$type,
  levels = c("Non-DEG", "Female-upregulated", "Male-upregulated")
)

# top panel: DEG counts (x labels hidden; shown on the bottom panel)
# text layer is drawn last so labels sit on top of the bars
p_top <- ggplot(df_long2, aes(x = cell_type, y = count, fill = type)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.5) +   # alpha to match the volcano plots
  scale_fill_manual(
    values = c(
      "Female-upregulated" = "red",
      "Male-upregulated"   = "blue",
      "Non-DEG"            = "grey80"
    ),
    breaks = c("Female-upregulated", "Male-upregulated", "Non-DEG")
  ) +
  scale_y_continuous(expand = expansion(mult = c(count_pad_bottom, 0.05))) +
  coord_cartesian(clip = "off") +                          # do not clip labels at the panel edge
  labs(x = NULL, y = "Number of genes", fill = NULL) +
  theme_classic(base_size = 14) +
  theme(
    axis.text    = element_text(size = 12, color = "black"),
    axis.title   = element_text(size = 14),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    legend.text  = element_text(size = 12),
    legend.position = "right",
    plot.margin  = margin(t = 5, r = 10, b = 0, l = 10)
  ) +
  geom_text(
    aes(label = ifelse(count > 0, count, "")),
    position = position_stack(vjust = 0.5),
    size = count_text_size,
    color = "black"
  )

# bottom panel: cell-type proportion, y axis pointing downwards
# transparent background so it never covers the top panel's labels
prop_df <- data.frame(cell_type = factor(cell_types, levels = cell_types),
                      prop      = cell_type_prop[cell_types])
p_bottom <- ggplot(prop_df, aes(x = cell_type, y = prop)) +
  geom_bar(stat = "identity", width = 0.7, fill = "grey70") +
  { if (show_prop_label) geom_text(aes(label = paste0(prop, "%")), vjust = -0.6, size = 4, color = "black") } +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_y_reverse(breaks = c(0, 10, 20, 30, 40),
                  labels = function(x) ifelse(x == 0, "", paste0(x, "%")),
                  expand = expansion(mult = c(0.15, 0))) +
  labs(x = NULL, y = "Cell-type proportion") +
  theme_classic(base_size = 14) +
  theme(
    axis.text   = element_text(size = 12, color = "black"),
    axis.title  = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13, color = "black"),
    axis.line.x = element_blank(),
    plot.background  = element_blank(),
    panel.background = element_blank(),
    plot.margin = margin(t = 0, r = 10, b = 5, l = 10)
  )

# stack: top ~2/3, bottom ~1/3; shared legend on the right
p <- p_top / p_bottom + plot_layout(heights = c(2, 1), guides = "collect")
print(p)
ggsave("barplot_DEG_across_model_split_direction.pdf", plot = p,
       width = 7.8, height = 5, dpi = 600)

# ==============================================================================
# 6. Percentage of female- vs male-biased DEGs per cell type + chi-square
#    [OUTPUT: barplot_DEG_female_male_percentage_by_cell_type.pdf,
#             chisq_Mic_vs_others_DEG_direction.xlsx]
# ==============================================================================
DEG_direction <- DEG_all %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::mutate(
    DEG_direction = dplyr::case_when(
      logFC > 0 ~ "Female-upregulated",
      logFC < 0 ~ "Male-upregulated",
      TRUE ~ "No direction"
    )
  ) %>%
  dplyr::filter(DEG_direction %in% c("Female-upregulated", "Male-upregulated"))
DEG_direction_summary <- DEG_direction %>%
  dplyr::group_by(cell_type, DEG_direction) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(cell_type) %>%
  dplyr::mutate(
    percent = n / sum(n) * 100,
    label = paste0(round(percent, 1), "%\n(n=", n, ")")
  ) %>%
  dplyr::ungroup()
lv <- cell_types
DEG_direction_summary$cell_type <- factor(DEG_direction_summary$cell_type, levels = lv)

### Chi-square test: reference cell type vs each of the others (female:male DEG ratio)
ref_cell     <- "Mic"
use_adjusted <- T   # if TRUE, star labels use FDR-adjusted p-values
get_stars <- function(p) {
  if (is.na(p))  return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}
# wide table: one row per cell_type, Female / Male count columns (chi-square on counts, not percentages)
DEG_wide <- DEG_direction_summary %>%
  dplyr::select(cell_type, DEG_direction, n) %>%
  tidyr::pivot_wider(names_from = DEG_direction, values_from = n, values_fill = 0) %>%
  as.data.frame()
# fill missing columns / cell types so a missing direction doesn't error
for (col in c("Female-upregulated", "Male-upregulated")) {
  if (!col %in% colnames(DEG_wide)) DEG_wide[[col]] <- 0
}
miss_ct <- setdiff(lv, as.character(DEG_wide$cell_type))
if (length(miss_ct) > 0) {
  DEG_wide <- dplyr::bind_rows(
    DEG_wide,
    data.frame(cell_type = miss_ct,
               `Female-upregulated` = 0,
               `Male-upregulated`   = 0,
               check.names = FALSE)
  )
}
rownames(DEG_wide) <- as.character(DEG_wide$cell_type)
exc_F <- DEG_wide[ref_cell, "Female-upregulated"]
exc_M <- DEG_wide[ref_cell, "Male-upregulated"]
other_cells <- lv[lv != ref_cell]
chisq_results <- lapply(other_cells, function(ct) {
  ct_F <- DEG_wide[ct, "Female-upregulated"]
  ct_M <- DEG_wide[ct, "Male-upregulated"]
  m <- matrix(
    c(exc_F, exc_M,
      ct_F,  ct_M),
    nrow = 2, byrow = TRUE,
    dimnames = list(c(ref_cell, ct), c("Female", "Male"))
  )
  test <- chisq.test(m)   # 2x2 with Yates correction; use fisher.test(m) for expected counts < 5
  data.frame(
    comparison     = paste0(ref_cell, "_vs_", ct),
    cell_type      = ct,
    exc_F          = exc_F,
    exc_M          = exc_M,
    exc_female_pct = exc_F / (exc_F + exc_M) * 100,
    ct_F           = ct_F,
    ct_M           = ct_M,
    ct_female_pct  = ct_F / (ct_F + ct_M) * 100,
    chisq          = unname(test$statistic),
    df             = unname(test$parameter),
    p_value        = test$p.value
  )
}) %>% bind_rows()
chisq_results$p_adj <- p.adjust(chisq_results$p_value, method = "fdr")
print(chisq_results)
writexl::write_xlsx(chisq_results, "chisq_Mic_vs_others_DEG_direction.xlsx")

### Prepare brackets (ggsignif manual mode) - uniform size/spacing
pos <- setNames(seq_along(lv), lv)
p_for_star <- if (use_adjusted) chisq_results$p_adj else chisq_results$p_value
star_labels <- vapply(p_for_star, get_stars, character(1))
max_h  <- 100                 # percentage stack tops out at 100
base_y <- max_h * 1.02        # just above the bar top (102)
step_y <- max_h * 0.05        # fixed spacing between brackets
bracket_df <- data.frame(
  cell_type  = other_cells,
  xmin       = pos[ref_cell],
  xmax       = pos[other_cells],
  annotation = star_labels,
  stringsAsFactors = FALSE
)
# order by span (small to large), then assign increasing heights (uniform spacing)
bracket_df <- bracket_df[order(bracket_df$xmax), ]
bracket_df$y_position <- base_y + step_y * (seq_len(nrow(bracket_df)) - 1)

### Plot
p <- ggplot(DEG_direction_summary, aes(x = cell_type, y = percent, fill = DEG_direction)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.5) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "black"
  ) +
  geom_signif(
    annotations = bracket_df$annotation,
    xmin        = bracket_df$xmin,
    xmax        = bracket_df$xmax,
    y_position  = bracket_df$y_position,
    tip_length  = 0.005,
    size        = 0.3,
    textsize    = 2.5,
    vjust       = 0.7
  ) +
  scale_fill_manual(
    values = c(
      "Female-upregulated" = "red",
      "Male-upregulated" = "blue"
    )
  ) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
  coord_cartesian(
    ylim = c(0, max(bracket_df$y_position) + step_y),  # leave room for brackets
    clip = "off"
  ) +
  labs(x = NULL, y = "Percentage of significant DEGs", fill = NULL) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    legend.position = "top",
    legend.direction = "horizontal",
    plot.margin = margin(t = 40, r = 10, b = 10, l = 10)   # top padding so stars are not clipped
  )
print(p)
ggsave(
  "barplot_DEG_female_male_percentage_by_cell_type.pdf",
  plot = p,
  width = 4.5,
  height = 5.4,
  dpi = 300
)
# ==============================================================================
# 7. Cell-type-specific vs shared X-chromosome DEGs
#    [OUTPUT: DEG_gene_sharing_combined.pdf]
#    A) histogram of gene sharing (purple = concordant / gold = discordant / grey = specific)
#    B) heatmap of shared genes (>= 2 cell types)
#    C) heatmap of cell-type-specific genes (significant in exactly 1 cell type)
#    colours: red = Female-up (sig) / blue = Male-up (sig) / mid grey = tested not sig / light grey = not tested
# ==============================================================================
gene_col <- "updated_Gene_symbol"
stopifnot(gene_col %in% colnames(DEG_all))

# significant DEG, direction, de-dup per (gene, cell_type) keeping most significant
DEG_sig <- DEG_all %>%
  dplyr::mutate(gene = .data[[gene_col]]) %>%
  dplyr::filter(adj.P.Val < 0.05, !is.na(logFC), logFC != 0) %>%
  dplyr::mutate(direction = ifelse(logFC > 0, "Female", "Male")) %>%
  dplyr::group_by(gene, cell_type) %>%
  dplyr::slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

# gene classification: cell-type-specific / shared-concordant / shared-discordant
gene_class <- DEG_sig %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    n_celltype = dplyr::n_distinct(cell_type),
    celltypes  = paste(cell_type[order(match(cell_type, lv))], collapse = ", "),
    n_female   = sum(direction == "Female"),
    n_male     = sum(direction == "Male"),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    category = dplyr::case_when(
      n_celltype == 1             ~ "Cell-type-specific",
      n_female == 0 | n_male == 0 ~ "Shared (concordant)",
      TRUE                        ~ "Shared (discordant)"
    )
  ) %>%
  dplyr::arrange(dplyr::desc(n_celltype), category, gene)
cat("\n===== gene classification summary =====\n"); print(table(gene_class$category))

# proportions: specific vs shared, and direction concordance within shared
prop_specific_shared <- gene_class %>%
  dplyr::mutate(group = ifelse(n_celltype == 1, "Cell-type-specific", "Shared")) %>%
  dplyr::count(group) %>%
  dplyr::mutate(percent = round(n / sum(n) * 100, 1))
prop_shared_dir <- gene_class %>%
  dplyr::filter(n_celltype >= 2) %>%
  dplyr::count(category) %>%
  dplyr::mutate(percent_within_shared = round(n / sum(n) * 100, 1))
cat("\n===== specific vs shared =====\n"); print(prop_specific_shared)
cat("\n===== shared: concordant vs discordant =====\n"); print(prop_shared_dir)

# colours + status-matrix helper (shared by B and C)
col_map <- c("Female"    = "red",
             "Male"      = "blue",
             "NS"        = "#9E9E9E",   # tested but not significant
             "NotTested" = "#ECECEC")   # gene x cell with no result at all
legend_param <- list(
  at     = c("Female", "Male", "NS", "NotTested"),
  labels = c("Female-upregulated DEGs", "Male-upregulated DEGs",
             "No differential expression", "Not tested")
)
legend_param_top <- c(legend_param, list(direction = "horizontal", nrow = 1))
# build a gene x 7-cell-type 4-state status matrix for a given gene set
build_status_mat <- function(genes) {
  if (length(genes) == 0) return(NULL)
  long <- DEG_all %>%
    dplyr::mutate(gene = .data[[gene_col]]) %>%
    dplyr::filter(gene %in% genes) %>%
    dplyr::group_by(gene, cell_type) %>%
    dplyr::slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(status = dplyr::case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "Female",
      adj.P.Val < 0.05 & logFC < 0 ~ "Male",
      TRUE                         ~ "NS"
    )) %>%
    dplyr::select(gene, cell_type, status) %>%
    tidyr::complete(gene = genes, cell_type = lv) %>%            # fill all combos; missing = not tested
    dplyr::mutate(status = ifelse(is.na(status), "NotTested", status),
                  cell_type = factor(cell_type, levels = lv))
  m <- long %>%
    tidyr::pivot_wider(names_from = cell_type, values_from = status) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  m[, lv, drop = FALSE]
}

## A) histogram of gene sharing; non-specific gene names on top of bars
col_purple <- "#6A3D9A"   # concordant
col_gold   <- "#C8920A"   # discordant
col_grey   <- "grey75"    # specific
gene_fontsize <- 9        # gene-name font on bar tops (pt)
y_axis_max    <- 60       # y-axis upper limit
share_summary <- gene_class %>%
  dplyr::mutate(concord = dplyr::case_when(
    n_celltype == 1             ~ "Specific",
    n_female == 0 | n_male == 0 ~ "Concordant",
    TRUE                        ~ "Discordant"
  )) %>%
  dplyr::count(n_celltype, concord)
share_summary$concord <- factor(share_summary$concord,
                                levels = c("Specific", "Concordant", "Discordant"))
bar_top <- share_summary %>%
  dplyr::group_by(n_celltype) %>%
  dplyr::summarise(bar_h = sum(n), .groups = "drop")
off <- max(bar_top$bar_h) * 0.02
# non-specific genes -> one richtext block per n_celltype (same colour grouped, <br> line breaks, italic)
label_rich <- gene_class %>%
  dplyr::filter(n_celltype >= 2) %>%
  dplyr::mutate(
    concord = ifelse(grepl("discordant", category), "Discordant", "Concordant"),
    col     = ifelse(concord == "Discordant", col_gold, col_purple),
    tag     = paste0("<span style='color:", col, "'><i>", gene, "</i></span>")
  ) %>%
  dplyr::arrange(n_celltype, concord, gene) %>%
  dplyr::group_by(n_celltype) %>%
  dplyr::summarise(label = paste(tag, collapse = "<br>"), .groups = "drop") %>%
  dplyr::left_join(bar_top, by = "n_celltype") %>%
  dplyr::mutate(y_anchor = bar_h + off)
y_max <- y_axis_max
p_A <- ggplot(share_summary, aes(x = factor(n_celltype), y = n, fill = concord)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            size = 3, color = "black") +
  ggtext::geom_richtext(
    data = label_rich,
    aes(x = factor(n_celltype), y = y_anchor, label = label),
    inherit.aes = FALSE, vjust = 0, hjust = 0.5,
    size = gene_fontsize / ggplot2::.pt, lineheight = 1.0,
    label.color = NA, fill = NA
  ) +
  scale_fill_manual(
    values = c("Specific"   = col_grey,
               "Concordant" = col_purple,
               "Discordant" = col_gold),
    breaks = c("Specific", "Concordant", "Discordant"),
    labels = c("Specific", "Shared (same direction)", "Shared (opposite direction)")
  ) +
  scale_y_continuous(breaks = seq(0, y_axis_max, 20)) +
  coord_cartesian(ylim = c(0, y_max), clip = "off") +
  labs(x = "Number of cell types in which a gene is DEG",
       y = "Number of genes", fill = NULL) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",
        plot.margin = margin(t = 15, r = 12, b = 10, l = 12))
n_label_max <- gene_class %>% dplyr::filter(n_celltype >= 2) %>%
  dplyr::count(n_celltype) %>% dplyr::pull(n)
n_label_max <- if (length(n_label_max) == 0) 1 else max(n_label_max)
fig_h_A <- max(5, n_label_max * 0.16 + 3.5)

## B) heatmap of all shared genes (>= 2 cell types)
shared_info  <- gene_class %>% dplyr::filter(n_celltype >= 2)
shared_genes <- shared_info$gene
if (length(shared_genes) >= 1) {
  mat_shared <- build_status_mat(shared_genes)
  # order genes (columns) by sharing breadth n_celltype, descending
  ord <- order(-shared_info$n_celltype[match(rownames(mat_shared), shared_info$gene)])
  mat_shared   <- mat_shared[ord, , drop = FALSE]
  mat_shared_t <- t(mat_shared)   # rows = cell types, cols = genes
  
  nB    <- ncol(mat_shared_t)
  showB <- nB <= 120
  wB <- if (showB) max(6, nB * 0.16 + 2.5) else min(16, max(8, nB * 0.035 + 3))
  hB <- if (showB) 4.2 else 3.2
  
  htB <- ComplexHeatmap::Heatmap(
    mat_shared_t, name = "Status", col = col_map,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_column_names = showB,
    column_names_gp  = grid::gpar(fontsize = 8, fontface = "italic"),
    column_names_rot = 45,
    row_names_side   = "left", row_names_gp = grid::gpar(fontsize = 10),
    rect_gp = grid::gpar(col = "white", lwd = 1),
    heatmap_legend_param = legend_param_top
  )
}

## C) heatmap of cell-type-specific genes (significant in exactly 1 cell type)
specific_info <- gene_class %>%
  dplyr::filter(n_celltype == 1) %>%
  dplyr::mutate(sig_ct  = celltypes,
                sig_dir = ifelse(n_female > 0, "Female", "Male"))
specific_genes <- specific_info$gene
if (length(specific_genes) >= 1) {
  mat_spec <- build_status_mat(specific_genes)
  # order genes by (significant cell type in lv order, direction, gene); split by cell type
  idx       <- match(rownames(mat_spec), specific_info$gene)
  sig_ct_v  <- specific_info$sig_ct[idx]
  sig_dir_v <- specific_info$sig_dir[idx]
  ordc      <- order(match(sig_ct_v, lv), sig_dir_v, rownames(mat_spec))
  mat_spec  <- mat_spec[ordc, , drop = FALSE]
  
  idx2        <- match(rownames(mat_spec), specific_info$gene)
  col_split_c <- droplevels(factor(specific_info$sig_ct[idx2], levels = lv))
  mat_spec_t  <- t(mat_spec)
  
  nC    <- ncol(mat_spec_t)
  showC <- nC <= 120
  wC <- if (showC) max(6, nC * 0.16 + 2.5) else min(18, max(8, nC * 0.03 + 3))
  hC <- if (showC) 4.5 else 3.5
  
  htC <- ComplexHeatmap::Heatmap(
    mat_spec_t, name = "Status", col = col_map,
    cluster_rows = FALSE, cluster_columns = FALSE,
    column_split    = col_split_c,
    column_title_gp = grid::gpar(fontsize = 11),
    column_gap      = grid::unit(1.5, "mm"),
    show_column_names = showC,
    column_names_gp  = grid::gpar(fontsize = 7, fontface = "italic"),
    column_names_rot = 45,
    row_names_side   = "left", row_names_gp = grid::gpar(fontsize = 10),
    rect_gp = grid::gpar(col = "white", lwd = 0.5),
    show_heatmap_legend = FALSE
  )
}

## combined figure: row 1 = A | B, row 2 = C (B legend on top, C no legend)
if (exists("htB") && exists("htC")) {
  wA <- 6
  gB <- grid::grid.grabExpr(ComplexHeatmap::draw(htB, heatmap_legend_side = "top"))
  gC <- grid::grid.grabExpr(ComplexHeatmap::draw(htC))
  
  row1 <- cowplot::plot_grid(p_A, gB, ncol = 2, rel_widths = c(wA, wB),
                             labels = c("A", "B"), label_size = 16)
  combo <- cowplot::plot_grid(row1, gC, ncol = 1, rel_heights = c(fig_h_A, hC),
                              labels = c("", "C"), label_size = 16)
  
  ggsave("DEG_gene_sharing_combined.pdf", combo,
         width = wA + wB, height = fig_h_A + hC, dpi = 300, limitsize = FALSE)
}

# ==============================================================================
# 8. Pairwise cell-type logFC scatter plots
#    Genes: panel B + panel C genes (shared + cell-type-specific DEGs),
#           excluding XIST, TSIX, JPX, FTX, RLIM, XACT
#    For each cell-type pair, only genes tested in BOTH cell types are kept,
#    and genes not significant in either cell type are dropped (not shown,
#    not used for statistics).
#    Significance = FDR (BH within input X-chr genes per cell type) < 0.05.
#    Stats per pair: correlation (r), R^2 and slope from lm(y ~ x);
#    only r, R^2 and slope are shown on the plot.
#    One page per cell-type pair (21 pages); legend inside each plot.
#    [OUTPUT: DEG_pairwise_logFC_scatter_panelBC_noXCI.pdf]
# ==============================================================================
# ============================== Tunable parameters ==============================
cor_method_pair <- "pearson"   # "pearson" or "spearman"
fdr_cutoff      <- 0.05        # FDR threshold for point colouring / filtering
exclude_genes   <- c("XIST", "TSIX", "JPX", "FTX", "RLIM", "XACT")   # XCI-related genes removed
label_mode      <- "all"       # gene labels: "all" / "both" (sig in both) / "none"
show_fit_line   <- TRUE        # draw the lm(y ~ x) fitted line
fit_line_col    <- "grey30"    # colour of the fitted line
fit_line_alpha  <- 0.5         # transparency of the fitted line (0 = invisible, 1 = opaque)
fit_line_width  <- 0.6         # width of the fitted line
# ---- page size (smaller page -> text looks relatively larger) ----
page_w          <- 4.2         # page width (inches)
page_h          <- 4.2         # page height (inches)
# ---- font / point sizes ----
base_font_size    <- 11        # theme base size (axis titles / ticks)
title_size        <- 12        # plot title size
stat_text_size    <- 3.2       # stats annotation size (mm, ggplot units)
gene_label_size   <- 2.6       # gene label size (mm, ggplot units)
point_size        <- 1.8       # scatter point size
# ---- legend inside the plot ----
legend_pos        <- c(0.98, 0.02)   # (x, y) in npc; bottom-right corner
legend_just       <- c(1, 0)         # anchor = legend's bottom-right
legend_text_size  <- 7.5             # legend item text (pt)
legend_title_size <- 8               # legend title text (pt)
# colours keyed by internal group codes (legend text is built per pair)
pair_cols <- c("both" = "#6A3D9A",   # significant in both cell types
               "x"    = "#E66101",   # significant in x-axis cell type only
               "y"    = "#1B9E77")   # significant in y-axis cell type only
# ===============================================================================

# one row per (gene, cell_type): logFC + FDR significance (most significant kept)
logfc_long <- DEG_all %>%
  dplyr::mutate(gene = .data[[gene_col]]) %>%
  dplyr::filter(!is.na(logFC)) %>%
  dplyr::group_by(gene, cell_type) %>%
  dplyr::slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sig = adj.P.Val < fdr_cutoff) %>%
  dplyr::select(gene, cell_type, logFC, sig)

# scatter for one cell-type pair:
# genes tested in both cell types AND significant in at least one of them
make_pair_plot <- function(genes, ct_x, ct_y, label_mode = "none") {
  dx <- logfc_long %>% dplyr::filter(gene %in% genes, cell_type == ct_x) %>%
    dplyr::select(gene, x = logFC, sig_x = sig)
  dy <- logfc_long %>% dplyr::filter(gene %in% genes, cell_type == ct_y) %>%
    dplyr::select(gene, y = logFC, sig_y = sig)
  d <- dplyr::inner_join(dx, dy, by = "gene") %>%
    dplyr::filter(sig_x | sig_y) %>%                      # drop "not significant in either"
    dplyr::mutate(group = dplyr::case_when(
      sig_x & sig_y ~ "both",
      sig_x         ~ "x",
      TRUE          ~ "y"
    ),
    group = factor(group, levels = names(pair_cols)))
  
  # legend text with actual cell-type names; FDR threshold goes in the legend title
  pair_labels <- c(
    "both" = paste0("Significant in both ", ct_x, " and ", ct_y),
    "x"    = paste0("Significant in ", ct_x, " only"),
    "y"    = paste0("Significant in ", ct_y, " only")
  )
  legend_title <- paste0("FDR < ", fdr_cutoff)
  
  # correlation + linear fit (needs >= 3 genes)
  n <- nrow(d)
  r <- pv <- r2 <- slope <- intercept <- NA_real_
  if (n >= 3) {
    ctest <- cor.test(d$x, d$y, method = cor_method_pair, exact = FALSE)
    r  <- unname(ctest$estimate)
    pv <- ctest$p.value
    fit       <- lm(y ~ x, data = d)
    r2        <- summary(fit)$r.squared
    slope     <- unname(coef(fit)["x"])
    intercept <- unname(coef(fit)["(Intercept)"])
  }
  r_name <- if (cor_method_pair == "pearson") "r" else "rho"
  lab <- if (is.na(r)) "" else
    sprintf("%s = %.2f\nR\u00B2 = %.2f\nslope = %.2f", r_name, r, r2, slope)
  
  # symmetric limits so the y = x line is meaningful
  lim <- if (n > 0) max(abs(c(d$x, d$y)), na.rm = TRUE) * 1.08 else 1
  
  # layer order (bottom -> top): reference lines, fitted line, points, text
  p <- ggplot(d, aes(x, y)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey60", linewidth = 0.3)
  if (show_fit_line && !is.na(slope)) {
    p <- p + geom_abline(slope = slope, intercept = intercept,
                         colour = fit_line_col, alpha = fit_line_alpha,
                         linewidth = fit_line_width)
  }
  p <- p +
    geom_point(aes(colour = group), size = point_size, alpha = 0.85) +
    scale_colour_manual(values = pair_cols, labels = pair_labels,
                        breaks = names(pair_cols), drop = FALSE,
                        name = legend_title) +
    annotate("text", x = -lim, y = lim, hjust = 0, vjust = 1,
             size = stat_text_size, label = lab, lineheight = 0.95) +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(x = paste0(ct_x, " sex logFC"),
         y = paste0(ct_y, " sex logFC")) +
    theme_bw(base_size = base_font_size) +
    theme(panel.grid = element_blank())
  
  # gene labels (top layer)
  d_lab <- switch(label_mode,
                  all  = d,
                  both = dplyr::filter(d, group == "both"),
                  d[0, ])
  if (nrow(d_lab) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = d_lab, aes(label = gene), size = gene_label_size, fontface = "italic",
      max.overlaps = Inf, segment.size = 0.2, min.segment.length = 0,
      box.padding = 0.25, seed = 1
    )
  }
  
  list(plot = p,
       stat = data.frame(x_cell_type = ct_x, y_cell_type = ct_y,
                         n_gene = n, method = cor_method_pair,
                         estimate = r, p_value = pv,
                         r_squared = r2, slope = slope, intercept = intercept,
                         n_sig_both = sum(d$group == "both"),
                         n_sig_x_only = sum(d$group == "x"),
                         n_sig_y_only = sum(d$group == "y")))
}

# one page per cell-type pair (x axis = earlier cell type in lv order)
build_pair_pages <- function(genes, label_mode) {
  plots <- list(); stats <- list()
  for (pr in combn(lv, 2, simplify = FALSE)) {
    res <- make_pair_plot(genes, pr[1], pr[2], label_mode)
    plots[[paste0(pr[1], "_vs_", pr[2])]] <- res$plot +
      labs(title = paste0(pr[1], " vs ", pr[2])) +
      guides(colour = guide_legend(ncol = 1, override.aes = list(size = 2.2, alpha = 1))) +
      theme(
        legend.position        = "inside",
        legend.position.inside = legend_pos,
        legend.justification   = legend_just,
        legend.title           = element_text(size = legend_title_size),
        legend.text            = element_text(size = legend_text_size),
        legend.key.size        = grid::unit(0.35, "cm"),
        legend.spacing.y       = grid::unit(0.05, "cm"),
        legend.margin          = margin(2, 4, 2, 4),
        legend.background      = element_rect(fill = scales::alpha("white", 0.85),
                                              colour = "grey70", linewidth = 0.2),
        plot.title             = element_text(hjust = 0.5, face = "bold", size = title_size),
        plot.margin            = margin(4, 6, 4, 4)
      )
    stats[[length(stats) + 1]] <- res$stat
  }
  list(plots = plots, stats = dplyr::bind_rows(stats))
}

# write a list of ggplots to a multi-page PDF (one plot per page)
save_pages <- function(plots, file, w = page_w, h = page_h) {
  pdf(file, width = w, height = h, onefile = TRUE)
  for (p in plots) print(p)
  dev.off()
}

# ---- gene set: panel B + panel C, XCI genes removed ----
genes_BC <- union(shared_genes, specific_genes)
cat("Excluded XCI genes present in panel B/C:",
    paste(intersect(genes_BC, exclude_genes), collapse = ", "), "\n")
genes_use <- setdiff(genes_BC, exclude_genes)
cat("Genes used:", length(genes_use), "\n")

# ---- run ----
pair_res <- build_pair_pages(genes_use, label_mode)
save_pages(pair_res$plots, "DEG_pairwise_logFC_scatter_panelBC_noXCI.pdf")
print(pair_res$stats)

# ------------------------------------------------------------------------------
# 8b. Lower-triangle scatter-plot matrix of the same 21 pairs
#     Columns = x cell type (Exc..Mic), rows = y cell type (Inh..End);
#     same genes / colours / stats as the single-page plots above.
#     Shared legend in the empty upper-right area.
#     [OUTPUT: DEG_pairwise_logFC_scatter_matrix_panelBC_noXCI.pdf]
# ------------------------------------------------------------------------------
# ============================== Tunable parameters ==============================
mat_label_mode      <- "none"  # gene labels in the matrix: "none" / "both" / "all"
mat_panel_in        <- 2.2     # size of each panel (inches)
mat_stat_text_size  <- 2.6     # r / R^2 / slope annotation size (mm)
mat_point_size      <- 1.2     # point size
mat_axis_text_size  <- 7       # tick label size (pt)
mat_axis_title_size <- 11      # cell-type axis title size (pt; outer panels only)
mat_legend_title    <- 11      # legend title size (pt)
mat_legend_text     <- 10      # legend item size (pt)
# generic legend labels (one legend shared by all panels)
mat_pair_labels <- c("both" = "Significant in both cell types",
                     "x"    = "Significant in column cell type only",
                     "y"    = "Significant in row cell type only")
# ===============================================================================

# temporarily use smaller text / points inside make_pair_plot(), restore afterwards
old_stat_text_size <- stat_text_size; old_point_size <- point_size
stat_text_size <- mat_stat_text_size; point_size <- mat_point_size

row_ct <- lv[-1]               # y axis: Inh .. End (top to bottom)
col_ct <- lv[-length(lv)]      # x axis: Exc .. Mic (left to right)
mat_plots <- list()
for (i in seq_along(row_ct)) {
  for (j in seq_len(i)) {      # lower triangle incl. the sub-diagonal
    ct_x <- col_ct[j]; ct_y <- row_ct[i]
    p <- make_pair_plot(genes_use, ct_x, ct_y, mat_label_mode)$plot
    p <- suppressMessages(
      p +
        scale_colour_manual(values = pair_cols, labels = mat_pair_labels,
                            breaks = names(pair_cols), drop = FALSE,
                            name = paste0("FDR < ", fdr_cutoff)) +
        labs(x = if (i == length(row_ct)) ct_x else NULL,   # x title on bottom row only
             y = if (j == 1)              ct_y else NULL) +  # y title on first column only
        guides(colour = guide_legend(ncol = 1, override.aes = list(size = 3, alpha = 1))) +
        theme(
          axis.text   = element_text(size = mat_axis_text_size),
          axis.title  = element_text(size = mat_axis_title_size, face = "bold"),
          plot.margin = margin(2, 2, 2, 2)
        )
    )
    mat_plots[[length(mat_plots) + 1]] <- p
  }
}
stat_text_size <- old_stat_text_size; point_size <- old_point_size   # restore

# legend area spans the empty upper-right block
mat_plots[[length(mat_plots) + 1]] <- patchwork::guide_area()
# letters A-U = 21 panels (row by row), V = legend, # = empty
mat_design <- "
A##VVV
BC#VVV
DEF###
GHIJ##
KLMNO#
PQRSTU
"
p_matrix <- wrap_plots(mat_plots, design = mat_design) +
  plot_layout(guides = "collect") &
  theme(
    legend.title = element_text(size = mat_legend_title),
    legend.text  = element_text(size = mat_legend_text),
    legend.justification = c(0, 1)
  )
ggsave("DEG_pairwise_logFC_scatter_matrix_panelBC_noXCI.pdf", p_matrix,
       width  = mat_panel_in * length(col_ct),
       height = mat_panel_in * length(row_ct), limitsize = FALSE)

# ==============================================================================
# 9. Lower-triangle heatmaps across cell-type pairs (from step 8)
#    (a) R^2 and (b) slope of lm(logFC_y ~ logFC_x),
#    y = later cell type in lv order
#    Genes: panel B + panel C genes excluding XCI genes; per pair only genes
#    tested in both cell types and FDR-significant in at least one (as in step 8).
#    Colour legend on top (horizontal bar); plot title above the legend.
#    [OUTPUT: DEG_pairwise_R2_heatmap_panelBC_noXCI.pdf,
#             DEG_pairwise_slope_heatmap_panelBC_noXCI.pdf]
# ==============================================================================
# ============================== Tunable parameters ==============================
include_diag    <- FALSE       # TRUE: 7 x 7 with a blank diagonal; FALSE: 6 x 6 lower triangle
r2_digits       <- 2           # decimals for R^2 labels
slope_digits    <- 2           # decimals for slope labels
tile_text_size  <- 3.8         # tile label size (mm, ggplot units)
# shared sequential palette for all three heatmaps
# (viridis family; options: "viridis", "magma", "plasma", "mako", "rocket")
heat_option     <- "viridis"
r2_limits       <- c(0, 1)     # colour range for R^2
slope_limits    <- NULL        # NULL: c(0, max slope); or fix, e.g. c(0, 1.5)
# ---- plot titles ----
title_r2        <- "Sex logFC R\u00B2"
title_slope     <- "Sex logFC slope"
heat_title_size <- 14          # plot title size (pt)
title_gap       <- 20          # space between title and legend (pt); larger = legend lower
# ---- legend on top ----
legend_bar_width   <- 5        # colour bar length (cm)
legend_bar_height  <- 0.45     # colour bar thickness (cm)
legend_title_size  <- 13       # legend title size (pt)
legend_text_size   <- 11.5     # legend tick label size (pt)
# ---- pdf size (title + legend on top -> slightly taller) ----
heat_w          <- 4.6         # pdf width (inches)
heat_h          <- 5.3         # pdf height (inches)
# ===============================================================================

# long table of the 21 pairs (+ optional blank diagonal)
heat_df <- pair_res$stats %>%
  dplyr::select(x_cell_type, y_cell_type, r = estimate, r_squared, slope)

# sanity check: palettes assume non-negative r and slope
if (any(heat_df$r < 0, na.rm = TRUE) || any(heat_df$slope < 0, na.rm = TRUE)) {
  warning("Negative r or slope found; those tiles are squished to the lower colour limit.")
}
if (is.null(slope_limits)) slope_limits <- c(0, max(heat_df$slope, na.rm = TRUE))

if (include_diag) {
  heat_df <- dplyr::bind_rows(
    heat_df,
    data.frame(x_cell_type = lv, y_cell_type = lv,
               r = NA_real_, r_squared = NA_real_, slope = NA_real_)
  )
  x_lv <- lv; y_lv <- rev(lv)
} else {
  x_lv <- lv[-length(lv)]; y_lv <- rev(lv[-1])
}
heat_df <- heat_df %>%
  dplyr::mutate(
    is_diag = x_cell_type == y_cell_type,
    # x = earlier cell type, y = later cell type (top row = earliest)
    x = factor(x_cell_type, levels = x_lv),
    y = factor(y_cell_type, levels = y_lv)
  )

# shared fill scale (viridis, reversed: low = light, high = dark); horizontal bar on top
heat_fill <- function(limits, name) {
  scale_fill_viridis_c(option = heat_option, direction = -1,
                       limits = limits, oob = scales::squish,
                       na.value = "grey92", name = name,
                       guide = guide_colourbar(
                         title.position = "left",
                         title.vjust    = 0.85,
                         barwidth       = grid::unit(legend_bar_width, "cm"),
                         barheight      = grid::unit(legend_bar_height, "cm")
                       ))
}

# build one triangle heatmap; text colour (black/white) chosen from tile luminance
make_tri_heat <- function(df, value_col, fill_scale, digits, title = NULL) {
  df$value <- df[[value_col]]
  df$label <- ifelse(df$is_diag, "",
                     ifelse(is.na(df$value), "NA",
                            formatC(df$value, format = "f", digits = digits)))
  base <- ggplot(df, aes(x = x, y = y)) +
    geom_tile(aes(fill = value), colour = "white", linewidth = 0.8) +
    fill_scale
  # luminance of the rendered fill -> white text on dark tiles
  fills <- ggplot2::layer_data(base, 1)$fill
  rgbm  <- grDevices::col2rgb(fills) / 255
  lum   <- 0.299 * rgbm[1, ] + 0.587 * rgbm[2, ] + 0.114 * rgbm[3, ]
  df$txt_col <- ifelse(lum < 0.5, "white", "black")
  
  base +
    geom_text(data = df, aes(label = label, colour = txt_col),
              size = tile_text_size) +
    scale_colour_identity() +
    coord_equal() +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = heat_title_size,
                                      margin = margin(b = title_gap)),
      plot.title.position = "plot",
      panel.grid       = element_blank(),
      axis.text.x      = element_text(angle = 45, hjust = 1, colour = "black", size = 12),
      axis.text.y      = element_text(colour = "black", size = 12),
      legend.position  = "top",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.title     = element_text(size = legend_title_size),
      legend.text      = element_text(size = legend_text_size),
      legend.margin    = margin(0, 0, 0, 0),
      legend.box.spacing = grid::unit(0.2, "cm"),
      plot.margin      = margin(5, 5, 5, 5)
    )
}

# (a) R^2
p_r2 <- make_tri_heat(heat_df, "r_squared", heat_fill(r2_limits, expression(R^2)), r2_digits,
                      title = title_r2)
print(p_r2)
ggsave("DEG_pairwise_R2_heatmap_panelBC_noXCI.pdf", p_r2,
       width = heat_w, height = heat_h, dpi = 300)

# (b) slope
p_slope <- make_tri_heat(heat_df, "slope", heat_fill(slope_limits, "Slope"), slope_digits,
                         title = title_slope)
print(p_slope)
ggsave("DEG_pairwise_slope_heatmap_panelBC_noXCI.pdf", p_slope,
       width = heat_w, height = heat_h, dpi = 300)

# ==============================================================================
# 10. logFC correlation: adjusted vs non-adjusted case/control
#    (adjusted model computed in-memory only; not written to disk)
#    [OUTPUT: DEG_logFC_correlation_adj_vs_nonadj.pdf]
# ==============================================================================
DEG_all_non_adj <- DEG_all
DEG_all_adj     <- run_deg(adjust_case_control = TRUE)   # in-memory only, no file

identical(DEG_all_adj$gene_cell, DEG_all_non_adj$gene_cell)

x <- DEG_all_adj$logFC
y <- DEG_all_non_adj$logFC
# auto: Pearson if normal, else Spearman
norm_ok <- if (length(x) <= 5000) {
  shapiro.test(x)$p.value > 0.05 && shapiro.test(y)$p.value > 0.05
} else {
  # n>5000 can't run Shapiro; use skewness (|skew|<1 ~ approx symmetric)
  sk <- function(v) mean((v - mean(v))^3) / sd(v)^3
  abs(sk(x)) < 1 && abs(sk(y)) < 1
}
method <- if (norm_ok) "pearson" else "spearman"
ct_cor <- cor.test(x, y, method = method, exact = FALSE)
r_val <- as.numeric(ct_cor$estimate)
label_method <- if (method == "pearson") "Pearson r" else "Spearman rho"
cat(sprintf("Normality passed: %s -> using %s\n", norm_ok, method))
cat(sprintf("%s = %.3f, p = %.3e\n", label_method, r_val, ct_cor$p.value))

lim <- range(c(x, y))
df_cor <- data.frame(x = x, y = y)
p <- ggplot(df_cor, aes(x, y)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.5, size = 1.6, color = "#2c7fb8") +
  geom_smooth(method = "lm", se = TRUE, color = "#d95f02", fill = "#d95f02", alpha = 0.15) +
  annotate("text", x = lim[1], y = lim[2], hjust = 0, vjust = 1, size = 3.8,
           label = sprintf("%s = %.6f\nn = %d", label_method, r_val, length(x))) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "logFC (adjusted case/control)",
       y = "logFC (non-adjusted case/control)",
       title = NULL) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
print(p)
ggsave("DEG_logFC_correlation_adj_vs_nonadj.pdf", p, width = 5.5, height = 5.5)

# ==============================================================================
# 11. Sensitivity analysis for steps 8/9: keep the XCI regulator genes
#     (XIST, TSIX, JPX, FTX, RLIM, XACT) in the gene set.
#     Same functions / parameters as steps 8 and 9; only the gene set differs.
#     Run after steps 8-9 (needs make_pair_plot, build_pair_pages, save_pages,
#     make_tri_heat, heat_fill, genes_BC and the step 9 parameters).
#     All outputs go to the sub-folder <sens_dir>.
#     [OUTPUT: DEG_pairwise_logFC_scatter_panelBC_withXCI.pdf,
#              DEG_pairwise_R2_heatmap_panelBC_withXCI.pdf,
#              DEG_pairwise_slope_heatmap_panelBC_withXCI.pdf]
# ==============================================================================
# ============================== Tunable parameters ==============================
sens_suffix       <- "withXCI"   # output file suffix
sens_dir          <- "sensitivity_withXCI"   # output folder (created if missing); main results stay in the working directory
sens_slope_limits <- NULL        # NULL: c(0, max slope of this run); or fix, e.g. c(0, 1.5)
# ===============================================================================

dir.create(sens_dir, showWarnings = FALSE)
cat("Sensitivity run: XCI genes kept:",
    paste(intersect(genes_BC, exclude_genes), collapse = ", "), "\n")
cat("Genes used:", length(genes_BC), "\n")

# ---- scatter plots (one page per pair) + stats ----
pair_res_sens <- build_pair_pages(genes_BC, label_mode)
save_pages(pair_res_sens$plots,
           file.path(sens_dir, paste0("DEG_pairwise_logFC_scatter_panelBC_", sens_suffix, ".pdf")))
print(pair_res_sens$stats)

# ---- lower-triangle heatmaps (R^2, slope) ----
heat_df_sens <- pair_res_sens$stats %>%
  dplyr::select(x_cell_type, y_cell_type, r = estimate, r_squared, slope)
if (any(heat_df_sens$r < 0, na.rm = TRUE) || any(heat_df_sens$slope < 0, na.rm = TRUE)) {
  warning("Negative r or slope found; those tiles are squished to the lower colour limit.")
}
if (is.null(sens_slope_limits)) sens_slope_limits <- c(0, max(heat_df_sens$slope, na.rm = TRUE))

if (include_diag) {
  heat_df_sens <- dplyr::bind_rows(
    heat_df_sens,
    data.frame(x_cell_type = lv, y_cell_type = lv,
               r = NA_real_, r_squared = NA_real_, slope = NA_real_)
  )
}
heat_df_sens <- heat_df_sens %>%
  dplyr::mutate(
    is_diag = x_cell_type == y_cell_type,
    x = factor(x_cell_type, levels = x_lv),
    y = factor(y_cell_type, levels = y_lv)
  )

p_r2_sens <- make_tri_heat(heat_df_sens, "r_squared", heat_fill(r2_limits, expression(R^2)),
                           r2_digits, title = title_r2)
print(p_r2_sens)
ggsave(file.path(sens_dir, paste0("DEG_pairwise_R2_heatmap_panelBC_", sens_suffix, ".pdf")), p_r2_sens,
       width = heat_w, height = heat_h, dpi = 300)

p_slope_sens <- make_tri_heat(heat_df_sens, "slope", heat_fill(sens_slope_limits, "Slope"),
                              slope_digits, title = title_slope)
print(p_slope_sens)
ggsave(file.path(sens_dir, paste0("DEG_pairwise_slope_heatmap_panelBC_", sens_suffix, ".pdf")), p_slope_sens,
       width = heat_w, height = heat_h, dpi = 300)
