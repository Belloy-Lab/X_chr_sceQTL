# GeneCards PD relevance score by gene group per cell type; two versions: raw and exclXCIreg (drops XCI regulator genes)
rm(list = ls())
library(tidyverse)
library(data.table)
library(HGNChelper)
library(openxlsx)
library(ggrepel)
# CONFIG
# --- disease tag (used in axis label and output file names) ---
disease_tag <- "PD"
# --- eQTL input genes (step 3) ---
eqtl_x_rds   <- "../03_merge_eQTL_res_X/df_eqtl_cis_all.rds"
eqtl_par_rds <- "../03_merge_eQTL_res_PAR/df_eqtl_cis_all.rds"
# --- reported eXCI genes (shared with 06.4_venn_plot) ---
exci_xlsx <- "../06.4_venn_plot/df_eXCI_gene_reported.xlsx"
# --- female-biased signal gene set (venn overlap: SEX HET FDR & DEG) ---
female_hetFDR_DEG_csv <- "../06.4_venn_plot/female_biased_female_up_adjP_overlap_genes.csv"
# --- GeneCards PD relevance table (NOT included in this repository) ---
# GeneCards data are copyrighted by LifeMap Sciences and cannot be redistributed,
# so this file must be downloaded by each user:
#   1. Log in to GeneCards (https://www.genecards.org/; free for academic users).
#   2. Search for: Parkinson's disease
#   3. On the results page, export the results table as CSV.
#   4. Save the file in this folder under the name below (or edit genecards_csv).
# Required columns: "Symbol" and "Relevance Score".
# The paper used GeneCards search results downloaded in June 2026. Scores change between
# GeneCards releases, so a newer download may give slightly different results.
genecards_csv <- "GeneCards - GeneCards - Search results for Parkinson's disease.csv"
if (!file.exists(genecards_csv))
  stop("GeneCards file not found: '", genecards_csv,
       "'. See the CONFIG comments above for how to download it.")
# --- parameters ---
score_threshold <- 1                       # eXCI escape-score threshold for "known eXCI"
bar2_flag       <- "female_hetFDR_DEG"      # evidence source for the signal group (SEX HET FDR & DEG)
selected_cell_types <- c("Exc","Inh","Ast","Oli","OPC","Mic")   # cell types to display, in this order (NULL = all)
ref_lines <- c(15, 50)                      # reference dashed lines on the y axis
# X-inactivation master-regulator genes (XIST and its cis regulators / antisense),
# dropped only in the exclXCIreg version.
genes_to_drop <- c("XIST", "TSIX", "JPX", "FTX", "RLIM", "XACT")
# the two versions to run
run_modes <- list(
  list(suffix = "raw",        drop_xci_regulators = FALSE),
  list(suffix = "exclXCIreg", drop_xci_regulators = TRUE)
)
# font sizes
base_font          <- 16   # axis / legend (pt)
n_font_pt          <- 11   # n= labels (pt)
sig_font_pt        <- 11   # significance stars (pt)
ns_font_pt         <- 7    # "ns" label (pt); smaller than stars
ns_vjust           <- 0    # vertical nudge for "ns" (negative = higher, e.g. -0.3)
num_font_pt        <- 9    # gene-name labels (pt)
effect_lab_font_pt <- 8    # delta effect-size labels (pt)
# 1. Shared inputs (loaded once; version-independent)
# eQTL input genes (full, before any dropping)
df_eqtl_all_X   <- read_rds(eqtl_x_rds)
df_eqtl_all_PAR <- read_rds(eqtl_par_rds)
df_input_full <- rbind(distinct(df_eqtl_all_X,   cell_type, gene),
                       distinct(df_eqtl_all_PAR, cell_type, gene))
# known eXCI genes
df_exci <- readxl::read_excel(exci_xlsx)
df_exci <- filter(df_exci, escape_score >= score_threshold)
df_exci <- dplyr::select(df_exci, gene = gene_name)
df_exci$known <- "Yes"
df_exci <- distinct(df_exci, gene, .keep_all = T)
# female-biased signal genes (SEX HET FDR & DEG)
df_female_hetFDR_DEG <- fread(female_hetFDR_DEG_csv)
# GeneCards PD relevance scores (base table; zero-fill for missing input genes done per version)
df_association_base <- fread(genecards_csv)
df_association_base$associationScore <- df_association_base$`Relevance Score`
df_hgnc_check <- checkGeneSymbols(df_association_base$Symbol)
df_association_base$gene <- df_hgnc_check$Suggested.Symbol
df_association_base <- dplyr::arrange(df_association_base, desc(associationScore))
df_association_base <- distinct(df_association_base, gene, .keep_all = T)
df_association_base <- dplyr::select(df_association_base, gene, associationScore)
# 2. Statistics helpers
## Effect size: Cliff's delta (non-parametric, -1..1, tie-corrected).
## Positive = group a tends to score higher than group b. Uses mid-ranks (U statistic).
cliffs_delta <- function(a, b) {
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  n_a <- length(a); n_b <- length(b)
  if (n_a < 1 || n_b < 1) return(NA_real_)
  r   <- rank(c(a, b))                       # mid-rank (ties averaged)
  R_a <- sum(r[seq_len(n_a)])
  U_a <- R_a - n_a * (n_a + 1) / 2           # Mann-Whitney U for group a
  (2 * U_a) / (n_a * n_b) - 1                # = P(a>b) - P(a<b)
}
## p-value + test used + effect size (Cliff's delta)
run_test <- function(scores_a, scores_b) {
  scores_a <- scores_a[!is.na(scores_a)]
  scores_b <- scores_b[!is.na(scores_b)]
  if (length(scores_a) < 3 | length(scores_b) < 3)
    return(list(p_value = NA_real_, test_used = "insufficient data",
                effect = NA_real_))
  sw_a <- if (length(scores_a) <= 5000) shapiro.test(scores_a)$p.value else NA
  sw_b <- if (length(scores_b) <= 5000) shapiro.test(scores_b)$p.value else NA
  both_normal <- !is.na(sw_a) & !is.na(sw_b) & sw_a > 0.05 & sw_b > 0.05
  if (both_normal) {
    res <- t.test(scores_a, scores_b); test_used <- "t-test"
  } else {
    res <- wilcox.test(scores_a, scores_b, exact = FALSE); test_used <- "Mann-Whitney"
  }
  list(p_value   = res$p.value,
       test_used = test_used,
       effect    = cliffs_delta(scores_a, scores_b))
}
## star label: * < 0.05, ** < 0.01, *** < 0.001, ns otherwise
make_label <- function(p) dplyr::case_when(
  is.na(p)  ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"
)
## numeric label: print the (adjusted) p/q value; prefix set at call site
make_label_num <- function(p, prefix = "p=") dplyr::case_when(
  is.na(p)   ~ "",
  p < 1e-300 ~ paste0(prefix, "<1e-300"),                          # guard against 0 / underflow
  TRUE       ~ paste0(prefix, formatC(p, format = "e", digits = 1))
)
## effect-size label: delta=0.42 (\u03B4 is the Greek letter; needs cairo_pdf to render)
make_label_effect <- function(d) dplyr::case_when(
  is.na(d) ~ "",
  TRUE     ~ paste0("\u03B4=", formatC(d, format = "f", digits = 2))
)
save_excel <- function(summary_df, stat_df, filename) {
  wb <- createWorkbook()
  addWorksheet(wb, "summary")
  addWorksheet(wb, "statistics")
  writeData(wb, "summary",    summary_df)
  writeData(wb, "statistics", stat_df)
  saveWorkbook(wb, filename, overwrite = TRUE)
}
## helper: add a flag column when (cell_type, gene) matches a reference set
add_match_flag <- function(target_df, ref_df,
                           key_cols = c("cell_type", "gene"),
                           new_col = "in_female_biased",
                           flag = "yes", default = NA_character_) {
  was_dt <- data.table::is.data.table(target_df)
  t_keys <- do.call(paste, c(as.data.frame(target_df)[key_cols], sep = "\r"))
  r_keys <- do.call(paste, c(as.data.frame(ref_df)[key_cols],   sep = "\r"))
  target_df <- as.data.frame(target_df)
  target_df[[new_col]] <- ifelse(t_keys %in% r_keys, flag, default)
  if (was_dt) data.table::as.data.table(target_df) else target_df
}
# 3. Plot constants + font conversions + plotting function (shared)
group_levels <- c("Female upregulated DEGs and/or Female biased eQTL genes",
                  "Non-overlapping known eXCI genes",
                  "Genes in neither group")
group_colors <- c(
  "Female upregulated DEGs and/or Female biased eQTL genes"  = "#E64B35",
  "Non-overlapping known eXCI genes" = "#7B4F9E",
  "Genes in neither group"            = "#B0B0B0"
)
n_font          <- n_font_pt          / ggplot2::.pt   # geom_text uses mm; convert pt -> mm
sig_font        <- sig_font_pt        / ggplot2::.pt
ns_font         <- ns_font_pt         / ggplot2::.pt   # "ns" label
num_font        <- num_font_pt        / ggplot2::.pt   # also used for the ggrepel gene labels
effect_lab_font <- effect_lab_font_pt / ggplot2::.pt
## Plotting function. stat_tbl needs sig_label (stars / "ns", lower line) and
## effect_label (delta, upper line). Bracket geometry set by y_base_mult / tier_gap / eff_gap.
build_plot <- function(plot_df, stat_tbl,
                       label_font    = sig_font,         # stars font
                       ns_label_font = ns_font,          # "ns" font (smaller)
                       effect_font   = effect_lab_font,  # delta font (smaller)
                       n_top         = 3,                # top N genes labelled per cell type
                       y_base_mult   = 1.05,             # overall bracket height
                       tier_gap      = 0.10,             # gap between the two brackets
                       eff_gap       = 0.07) {           # gap between delta and stars
  ## ---- base plot (dashed reference lines under the violins) ----
  p_base <- ggplot(plot_df, aes(x = cell_type, y = associationScore, fill = group)) +
    geom_hline(yintercept = ref_lines, linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_violin(position = position_dodge(0.85), width = 0.8, alpha = 0.8, scale = "width") +
    geom_boxplot(position = position_dodge(0.85), width = 0.12, outlier.shape = NA, alpha = 0.5) +
    scale_fill_manual(values = group_colors) +
    scale_y_continuous(
      breaks = sort(unique(c(
        scales::extended_breaks()(range(plot_df$associationScore, na.rm = TRUE)),
        ref_lines
      ))),
      expand = expansion(mult = c(0.05, 0.08))
    ) +
    theme_classic(base_size = base_font)
  ## ---- read the real dodged x-centres from the boxplot layer ----
  pb  <- ggplot_build(p_base)
  box <- pb$data[[which(sapply(pb$data, function(d) "ymin_final" %in% names(d)))[1]]]
  box$cell_type_num <- round(box$x)
  box$group_name    <- names(group_colors)[match(box$fill, group_colors)]
  ct_levels <- levels(factor(plot_df$cell_type))
  ct_map    <- tibble::tibble(cell_type = ct_levels, cell_type_num = seq_along(ct_levels))
  xc <- box %>% dplyr::transmute(cell_type_num, group = group_name, x_center = x)
  ## ---- n labels ----
  n_counts <- plot_df %>% dplyr::count(cell_type, group, name = "n_genes")
  y_pos_global <- max(plot_df$associationScore, na.rm = TRUE) * 1.02
  n_lab <- xc %>%
    dplyr::left_join(ct_map,   by = "cell_type_num") %>%
    dplyr::left_join(n_counts, by = c("cell_type", "group")) %>%
    dplyr::mutate(y_pos = y_pos_global)
  ## ---- significance brackets ----
  tier_map <- c("Non-overlapping known eXCI genes" = 1,
                "Female upregulated DEGs and/or Female biased eQTL genes"  = 2)
  y_base_global <- max(plot_df$associationScore, na.rm = TRUE) * y_base_mult
  none_x <- xc %>%
    dplyr::filter(group == "Genes in neither group") %>%
    dplyr::transmute(cell_type_num, x_none = x_center)
  signif_df <- stat_tbl %>%
    dplyr::mutate(grp = sub(" vs Genes in neither group$", "", comparison)) %>%
    dplyr::left_join(ct_map, by = "cell_type") %>%
    dplyr::left_join(xc, by = c("cell_type_num", "grp" = "group")) %>%
    dplyr::rename(x_test = x_center) %>%
    dplyr::left_join(none_x, by = "cell_type_num") %>%
    dplyr::mutate(
      x_left  = pmin(x_test, x_none),
      x_right = pmax(x_test, x_none),
      x_mid   = (x_left + x_right) / 2,
      tier    = tier_map[grp],
      y_bar   = y_base_global * (1 + tier_gap * tier),   # bracket bar
      y_tip   = y_bar - y_base_global * 0.03,            # bracket tips
      y_sig   = y_bar,                                   # stars (on the bracket)
      y_eff   = y_bar + y_base_global * eff_gap          # delta (slightly higher)
    )
  ## split labels: significant stars vs "ns" (drawn with different font sizes)
  sig_df <- dplyr::filter(signif_df, sig_label != "ns")
  ns_df  <- dplyr::filter(signif_df, sig_label == "ns")
  ## top-N genes of the first group: mark position + gene name
  biased_group <- group_levels[1]
  top_lab <- plot_df %>%
    dplyr::filter(group == biased_group, !is.na(associationScore)) %>%
    dplyr::mutate(cell_type = as.character(cell_type)) %>%
    dplyr::group_by(cell_type) %>%
    dplyr::slice_max(order_by = associationScore, n = n_top, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(ct_map, by = "cell_type") %>%
    dplyr::left_join(
      dplyr::filter(xc, group == biased_group) %>%
        dplyr::select(cell_type_num, x_center),
      by = "cell_type_num"
    )
  ## ---- overlay layers ----
  p_out <- p_base +
    geom_text(data = n_lab,
              aes(x = x_center, y = y_pos, label = paste0("n=", n_genes)),
              vjust = 0, size = n_font, inherit.aes = FALSE) +
    geom_segment(data = signif_df,
                 aes(x = x_left,  xend = x_right, y = y_bar, yend = y_bar),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = signif_df,
                 aes(x = x_left,  xend = x_left,  y = y_tip, yend = y_bar),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = signif_df,
                 aes(x = x_right, xend = x_right, y = y_tip, yend = y_bar),
                 inherit.aes = FALSE, linewidth = 0.4) +
    ## lower line: stars (significant comparisons, normal size)
    geom_text(data = sig_df,
              aes(x = x_mid, y = y_sig, label = sig_label),
              inherit.aes = FALSE, size = label_font, hjust = 0.5, vjust = 0) +
    ## lower line: "ns" (non-significant comparisons, smaller size)
    geom_text(data = ns_df,
              aes(x = x_mid, y = y_sig, label = sig_label),
              inherit.aes = FALSE, size = ns_label_font, hjust = 0.5, vjust = ns_vjust) +
    ## upper line: delta effect size (smaller font)
    geom_text(data = signif_df,
              aes(x = x_mid, y = y_eff, label = effect_label),
              inherit.aes = FALSE, size = effect_font, hjust = 0.5, vjust = 0) +
    ## top-N position points
    geom_point(data = top_lab,
               aes(x = x_center, y = associationScore),
               inherit.aes = FALSE,
               shape = 21, fill = "white", color = "black",
               size = 1.6, stroke = 0.5)
  ## top-N gene labels: grey rounded label boxes
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_out <- p_out +
      ggrepel::geom_label_repel(
        data = top_lab,
        aes(x = x_center, y = associationScore, label = gene),
        inherit.aes = FALSE,
        size = num_font, fontface = "italic",
        colour = "grey15",
        fill = "grey90", alpha = 0.7,
        label.size = 0.3,
        label.padding = unit(0.12, "lines"),
        label.r = unit(0.1, "lines"),
        segment.colour = "grey60", segment.size = 0.3,
        min.segment.length = 0, max.overlaps = Inf,
        nudge_x = -0.22, direction = "y",     # nudge labels left
        box.padding = 0.3, point.padding = 0.2, seed = 1
      )
  } else {
    stop("ggrepel is required for the label boxes; install.packages('ggrepel')")
  }
  p_out +
    labs(x = "Cell type",
         y = sprintf("GeneCards %s Relevance Score", disease_tag),
         fill = NULL) +
    theme(
      legend.position = "top",
      axis.title      = element_text(size = base_font),
      axis.text       = element_text(size = base_font),
      axis.text.x     = element_text(size = base_font, angle = 45, hjust = 1),
      legend.text     = element_text(size = base_font),
      legend.title    = element_text(size = base_font)
    )
}
## delta is a unicode character; base pdf() cannot encode it -> use cairo_pdf (fallback pdf)
pdf_dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
# 4. Run both versions (raw / exclXCIreg): star plot + stats xlsx
for (m in run_modes) {
  suffix <- paste0(disease_tag, "_", m$suffix)   # e.g. PD_raw / PD_exclXCIreg
  message("\n=== Version: ", suffix, " ===")
  
  # input genes (optionally drop X-inactivation regulators)
  df_input <- df_input_full
  if (m$drop_xci_regulators) {
    n_before <- nrow(df_input)
    df_input <- dplyr::filter(df_input,
                              !toupper(as.character(gene)) %in% toupper(genes_to_drop))
    message("Dropped ", n_before - nrow(df_input),
            " row(s) for gene(s): ", paste(genes_to_drop, collapse = ", "))
  }
  
  # GeneCards scores restricted to this input set (genes with no hit get 0)
  df_assoc_not_exist <- anti_join(
    data.frame(gene = unique(df_input$gene)),
    df_association_base, by = "gene"
  )
  df_assoc_not_exist$associationScore <- 0
  df_association <- rbind(df_association_base, df_assoc_not_exist)
  
  # evidence flags: female-biased signal set + known eXCI
  df_flagged <- add_match_flag(df_input, df_female_hetFDR_DEG, new_col = "female_hetFDR_DEG")
  df_flagged <- left_join(df_flagged, df_exci, by = "gene")
  
  # three-group plot data
  base_df <- df_flagged %>% dplyr::left_join(df_association, by = "gene")
  plot_df_3grp <- bind_rows(
    # 1. signal group (source controlled by bar2_flag)
    base_df %>%
      dplyr::filter(!is.na(.data[[bar2_flag]]) & .data[[bar2_flag]] == "yes") %>%
      dplyr::mutate(group = "Female upregulated DEGs and/or Female biased eQTL genes"),
    # 2. known eXCI genes not already in the signal group
    base_df %>%
      dplyr::filter(!is.na(known) & known == "Yes" &
                      (is.na(.data[[bar2_flag]]) | .data[[bar2_flag]]  != "yes")) %>%
      dplyr::mutate(group = "Non-overlapping known eXCI genes"),
    # 3. genes in neither group (not signal, not known eXCI)
    base_df %>%
      dplyr::filter((is.na(.data[[bar2_flag]]) | .data[[bar2_flag]]  != "yes") &
                      (is.na(known)            | known            != "Yes")) %>%
      dplyr::mutate(group = "Genes in neither group")
  ) %>%
    dplyr::mutate(group = factor(group, levels = group_levels))
  
  # select / order cell types (NULL = all)
  if (!is.null(selected_cell_types)) {
    missing_ct <- setdiff(selected_cell_types, unique(as.character(plot_df_3grp$cell_type)))
    if (length(missing_ct) > 0)
      warning("cell type name(s) not matched, check spelling: ",
              paste(missing_ct, collapse = ", "))
    plot_df_3grp <- plot_df_3grp %>%
      dplyr::filter(cell_type %in% selected_cell_types) %>%
      dplyr::mutate(cell_type = factor(cell_type, levels = selected_cell_types))
  } else {
    plot_df_3grp <- plot_df_3grp %>%
      dplyr::mutate(cell_type = factor(cell_type))
  }
  
  # statistics: each signal group vs "Genes in neither group", then global BH-FDR
  stat_3grp_raw <- plot_df_3grp %>%
    dplyr::group_by(cell_type) %>%
    dplyr::group_modify(~ {
      s <- split(.x$associationScore, .x$group)
      comparisons <- list(
        c("Female upregulated DEGs and/or Female biased eQTL genes",  "Genes in neither group"),
        c("Non-overlapping known eXCI genes", "Genes in neither group")
      )
      purrr::map_dfr(comparisons, function(pair) {
        r <- run_test(s[[pair[1]]], s[[pair[2]]])
        tibble(comparison = paste(pair[1], "vs", pair[2]),
               p_value    = r$p_value,
               test_used  = r$test_used,
               effect     = r$effect)
      })
    }) %>%
    dplyr::ungroup()
  
  summary_3grp <- plot_df_3grp %>%
    dplyr::group_by(cell_type, group) %>%
    dplyr::summarise(
      mean_score   = mean(associationScore, na.rm = TRUE),
      median_score = median(associationScore, na.rm = TRUE),
      sd_score     = sd(associationScore, na.rm = TRUE),
      n_genes      = n(), .groups = "drop"
    )
  
  stat_tbl <- stat_3grp_raw
  stat_tbl$p_value_adj  <- p.adjust(stat_tbl$p_value, method = "BH")
  stat_tbl$p_label_star <- make_label(stat_tbl$p_value_adj)
  stat_tbl$p_label_num  <- make_label_num(stat_tbl$p_value_adj, "q=")
  stat_tbl$effect_label <- make_label_effect(stat_tbl$effect)
  stat_tbl$adjustment   <- "BH-FDR across cell types"
  stat_tbl <- dplyr::select(stat_tbl, cell_type, comparison, test_used,
                            p_value, p_value_adj, effect, effect_label,
                            p_label_star, p_label_num, adjustment)
  
  # ---- star plot ----
  p_star <- build_plot(plot_df_3grp,
                       dplyr::mutate(stat_tbl, sig_label = p_label_star),
                       label_font = sig_font)
  ggsave(sprintf("female_3grp_violin_ct_specific_fdr_global_star_%s.pdf", suffix),
         plot = p_star, width = 15, height = 4.5, dpi = 300, device = pdf_dev)
  
  # ---- statistics workbook ----
  save_excel(summary_3grp, stat_tbl,
             sprintf("female_3grp_violin_ct_specific_fdr_global_%s.xlsx", suffix))
  
}