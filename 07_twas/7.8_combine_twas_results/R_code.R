# Combine TWAS results (non-PAR): long table, best eQTL source per gene, FDR-support set, stacked chrX Manhattan overview
rm(list = ls())
library(data.table)
library(tidyverse)
library(ggrepel)
# CONFIG
job_ref_file <- "../7.7_twas/TWAS_job_reference_all.txt"
tss_xlsx     <- "../../03_merge_eQTL_res_X/tss_gene_X_chr.xlsx"
fdr_cut      <- 0.05
# 1. Load TWAS outputs
df_analysis_ref <- fread(job_ref_file)
df_analysis_ref <- filter(df_analysis_ref, PAR_type == 'non_PAR')
res_list <- vector("list", nrow(df_analysis_ref))
for (i in seq_len(nrow(df_analysis_ref))) {
  df_twas_res <- fread(df_analysis_ref$out_file[i])
  df_twas_res <- filter(df_twas_res, !is.na(TWAS.Z))
  
  df_twas_res[, TWAS.fdr := p.adjust(TWAS.P, method = "fdr")]
  df_twas_res[, PAR_type := df_analysis_ref$PAR_type[i]]
  df_twas_res[, trait := df_analysis_ref$trait[i]]
  df_twas_res[, eqtl_sex := df_analysis_ref$sex[i]]
  df_twas_res[, cell_type := df_analysis_ref$cell_type[i]]
  df_twas_res[, gwas_model := df_analysis_ref$gwas_model[i]]
  
  res_list[[i]] <- df_twas_res
}
df_twas_all <- rbindlist(res_list, fill = TRUE)
df_tss <- readxl::read_excel(tss_xlsx)
df_tss <- dplyr::select(df_tss, gene, gene_start)
colnames(df_tss)[1] <- "ID"
df_twas_all <- left_join(df_twas_all, df_tss)
writexl::write_xlsx(df_twas_all, "df_twas_all.xlsx")
saveRDS(df_twas_all, "df_twas_all.rds")
# 2. Wide table across eqtl_sex; best source = smallest TWAS.P; eqtl_support = sources with FDR < cutoff
df_twas_all <- dplyr::select(df_twas_all, cell_type, eqtl_sex, trait, gwas_model, Gene = ID,
                             CHR, gene_start, TWAS.Z, TWAS.P, TWAS.fdr)
df_wide <- df_twas_all %>%
  group_by(cell_type, trait, gwas_model, Gene, eqtl_sex) %>%
  slice_min(TWAS.P, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  pivot_wider(
    names_from = eqtl_sex,
    values_from = c(TWAS.Z, TWAS.P, TWAS.fdr),
    names_sep = "_"
  )
df_wide <- df_wide %>%
  mutate(
    best_eqtl_source = c("both", "female", "male")[
      max.col(-cbind(
        ifelse(is.na(TWAS.P_both),   Inf, TWAS.P_both),
        ifelse(is.na(TWAS.P_female), Inf, TWAS.P_female),
        ifelse(is.na(TWAS.P_male),   Inf, TWAS.P_male)
      ), ties.method = "first")
    ],
    TWAS.Z_best = case_when(
      best_eqtl_source == "both"   ~ TWAS.Z_both,
      best_eqtl_source == "female" ~ TWAS.Z_female,
      best_eqtl_source == "male"   ~ TWAS.Z_male
    ),
    TWAS.P_best = case_when(
      best_eqtl_source == "both"   ~ TWAS.P_both,
      best_eqtl_source == "female" ~ TWAS.P_female,
      best_eqtl_source == "male"   ~ TWAS.P_male
    ),
    TWAS.fdr_best = case_when(
      best_eqtl_source == "both"   ~ TWAS.fdr_both,
      best_eqtl_source == "female" ~ TWAS.fdr_female,
      best_eqtl_source == "male"   ~ TWAS.fdr_male
    )
  ) %>%
  relocate(TWAS.Z_best, TWAS.P_best, TWAS.fdr_best, best_eqtl_source,
           .after = gene_start)
df_wide <- df_wide %>%
  rowwise() %>%
  mutate(
    eqtl_support = paste(
      c("both", "female", "male")[
        c(
          !is.na(TWAS.fdr_both)   & TWAS.fdr_both   < fdr_cut,
          !is.na(TWAS.fdr_female) & TWAS.fdr_female < fdr_cut,
          !is.na(TWAS.fdr_male)   & TWAS.fdr_male   < fdr_cut
        )
      ],
      collapse = ";"
    ),
    eqtl_support = ifelse(eqtl_support == "", NA_character_, eqtl_support)
  ) %>%
  ungroup() %>%
  relocate(eqtl_support, .after = gene_start)
saveRDS(df_wide, "twas_res_distinct.rds")
writexl::write_xlsx(df_wide, "twas_res_distinct.xlsx")
write.csv(df_wide, "twas_res_distinct.csv", row.names = FALSE)
# 3. chrX TWAS overview: one row per disease, cell types overlaid; colour = cell type, diamond = FDR < 0.05, labels for FDR hits
df_raw <- readRDS("twas_res_distinct.rds")
traits_order <- c("AD", "PD", "MSA", "LBD")
cell_types   <- c("Ast", "End", "Exc", "Inh", "Mic", "OPC", "Oli")
gene_col     <- "Gene"
pt_size   <- 0.7
pt_alpha  <- 0.43
sig_size  <- 1.5
lab_size  <- 2.5
fdr_cut   <- 0.05
ymax      <- 6   # -log10(P) cap; NA = none
x_break   <- 20   # Mb
celltype_col <- c(
  Ast = "#1B9E77", End = "#D95F02", Exc = "#7570B3", Inh = "#E7298A",
  Mic = "#66A61E", OPC = "#E6AB02", Oli = "#A6761D"
)
shape_ns  <- 21
shape_sig <- 23
label_model <- FALSE   # append best model to gene labels
chrX_len <- 156040895   # hg38
sig_levels <- c("FDR < 0.05", "FDR \u2265 0.05")
cap_y      <- function(x) if (is.na(ymax)) x else pmin(x, ymax)
df <- df_raw %>%
  filter(trait %in% traits_order,
         cell_type %in% cell_types,
         !is.na(gene_start), !is.na(TWAS.P_best)) %>%
  mutate(
    logp      = -log10(TWAS.P_best),
    logp_c    = cap_y(logp),
    pos_mb    = gene_start / 1e6,
    is_sig    = !is.na(TWAS.fdr_best) & TWAS.fdr_best < fdr_cut,
    cell_type = factor(cell_type, levels = cell_types),
    trait     = factor(trait,     levels = traits_order),
    sig_f     = factor(ifelse(is_sig, sig_levels[1], sig_levels[2]),
                       levels = sig_levels)
  )
y_top <- if (is.na(ymax)) ceiling(max(df$logp_c, na.rm = TRUE)) else ymax
# label each significant gene once per disease (strongest cell type)
ord   <- order(as.integer(df$cell_type), decreasing = TRUE)
df    <- df[ord, , drop = FALSE]
d_ns  <- filter(df, !is_sig)
d_sg  <- filter(df,  is_sig)
d_lab <- d_sg[order(-d_sg$logp), , drop = FALSE]
d_lab <- d_lab[!duplicated(d_lab[c("trait", gene_col)]), , drop = FALSE]
d_lab$.lab <- if (label_model) {
  paste0(d_lab[[gene_col]], " (", as.character(d_lab$gwas_model), ")")
} else {
  d_lab[[gene_col]]
}
# invisible layers driving the legends
ct_key  <- data.frame(pos_mb = 0, logp_c = 0,
                      cell_type = factor(cell_types, levels = cell_types))
sig_key <- data.frame(pos_mb = 0, logp_c = 0,
                      sig_f = factor(sig_levels, levels = sig_levels))
p <- ggplot(mapping = aes(pos_mb, logp_c)) +
  geom_point(data = d_ns, aes(fill = cell_type),
             shape = shape_ns, colour = "grey30", size = pt_size,
             stroke = 0.12, alpha = pt_alpha, show.legend = FALSE) +
  geom_point(data = d_sg, aes(fill = cell_type),
             shape = shape_sig, colour = "black", size = sig_size,
             stroke = 0.55, alpha = pt_alpha, show.legend = FALSE) +
  geom_point(data = ct_key,  aes(fill = cell_type),
             shape = 21, size = 0, stroke = 0, na.rm = TRUE) +
  geom_point(data = sig_key, aes(shape = sig_f),
             size = 0, stroke = 0, na.rm = TRUE) +
  geom_label_repel(data = d_lab, aes(label = .lab),
                   size = lab_size, colour = "grey15",
                   fill = "grey95", alpha = 0.7,
                   label.size = 0,
                   label.padding = unit(0.12, "lines"),
                   label.r = unit(0.1, "lines"),
                   segment.colour = "grey60", segment.size = 0.3,
                   min.segment.length = 0, max.overlaps = Inf,
                   box.padding = 0.3, point.padding = 0.2, seed = 1) +
  scale_fill_manual("Cell type", values = celltype_col,
                    limits = cell_types, drop = FALSE) +
  scale_shape_manual("Significance",
                     values = setNames(c(shape_sig, shape_ns), sig_levels),
                     limits = sig_levels, drop = FALSE) +
  scale_x_continuous(limits = c(0, chrX_len / 1e6),
                     breaks = seq(0, 160, by = x_break),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(0, y_top),
                     breaks = pretty(c(0, y_top), n = 4),
                     expand = expansion(mult = c(0.02, 0.12))) +
  facet_grid(trait ~ ., switch = "y") +
  labs(x = "ChrX position (Mb)", y = expression(-log[10](italic(P)))) +
  guides(
    fill  = guide_legend("Cell type", order = 1, override.aes = list(
      shape = 21, size = 2.8, colour = "grey25", stroke = 0.2, alpha = 1)),
    shape = guide_legend("Significance", order = 2, override.aes = list(
      fill   = c("grey75", "grey75"),
      colour = c("black",  "grey45"),
      size   = c(2.9, 1.7),
      stroke = c(0.6, 0.2),
      alpha  = 1))
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(linewidth = 0.2, colour = "grey92"),
    panel.spacing.y   = unit(0.5, "lines"),
    strip.background   = element_rect(fill = "grey95", colour = NA),
    strip.placement    = "outside",
    strip.text.y.left  = element_text(angle = 0, face = "bold", size = 12),
    legend.position    = "right",
    legend.key.size    = unit(0.9, "lines"),
    legend.title       = element_text(face = "bold"),
    axis.title         = element_text(face = "bold")
  )
# cairo_pdf renders the unicode >= in the legend
pdf_dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
W <- 9.5; H <- 6
ggsave("chrX_twas_overview.png", p, width = W, height = H, dpi = 300)
ggsave("chrX_twas_overview.pdf", p, width = W, height = H, device = pdf_dev)
message(sprintf("Done: chrX_twas_overview.png / .pdf  (%g x %g in)", W, H))