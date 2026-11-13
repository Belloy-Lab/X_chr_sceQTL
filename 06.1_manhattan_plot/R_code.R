# Manhattan plot of significant chrX (PAR + non-PAR) cis-eQTLs; top independent eGenes labelled
rm(list = ls())
library(tidyverse)
library(ggrepel)
# CONFIG
sig_par_dir <- "../04_significant_eQTL_PAR"
sig_x_dir   <- "../04_significant_eQTL_X"
indep_x_rds   <- "../05_independent_signal_X/df_eqtl_wide_after_second_clump.rds"
indep_par_rds <- "../05_independent_signal_PAR/df_eqtl_wide_after_second_clump.rds"

# hg38
chrX_len   <- 156040895
par1_start <- 10001
par1_end   <- 2781479
par2_start <- 155701383
par2_end   <- 156030895

# Gene labelling
#   n_label    : number of eGenes to label (most significant first)
#   min_gap_kb : genes within this distance of an already-selected (stronger) gene are not
#                labelled, so each independent locus gets one label (1 Mb = LD clumping window;
#                set to 0 to label by gene only)
n_label    <- 20
min_gap_kb <- 1000

cell_type_colors <- c(
  Ast = "#1b9e77",
  End = "#d95f02",
  Exc = "#7570b3",
  Inh = "#e7298a",
  Mic = "#66a61e",
  OPC = "#e6ab02",
  Oli = "#a6761d"
)
# 1. Significant eQTLs (2-step FDR), PAR + X
PAR_rds_path <- dir(sig_par_dir, pattern = "eqtl_res_PAR", full.names = TRUE)
df_list <- list()
for (i in seq_along(PAR_rds_path)) {
  cat("Reading:", PAR_rds_path[i], "\n")
  df <- readRDS(PAR_rds_path[i])
  df <- df %>%
    select("#CHROM", POS, rsid, cell_type, gene, par_type, setting,  BETA, LOG10_P, significant_by_2step_FDR) %>%
    filter(significant_by_2step_FDR == "Yes")
  df_list[[i]] <- df
}
df_PAR_pass_two_step_fdr <- bind_rows(df_list)

X_rds_path <- dir(sig_x_dir, pattern = "eqtl_res_X", full.names = TRUE)
df_list <- list()
for (i in seq_along(X_rds_path)) {
  cat("Reading:", X_rds_path[i], "\n")
  df <- readRDS(X_rds_path[i])
  df <- df %>%
    select("#CHROM", POS, rsid, cell_type, gene, par_type, setting, BETA, LOG10_P, significant_by_2step_FDR) %>%
    filter(significant_by_2step_FDR == "Yes")
  df_list[[i]] <- df
}
df_X_pass_two_step_fdr <- bind_rows(df_list)

df_X_PAR_pass_two_step_fdr <- rbind(df_PAR_pass_two_step_fdr,
                                    df_X_pass_two_step_fdr)
saveRDS(df_X_PAR_pass_two_step_fdr,"df_X_PAR_pass_two_step_fdr.rds")
# 2. Flag independent lead signals
df_eqtl_wide_X <- readRDS(indep_x_rds)
df_eqtl_wide_PAR <- readRDS(indep_par_rds)
df_eqtl_wide_X <- df_eqtl_wide_X %>%
  select("#CHROM", POS, ID, rsid, REF, ALT,
         cell_type, gene, par_type,
         setting = min_pvalue_model)
df_eqtl_wide_PAR <- df_eqtl_wide_PAR %>%
  select("#CHROM", POS, ID, rsid, REF, ALT,
         cell_type, gene, par_type,
         setting = min_pvalue_model)
df_eqtl_wide_X_PAR <- rbind(df_eqtl_wide_X,df_eqtl_wide_PAR)
df_eqtl_wide_X_PAR$independent_signal <- T

df_X_PAR_pass_two_step_fdr <- left_join(df_X_PAR_pass_two_step_fdr,df_eqtl_wide_X_PAR)
df_X_PAR_pass_two_step_fdr$MARKER <- paste0(df_X_PAR_pass_two_step_fdr$gene)
df_X_PAR_pass_two_step_fdr_plot <- select(df_X_PAR_pass_two_step_fdr,
                                          POS,
                                          cell_type,
                                          LOG10_P,
                                          MARKER,
                                          setting,
                                          independent_signal)
# 3. Plot
df <- df_X_PAR_pass_two_step_fdr_plot

# Labels: independent lead variants only; one row per gene (its strongest signal);
# one gene per locus (drop genes within min_gap_kb of a stronger selected gene); top n_label
df_label <- df %>%
  filter(independent_signal == TRUE) %>%
  group_by(MARKER) %>%
  slice_max(LOG10_P, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(LOG10_P))

keep <- rep(FALSE, nrow(df_label))
for (k in seq_len(nrow(df_label))) {
  keep[k] <- !any(keep & abs(df_label$POS[keep] - df_label$POS[k]) < min_gap_kb * 1000)
}
df_label <- df_label[keep, ] %>% slice_head(n = n_label)

p <- ggplot(df, aes(x = POS, y = LOG10_P, color = cell_type, shape = setting)) +

  annotate(
    "rect",
    xmin = par1_start, xmax = par1_end,
    ymin = -Inf, ymax = Inf,
    fill = "brown", alpha = 0.15
  ) +

  annotate(
    "rect",
    xmin = par2_start, xmax = par2_end,
    ymin = -Inf, ymax = Inf,
    fill = "brown", alpha = 0.15
  ) +

  geom_point(size = 1.2, alpha = 0.8) +

  geom_text_repel(
    data = df_label,
    aes(label = MARKER),
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +

  scale_color_manual(values = cell_type_colors) +

  scale_shape_manual(
    values = c(
      "both_par" = 16,
      "both_model_2" = 16,
      "both_model_1" = 17,
      "female_par" = 18,
      "female_xchr_model_2" = 18,
      "male_par" = 15,
      "male_xchr_model_1" = 15
    )
  ) +

  scale_x_continuous(
    limits = c(0, chrX_len),
    breaks = seq(0, 1.5e8, by = 3e7),
    labels = seq(0, 150, by = 30),
    expand = c(0, 0)
  ) +

  labs(
    x = "Position on ChrX",
    y = "-log10(P)",
    color = "Cell Type",
    shape = "Setting",
    title = "Cell-type specific eQTL Manhattan plot (ChrX)"
  ) +

  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
p
ggsave("manhattanplot.pdf",p,width = 16,height = 5, dpi=300)
