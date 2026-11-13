# Independent non-PAR X signals: combine clumped significant eQTLs across 4 settings (wide),
# second LD clumping, sex-heterogeneity test (BH-FDR on independent lead variants), Fig 2A-C
rm(list = ls())
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpattern)
library(ieugwasr)
library(purrr)
# CONFIG
cis_all_rds <- "../03_merge_eQTL_res_X/df_eqtl_cis_all.rds"
rsid_annotation_rds <- "../00.3_rsid_annotation/rsid_gene_annotation.rds"
bfile_prefix      <- "../00.1_X_chr_bfile_QC/chrX_QC_EU.0.G"
afreq_both_path   <- paste0(bfile_prefix, ".afreq")
afreq_female_path <- paste0(bfile_prefix, ".female.afreq")
afreq_male_path   <- paste0(bfile_prefix, ".male.afreq")
clump_dir              <- "../04_significant_eQTL_X/"
clump_both_model_1_rds <- paste0(clump_dir, "eqtl_sig_clumped_all_X_chr_2stepFDR_after_prune_both_model_1.rds")
clump_both_model_2_rds <- paste0(clump_dir, "eqtl_sig_clumped_all_X_chr_2stepFDR_after_prune_both_model_2.rds")
clump_female_rds       <- paste0(clump_dir, "eqtl_sig_clumped_all_X_chr_2stepFDR_after_prune_female_xchr_model_2.rds")
clump_male_rds         <- paste0(clump_dir, "eqtl_sig_clumped_all_X_chr_2stepFDR_after_prune_male_xchr_model_1.rds")
# Fig 2B: eQTLs marked with * (and labeled) at this sex-het FDR threshold
fig2b_label_fdr <- 0.1
# 1. Annotate cis-eQTL results with rsID, gene, allele frequency
df_eqtl_cis_all <- readRDS(cis_all_rds)
df_eqtl_cis_all$fdr <- NULL
df_rsid_gene_annotation <- readRDS(rsid_annotation_rds)
df_eqtl_cis_all <- left_join(df_eqtl_cis_all,df_rsid_gene_annotation,by='ID')
afreq_both <- fread(afreq_both_path)
afreq_both$OBS_CT <- NULL
colnames(afreq_both)[5] <- "EAF_both"
afreq_female <- fread(afreq_female_path)
afreq_female$OBS_CT <- NULL
colnames(afreq_female)[5] <- "EAF_female"
afreq_male <- fread(afreq_male_path)
afreq_male$OBS_CT <- NULL
colnames(afreq_male)[5] <- "EAF_male"
df_eqtl_cis_all <- left_join(df_eqtl_cis_all,afreq_both,by=c("#CHROM","ID", "REF", "ALT"))
df_eqtl_cis_all <- left_join(df_eqtl_cis_all,afreq_female,by=c("#CHROM","ID", "REF", "ALT"))
df_eqtl_cis_all <- left_join(df_eqtl_cis_all,afreq_male,by=c("#CHROM","ID", "REF", "ALT"))
# 2. Keep cell-type-gene-variant pairs significant in >= 1 setting; pivot wide across settings
df_clump_sig_both_model_1 <- readRDS(clump_both_model_1_rds)
df_clump_sig_both_model_2 <- readRDS(clump_both_model_2_rds)
df_clump_sig_female_model_2 <- readRDS(clump_female_rds)
df_clump_sig_male_model_1 <- readRDS(clump_male_rds)
df_clump_sig_all_model <- rbind(df_clump_sig_both_model_1,
                                df_clump_sig_both_model_2,
                                df_clump_sig_female_model_2,
                                df_clump_sig_male_model_1)

df_eqtl_cis_clump_sig_all_model <- semi_join(df_eqtl_cis_all,df_clump_sig_all_model,
                                             by=c("ID","cell_type","gene"))
df_eqtl_cis_clump_sig_all_model <- df_eqtl_cis_clump_sig_all_model %>%
  left_join(
    df_clump_sig_all_model %>%
      select(ID, cell_type, gene, setting, significant_by_2step_FDR),
    by = c("ID", "cell_type", "gene", "setting")
  )
df_eqtl_cis_clump_sig_all_model$TEST <- NULL
df_eqtl_cis_clump_sig_all_model$ERRCODE <- NULL
# A1 should equal ALT
identical(df_eqtl_cis_clump_sig_all_model$A1,df_eqtl_cis_clump_sig_all_model$ALT)
df_eqtl_cis_clump_sig_all_model$A1 <- NULL

id_cols <- c(
  "#CHROM", "POS", "ID", "rsid","REF", "ALT",
  "cell_type", "gene", "par_type",
  "EAF_both", "EAF_female", "EAF_male"
)
value_cols <- c(
  "OBS_CT", "BETA", "SE", "T_STAT", "LOG10_P", "pvalue",
  "significant_by_2step_FDR"
)
df_eqtl_wide <- df_eqtl_cis_clump_sig_all_model %>%
  pivot_wider(
    id_cols   = all_of(id_cols),
    names_from  = setting,
    values_from = all_of(value_cols),
    names_glue  = "{.value}_{setting}"
  )
# 3. Min p-value across settings, male/female beta ratio, sex-heterogeneity test
p_cols <- grep("^pvalue_", colnames(df_eqtl_wide), value = TRUE)
df_eqtl_wide$min_pvalue <- apply(
  df_eqtl_wide[, p_cols],
  1,
  function(x) min(x, na.rm = TRUE)
)
df_eqtl_wide$min_pvalue_model <- apply(
  df_eqtl_wide[, p_cols],
  1,
  function(x) {
    p_cols[which.min(x)]
  }
)
df_eqtl_wide$min_pvalue_model <- sub("pvalue_", "", df_eqtl_wide$min_pvalue_model)

df_eqtl_wide$male_female_ratio <- (df_eqtl_wide$BETA_male_xchr_model_1)/(df_eqtl_wide$BETA_female_xchr_model_2)
df_eqtl_wide$abs_male_female_ratio <- abs(df_eqtl_wide$male_female_ratio)

# classify sex-effect pattern given a (nominal or FDR) p-value
label_sex_het <- function(p, opposite_dir, delta_abs, alpha = 0.05) {
  case_when(
    is.na(p) ~ NA_character_,
    p >= alpha ~ "No sex heterogeneity",
    !opposite_dir & delta_abs > 0 ~ "Male-biased (concordant)",
    !opposite_dir & delta_abs < 0 ~ "Female-biased (concordant)",
    opposite_dir & delta_abs > 0 ~ "Male-dominant (opposite)",
    opposite_dir & delta_abs < 0 ~ "Female-dominant (opposite)",
    TRUE ~ "Unclassified"
  )
}

eps <- 1e-6
df_eqtl_wide <- df_eqtl_wide %>%
  mutate(
    # male beta/SE (0/2 coding) rescaled to the female 0/1/2 scale
    male_beta_for_het = BETA_male_xchr_model_1 / 2,
    male_se_for_het   = SE_male_xchr_model_1 / 2,

    Z_het = (male_beta_for_het - BETA_female_xchr_model_2) /
      sqrt(male_se_for_het^2 + SE_female_xchr_model_2^2),
    P_het = 2 * pnorm(-abs(Z_het)),
    # BH-FDR is computed in step 4b on the final independent lead variants
    P_het_FDR = NA_real_,

    male_sign   = sign(ifelse(abs(male_beta_for_het) < eps, 0, male_beta_for_het)),
    female_sign = sign(ifelse(abs(BETA_female_xchr_model_2) < eps, 0, BETA_female_xchr_model_2)),
    opposite_dir = (male_sign * female_sign) == -1,

    delta_abs = abs(male_beta_for_het) - abs(BETA_female_xchr_model_2)
  ) %>%
  mutate(
    sex_het_label_nominal = label_sex_het(P_het, opposite_dir, delta_abs),
    sex_het_label_fdr     = NA_character_
  )

saveRDS(df_eqtl_wide,"df_eqtl_wide_before_second_clump.rds")
# 4. Second LD clumping across settings
column_order <- colnames(df_eqtl_wide)

df_eqtl_wide$snp <- df_eqtl_wide$rsid
df_eqtl_wide$rsid <- df_eqtl_wide$ID
df_eqtl_wide$pval <- df_eqtl_wide$min_pvalue
df_eqtl_wide$id <- paste0(df_eqtl_wide$gene,"_",df_eqtl_wide$cell_type)

df_eqtl_wide <- ld_clump(
  df_eqtl_wide,
  clump_kb  = 1000,
  clump_r2  = 0.1,
  clump_p   = 1,
  bfile     = bfile_prefix,
  plink_bin = plinkbinr::get_plink_exe()
)

df_eqtl_wide$rsid <- df_eqtl_wide$snp
df_eqtl_wide$snp <- NULL
df_eqtl_wide$pval <- NULL
df_eqtl_wide$id <- NULL
df_eqtl_wide <- select(df_eqtl_wide,column_order)

# 4b. Sex-heterogeneity BH-FDR on the final independent lead variants
#     (all variants with both male and female results, pooled across cell types)
het_idx <- !is.na(df_eqtl_wide$P_het)
df_eqtl_wide$P_het_FDR[het_idx] <- p.adjust(df_eqtl_wide$P_het[het_idx], method = "BH")
df_eqtl_wide$sex_het_label_fdr <- label_sex_het(df_eqtl_wide$P_het_FDR,
                                                df_eqtl_wide$opposite_dir,
                                                df_eqtl_wide$delta_abs)

saveRDS(df_eqtl_wide,"df_eqtl_wide_after_second_clump.rds")
# Shared layers for Fig 2A/2B: reference lines for the expected log2(|beta_male| / |beta_female|)
# (male coded 0/1, female 0/1/2): 1 = rXCI expectation, 0 = eXCI expectation
xci_ref <- data.frame(y = c(1, 0), model = c("rXCI", "eXCI"))
xci_ref_layers <- list(
  geom_hline(data = xci_ref, aes(yintercept = y),
             linetype = "dashed", linewidth = 0.3, color = "grey20"),
  scale_y_continuous(
    breaks   = seq(-5, 10, by = 1),
    sec.axis = sec_axis(~ ., breaks = xci_ref$y, labels = xci_ref$model)
  ),
  theme(axis.text.y.right  = element_text(size = 11, face = "bold", color = "grey20"),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right  = element_blank())
)
cell_levels <- c("Exc","Inh","Ast","Oli","OPC","Mic","End")
# 5. Fig 2A: per-cell-type bar plot of log2(|male/female beta ratio|), colored by best model
df_eqtl_wide_plot <- readRDS("df_eqtl_wide_after_second_clump.rds") %>%
  filter(!is.na(sex_het_label_nominal)) %>%
  arrange(desc(abs_male_female_ratio)) %>%
  mutate(log_ratio = log2(abs_male_female_ratio),
         cell_type = factor(cell_type, levels = cell_levels)) %>%
  arrange(cell_type) %>%
  mutate(idx = row_number())

cell_breaks <- df_eqtl_wide_plot %>%
  group_by(cell_type) %>%
  summarise(max_idx = max(idx)) %>%
  pull(max_idx)
cell_breaks <- c(0.5, cell_breaks + 0.5)

ggplot(df_eqtl_wide_plot, aes(x = idx, y = log_ratio, fill = min_pvalue_model)) +
  geom_col(width = 1, color = "black", linewidth = 0.2) +
  xci_ref_layers +
  scale_fill_manual(
    values = c(
      "male_xchr_model_1"   = "#1f78b4",
      "female_xchr_model_2" = "#e31a1c",
      "both_model_2"        = "#33a02c",
      "both_model_1"        = "#ff7f00"
    ),
    breaks = c("male_xchr_model_1", "female_xchr_model_2", "both_model_1", "both_model_2"),
    labels = c("Male", "Female", "eXCI", "rXCI"),
    name = "Most significant model"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    legend.position = c(0.02, 0),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2)
  ) +
  labs(title = " ", x = NULL, y = "log2 absolute value of\nMale / Female beta ratio") +
  geom_vline(xintercept = cell_breaks, linetype = "dashed", linewidth = 0.15, color = "grey50")

# variants without a male result are not shown
ggsave("Barplot_best_model_7_cell_type_panel_WIDE.pdf",width = 11,height = 3.9,dpi = 300)
# 6. Fig 2B inset: variants per sex-het category, striped overlay = MAF >= 0.05 in all strata
# NA EAF_male counts as not common
df_eqtl_wide_plot <- df_eqtl_wide_plot %>%
  mutate(
    maf_both    = pmin(EAF_both,   1 - EAF_both),
    maf_female  = pmin(EAF_female, 1 - EAF_female),
    maf_male    = pmin(EAF_male,   1 - EAF_male),
    common_all3 = (maf_both >= 0.05 & maf_female >= 0.05 & maf_male >= 0.05) %in% TRUE
  )

sex_het_colors_5 <- c(
  "Male-biased (concordant)"   = "#1f78b4",
  "Female-biased (concordant)" = "#e31a1c",
  "Male-dominant (opposite)"   = "#08306b",
  "Female-dominant (opposite)" = "#67000d",
  "No sex heterogeneity"       = "grey70"
)

# levels reversed: y axis runs bottom-up
sexhet_levels <- rev(c(
  "Male-biased (concordant)",
  "Male-dominant (opposite)",
  "Female-biased (concordant)",
  "Female-dominant (opposite)",
  "No sex heterogeneity"
))

striped_hbar <- function(bar_df, y_var, fill_values) {
  ggplot(bar_df, aes(y = .data[[y_var]])) +
    geom_col(aes(x = n_total, fill = .data[[y_var]]),
             width = 0.7, color = "black", linewidth = 0.2) +
    geom_col_pattern(
      aes(x = n_common),
      fill            = NA,
      color           = "black",
      linewidth       = 0.2,
      width           = 0.7,
      pattern         = "stripe",
      pattern_fill    = "black",
      pattern_color   = "black",
      pattern_angle   = 45,
      pattern_density = 0.1,
      pattern_spacing = 0.02
    ) +
    geom_text(aes(x = n_total, label = n_total), hjust = -0.3, size = 4, fontface = "bold") +
    scale_fill_manual(values = fill_values) +
    scale_y_discrete(drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    theme_classic(base_size = 14) +
    labs(x = "Number of independent lead variants", y = NULL) +
    theme(legend.position = "none",
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank())
}

bar_df2 <- df_eqtl_wide_plot %>%
  group_by(sex_het_label_nominal) %>%
  summarise(n_total = n(), n_common = sum(common_all3), .groups = "drop") %>%
  mutate(sex_het_label_nominal = factor(sex_het_label_nominal, levels = sexhet_levels))

p_bar2 <- striped_hbar(bar_df2, "sex_het_label_nominal", sex_het_colors_5)
ggsave("barplot_sex_het_horizontal_MAF.pdf", p_bar2, width = 5, height = 2.7, dpi = 300)
# 7. Fig 2A inset: variants per best model, striped overlay = MAF >= 0.05 in all strata
bar_df <- df_eqtl_wide_plot %>%
  group_by(min_pvalue_model) %>%
  summarise(n_total = n(), n_common = sum(common_all3), .groups = "drop") %>%
  mutate(
    model_label = recode(min_pvalue_model,
                         "both_model_1"        = "eXCI",
                         "both_model_2"        = "rXCI",
                         "female_xchr_model_2" = "Female",
                         "male_xchr_model_1"   = "Male"),
    model_label = factor(model_label, levels = c("rXCI", "eXCI", "Female", "Male"))
  )

p_bar <- striped_hbar(bar_df, "model_label",
                      c("Male" = "#1f78b4", "Female" = "#e31a1c",
                        "rXCI" = "#33a02c", "eXCI" = "#ff7f00"))
ggsave("barplot_best_model_variant_horizontal_MAF.pdf", p_bar, width = 5, height = 2.7, dpi = 300)
# 8. Fig 2B: same bars as Fig 2A colored by sex-het category (P_het < 0.05);
#    * = P_het FDR < fig2b_label_fdr (these eQTLs are labeled with their eGene)
ggplot(df_eqtl_wide_plot, aes(x = idx, y = log_ratio, fill = sex_het_label_nominal)) +
  geom_col(width = 1, color = "black", linewidth = 0.2) +
  xci_ref_layers +
  scale_fill_manual(values = sex_het_colors_5, name = "Sex-effect pattern (P-HET<0.05)") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    legend.position = c(0.02, 0),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2)
  ) +
  labs(title = " ", x = NULL, y = "log2 absolute value of\nMale / Female beta ratio") +
  geom_vline(xintercept = cell_breaks, linetype = "dashed", linewidth = 0.15, color = "grey50") +
  geom_text(
    data = df_eqtl_wide_plot %>%
      filter(P_het_FDR < fig2b_label_fdr) %>%
      mutate(y_star = ifelse(log_ratio >= 0, log_ratio + 0.05, log_ratio - 0.15)),
    aes(x = idx, y = y_star, label = "*"),
    size = 2
  )

ggsave("Barplot_sex_het_7_cell_type_panel_WIDE.pdf",width = 11,height = 3.9,dpi = 300)

# eGenes labeled in Fig 2B
fig2b_label_genes <- df_eqtl_wide_plot %>%
  filter(P_het_FDR < fig2b_label_fdr) %>%
  group_by(cell_type, gene) %>%
  summarise(min_P_het_FDR = min(P_het_FDR), .groups = "drop") %>%
  arrange(cell_type, min_P_het_FDR)
print(fig2b_label_genes, n = Inf)
# 9. Fig 2C: unique sex-het eGenes per cell type (male-/female-biased) at 3 thresholds
df_sexhet <- df_eqtl_wide_plot %>%
  mutate(direction = ifelse(delta_abs > 0, "Male", "Female"))

thresholds <- list(
  "P_het<0.05"   = df_sexhet$P_het     < 0.05,
  "FDR_het<0.1"  = df_sexhet$P_het_FDR < 0.1,
  "FDR_het<0.05" = df_sexhet$P_het_FDR < 0.05
)

sexhet_rows <- imap_dfr(thresholds, function(sel, nm) {
  df_sexhet[sel %in% TRUE, ] %>%
    transmute(threshold = nm, cell_type, gene, direction, rsid, P_het, P_het_FDR)
}) %>%
  mutate(threshold = factor(threshold, levels = names(thresholds)),
         direction = factor(direction, levels = c("Male", "Female")))

fig2c_counts <- sexhet_rows %>%
  distinct(threshold, cell_type, direction, gene) %>%
  count(threshold, cell_type, direction, .drop = FALSE) %>%
  pivot_wider(names_from = c(cell_type, direction), values_from = n, names_sep = "_")
print(fig2c_counts)
writexl::write_xlsx(list(Fig2B_labels = fig2b_label_genes, Fig2C_counts = fig2c_counts),
                    "Fig2_sex_het_eGenes.xlsx")
