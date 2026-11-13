# MAF distribution of all independent lead (index) eQTL variants (non-PAR + PAR, all settings)
rm(list = ls())
library(dplyr)
library(ggplot2)
library(ggpubr)
# CONFIG
indep_x_rds   <- "../05_independent_signal_X/df_eqtl_wide_after_second_clump.rds"
indep_par_rds <- "../05_independent_signal_PAR/df_eqtl_wide_after_second_clump.rds"
# MAF reference lines
maf_cut <- data.frame(
  maf   = c(0.01, 0.05, 0.10),
  label = c("MAF = 1%", "MAF = 5%", "MAF = 10%"),
  vjust = c(1.5, 3.5, 5.5)
)
# bin width on the log10 scale; density * bin width = proportion of variants per bin
bw_log <- 0.1
fill_col <- "grey50"
# 1. Independent lead variants after the second clumping (unique IDs; MAF from EAF_both)
df_index <- bind_rows(
  readRDS(indep_x_rds)   %>% select(ID, EAF_both),
  readRDS(indep_par_rds) %>% select(ID, EAF_both)
) %>%
  distinct(ID, .keep_all = TRUE) %>%
  mutate(MAF = pmin(EAF_both, 1 - EAF_both)) %>%
  filter(!is.na(MAF), MAF > 0)
cat(nrow(df_index), "index variants;",
    sprintf("MAF >= 5%%: %.1f%%", 100 * mean(df_index$MAF >= 0.05)), "\n")
# 2. Density plot
p <- ggplot(df_index, aes(x = MAF)) +
  geom_density(aes(y = after_stat(density * bw_log)),
               fill = fill_col, color = fill_col, alpha = 0.4) +
  scale_x_log10(
    breaks = c(0.005, 0.01, 0.05, 0.1, 0.5),
    labels = c("0.5", "1", "5", "10", "50")
  ) +
  geom_vline(data = maf_cut, aes(xintercept = maf),
             linetype = "dashed", color = "grey30") +
  annotate("text", x = maf_cut$maf, y = Inf, label = maf_cut$label,
           vjust = maf_cut$vjust, hjust = -0.05, size = 4.5, color = "black") +
  labs(x = "MAF (%, log10 scale)",
       y = "Proportion of variants") +
  theme_pubr()
ggsave("ggdensity_index.pdf", p, width = 3.675, height = 3.55, dpi = 300)
