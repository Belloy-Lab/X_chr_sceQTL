# Combine SMR results; keep eQTL settings valid for each XWAS model; most significant hit per trait x cell type x model x gene
rm(list = ls())
library(tidyverse)
# CONFIG
smr_res_dir <- "../8.3_smr_analysis/"
out_rds     <- "smr_res_distinct.rds"
# 1. Result files: smr_res_<trait>_<gwas_model>.smr
files <- list.files(smr_res_dir, pattern = "\\.smr$", full.names = TRUE)
name <- sub("\\.smr$", "", sub("^smr_res_", "", basename(files)))
smr_res_path <- data.frame(
  path       = files,
  trait      = sub("_[^_]+$", "", name),
  gwas_model = sub(".*_", "", name),
  stringsAsFactors = FALSE
)
# 2. Combine
smr_res <- purrr::map_dfr(seq_len(nrow(smr_res_path)), function(i) {
  df <- read.table(smr_res_path$path[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df$trait      <- smr_res_path$trait[i]
  df$gwas_model <- smr_res_path$gwas_model[i]
  df
})
# 3. probeID = <gene>_<cell_type>_<setting>
smr_res$cell_type        <- sub("^[^_]+_([^_]+)_.*", "\\1", smr_res$probeID)
smr_res$eqtl_setting_smr <- sub("^[^_]+_[^_]+_", "", smr_res$probeID)
# 4. Valid eQTL settings per XWAS model
smr_res <- smr_res %>%
  filter(
    (gwas_model == "eXCI"   & eqtl_setting_smr %in% c("female_xchr_model_2", "both_model_1")) |
      (gwas_model == "female" & eqtl_setting_smr %in% c("female_xchr_model_2", "both_model_2", "both_model_1")) |
      (gwas_model == "male"   & eqtl_setting_smr %in% c("male_xchr_model_1", "both_model_2")) |
      (gwas_model == "rXCI"   & eqtl_setting_smr %in% c("both_model_2"))
  )
# 5. Smallest p_SMR per trait x cell type x model x gene
smr_res_distinct <- smr_res %>%
  arrange(trait, cell_type, p_SMR) %>%
  group_by(trait, cell_type) %>%
  distinct(gwas_model, Gene, .keep_all = TRUE) %>%
  ungroup()
saveRDS(smr_res_distinct, out_rds)