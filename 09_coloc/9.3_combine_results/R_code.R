# Combine coloc results; keep eQTL settings valid for each XWAS model; highest PP4 per trait x cell type x model x gene
rm(list = ls())
library(data.table)
library(tidyverse)

coloc_res_path <- dir("../9.2_coloc_analysis/",pattern = "coloc_summary_job")
coloc_res_path <- paste0("../9.2_coloc_analysis/",coloc_res_path)

coloc_res <- rbindlist(lapply(coloc_res_path, fread), fill = TRUE)
setnames(coloc_res, "gene", "Gene")
setnames(coloc_res, "eQTL_type", "eqtl_setting_coloc")
setnames(coloc_res, "GWAS_type", "gwas_model")


coloc_res <- coloc_res %>%
  filter(
    (gwas_model == "eXCI"   & eqtl_setting_coloc %in% c("female_xchr_model_2", "both_model_1")) |
      (gwas_model == "female" & eqtl_setting_coloc %in% c("female_xchr_model_2", "both_model_2", "both_model_1")) |
      (gwas_model == "male"   & eqtl_setting_coloc %in% c("male_xchr_model_1", "both_model_2")) |
      (gwas_model == "rXCI"   & eqtl_setting_coloc %in% c("both_model_2"))
  )

coloc_res$trait <- sub("_(?=[^_]+$).*", "", coloc_res$trait, perl = TRUE)


coloc_res_distinct <- coloc_res %>%
  arrange(trait, cell_type, desc(PP4)) %>%
  group_by(trait, cell_type) %>%
  distinct(gwas_model, Gene, .keep_all = TRUE) %>%
  ungroup()

saveRDS(coloc_res_distinct,"coloc_res_distinct.rds")
writexl::write_xlsx(coloc_res_distinct,"coloc_res_distinct.xlsx")