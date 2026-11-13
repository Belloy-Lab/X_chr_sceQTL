# TWAS job table: one row per weights panel x XWAS; combined plus per-trait files
library(data.table)
library(dplyr)
# CONFIG
traits <- c("AD", "PD", "LBD", "MSA")
models <- c("rXCI", "eXCI", "female", "male")
sumstats_dir <- "../7.6_xwas_clean"
weights_dir  <- "../7.5_generate_weights/weights"
out_root     <- "."
ref_ld_chr   <- "../7.1_input_bfile/X_bfile"
log_root     <- file.path(out_root, "log")
# 1. XWAS map and weight-sex to XWAS-model pairing
df_gwas_map <- expand.grid(
  trait      = traits,
  gwas_model = models,
  stringsAsFactors = FALSE
)
df_gwas_map$sumstats <- file.path(
  sumstats_dir,
  paste0(df_gwas_map$trait, "_", df_gwas_map$gwas_model, ".txt")
)
df_rule <- data.frame(
  sex        = c("both", "both", "female", "both", "female", "both", "male"),
  gwas_model = c("rXCI", "eXCI", "eXCI",   "female", "female", "male", "male"),
  stringsAsFactors = FALSE
)
# 2. Panel table from .pos files
base_dirs <- data.frame(
  PAR_type   = "non_PAR",
  base_dir   = weights_dir,
  out_root   = out_root,
  REF_LD_CHR = ref_ld_chr,
  stringsAsFactors = FALSE
)
df_pos_all <- rbindlist(lapply(seq_len(nrow(base_dirs)), function(i) {
  base_dir   <- base_dirs$base_dir[i]
  out_root   <- base_dirs$out_root[i]
  par_type   <- base_dirs$PAR_type[i]
  ref_ld_chr <- base_dirs$REF_LD_CHR[i]
  
  pos_files <- list.files(
    path       = base_dir,
    pattern    = "\\.pos$",
    recursive  = TRUE,
    full.names = TRUE
  )
  
  if (length(pos_files) == 0) return(NULL)
  
  rel_path   <- sub(paste0("^", base_dir, "/"), "", pos_files)
  path_parts <- strsplit(rel_path, "/")
  
  data.frame(
    PAR_type    = par_type,
    pos_path    = pos_files,
    sex         = sapply(path_parts, `[`, 1),
    cell_type   = sapply(path_parts, `[`, 2),
    file_name   = sapply(path_parts, function(x) tail(x, 1)),
    weights_dir = dirname(pos_files),
    out_root    = out_root,
    REF_LD_CHR  = ref_ld_chr,
    stringsAsFactors = FALSE
  )
}))
# 3. Job table
df_job_all <- df_pos_all %>%
  inner_join(df_rule,     by = "sex",        relationship = "many-to-many") %>%
  inner_join(df_gwas_map, by = "gwas_model", relationship = "many-to-many") %>%
  mutate(
    chr        = 23,
    out_dir    = file.path(out_root, trait, gwas_model, sex, cell_type),
    out_prefix = paste0(trait, "-", cell_type, "-", sex, "-", gwas_model, "-WEIGHTS"),
    out_file   = file.path(out_dir, paste0(out_prefix, ".", chr, ".dat")),
    log_file   = file.path(
      log_root,
      paste0(PAR_type, "_", trait, "_", cell_type, "_", sex, "_", gwas_model, "_chr", chr, ".log")
    )
  ) %>%
  select(
    PAR_type, trait, sex, cell_type, gwas_model,
    sumstats, pos_path, weights_dir, REF_LD_CHR,
    chr, out_dir, out_prefix, out_file, log_file
  )
# no XWAS has PAR variants
df_job_all <- df_job_all %>% filter(PAR_type == "non_PAR")
write.table(df_job_all,
            "TWAS_job_reference_all.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
# 4. Split by trait
for (tr in unique(df_job_all$trait)) {
  df_sub   <- df_job_all %>% filter(trait == tr)
  out_file <- paste0("TWAS_job_reference_", tr, ".txt")
  write.table(df_sub, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Written:", out_file, " (n =", nrow(df_sub), ")\n")
}