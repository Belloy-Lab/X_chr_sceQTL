# Apply remove-lists (f3-f5), sex-aware MAF >= 1% filter, final allele frequencies
rm(list = ls())
library(data.table)
# CONFIG
out_file <- "chrX_QC"
maf_cut  <- 0.01
# 1. Apply remove-lists
command <- paste("plink2",
                 "--bfile", paste0(out_file, "_EU.0.E.2"),
                 "--exclude", paste0(out_file, "_EU.0.E_cn_fish_test_sex__remove.txt"),
                 "--make-bed --out", paste0(out_file, "_EU.0.F.1"),
                 sep = " ")
system(command)
command <- paste("plink2",
                 "--bfile", paste0(out_file, "_EU.0.F.1"),
                 "--exclude", paste0(out_file, "_EU.0.E.Gnomad3_LCR__remove.txt"),
                 "--make-bed --out", paste0(out_file, "_EU.0.F.2"),
                 sep = " ")
system(command)
command <- paste("plink2",
                 "--bfile", paste0(out_file, "_EU.0.F.2"),
                 "--exclude", paste0(out_file, "_EU.0.E.Gnomad3_nonPASS__remove.txt"),
                 "--make-bed --out", paste0(out_file, "_EU.0.F.3"),
                 sep = " ")
system(command)
command <- paste("plink2",
                 "--bfile", paste0(out_file, "_EU.0.F.3"),
                 "--exclude", paste0(out_file, "_EU.0.E.Gnomad3_10%_MAF_diff__remove.txt"),
                 "--make-bed --out", paste0(out_file, "_EU.0.F.4"),
                 sep = " ")
system(command)
# 2. Sex-aware MAF filter: MAF >= cutoff in one stratum and > 0 in the others
cmd_both <- paste(
  "plink2",
  "--bfile", paste0(out_file, "_EU.0.F.4"),
  "--freq",
  "--out", paste0(out_file, "_both")
)
system(cmd_both)
cmd_male <- paste(
  "plink2",
  "--bfile", paste0(out_file, "_EU.0.F.4"),
  "--keep-males",
  "--freq",
  "--out", paste0(out_file, "_male")
)
system(cmd_male)
cmd_female <- paste(
  "plink2",
  "--bfile", paste0(out_file, "_EU.0.F.4"),
  "--keep-females",
  "--freq",
  "--out", paste0(out_file, "_female")
)
system(cmd_female)
both_freq   <- fread(paste0(out_file, "_both.afreq"))
male_freq   <- fread(paste0(out_file, "_male.afreq"))
female_freq <- fread(paste0(out_file, "_female.afreq"))
df_keep <- merge(
  merge(
    both_freq[, .(ID, MAF_both = pmin(ALT_FREQS, 1 - ALT_FREQS))],
    male_freq[, .(ID, MAF_male = pmin(ALT_FREQS, 1 - ALT_FREQS))],
    by = "ID",
    all = TRUE
  ),
  female_freq[, .(ID, MAF_female = pmin(ALT_FREQS, 1 - ALT_FREQS))],
  by = "ID",
  all = TRUE
)
df_keep <- df_keep[
  MAF_both   >= maf_cut &
  MAF_male   >= maf_cut &
  MAF_female >= maf_cut
]
fwrite(
  df_keep[, .(ID)],
  paste0(out_file, "_maf_both_sex.txt"),
  sep = "\t",
  col.names = FALSE
)
cmd_extract <- paste(
  "plink2",
  "--bfile", paste0(out_file, "_EU.0.F.4"),
  "--extract", paste0(out_file, "_maf_both_sex.txt"),
  "--make-bed",
  "--out", paste0(out_file, "_EU.0.G")
)
system(cmd_extract)
fam_file <- paste0(out_file, "_EU.0.G.fam")
fam <- fread(fam_file, header = FALSE)
fam$V6 <- -9
fwrite(fam, fam_file, sep = " ", col.names = FALSE)
# 3. Final allele frequencies
command <- paste(
  "plink2",
  "--bfile", paste0(out_file, "_EU.0.G"),
  "--freq",
  "--out", paste0(out_file, "_EU.0.G")
)
system(command)
command <- paste(
  "plink2",
  "--bfile", paste0(out_file, "_EU.0.G"),
  "--filter-females",
  "--freq",
  "--out", paste0(out_file, "_EU.0.G.female")
)
system(command)
command <- paste(
  "plink2",
  "--bfile", paste0(out_file, "_EU.0.G"),
  "--filter-males",
  "--freq",
  "--out", paste0(out_file, "_EU.0.G.male")
)
system(command)
