# Build chrX bfile: intersect snRNA donors with WGS, extract chrX, re-ID variants, drop duplicates/flips, write sample lists
rm(list = ls())
library(dplyr)
library(data.table)
library(readxl)
# CONFIG
sn_cov_file   <- "path/to/single_nucleus_covariates.tsv"   # used only for its sample columns
wgs_bfile     <- "path/to/wgs_genotypes"
metadata_xlsx <- "path/to/metadata.xlsx"
plink2        <- "plink2"
out_prefix    <- "chrX_QC"
# 1. Sample IDs shared by snRNA and WGS
sn_cov_df <- fread(sn_cov_file)
sn_id <- data.frame(V1 = colnames(sn_cov_df)[-1])
fam <- fread(paste0(wgs_bfile, ".fam"))
overlapping_id <- semi_join(fam, sn_id)
overlapping_id <- overlapping_id[, c("V1", "V2")]
fwrite(overlapping_id, "extracted_overlapping_ids.txt", sep = "\t", quote = FALSE, col.names = FALSE)
# 2. Extract chrX; re-ID variants as chr:bp:REF:ALT (source IDs not unique for multi-allelics)
output_file <- file.path(getwd(), out_prefix)
cmd <- paste0(
  plink2, " ",
  "--bfile ", wgs_bfile, " ",
  "--chr X ",
  "--keep extracted_overlapping_ids.txt ",
  "--set-all-var-ids 'chr@:#:$r:$a' ",
  "--new-id-max-allele-len 10000 ",
  "--make-bed ",
  "--out ", paste0(output_file, "_0.0.1")
)
cat(cmd, "\n")
system(cmd)
# 3a. Remove duplicate variant IDs
df_bim <- fread(paste0(output_file, "_0.0.1.bim"))
num_variant <- as.data.frame(table(df_bim$V2))
num_variant_duplicate <- filter(num_variant, Freq > 1)
write.table(as.character(num_variant_duplicate$Var1),
            paste0(output_file, "_0.0.1_duplicates_to_be_removed.txt"),
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
cmd <- paste0(
  plink2, " ",
  "--bfile ", paste0(output_file, "_0.0.1"), " ",
  "--exclude ", paste0(output_file, "_0.0.1_duplicates_to_be_removed.txt"), " ",
  "--make-bed ",
  "--out ", paste0(output_file, "_0.0.2")
)
cat(cmd, "\n")
system(cmd)
# 3b. Remove same-position REF/ALT-flipped variants
df_bim <- fread(paste0(output_file, "_0.0.2.bim"))
bad_rows <- df_bim$V5 == df_bim$V6
if (any(bad_rows)) {
  stop(
    paste0(
      "ERROR: Found ", sum(bad_rows),
      " variants with ALT == REF. ",
      "Please clean the BIM file before proceeding."
    )
  )
}
key_forward <- with(df_bim, paste(V4, V5, V6, sep = ":"))
key_reverse <- with(df_bim, paste(V4, V6, V5, sep = ":"))
flipped_variants <- df_bim[key_reverse %in% key_forward, ]
write.table(flipped_variants$V2,
            paste0(output_file, "_0.0.2_variants_with_flipped_allele_to_be_removed.txt"),
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
cmd <- paste0(
  plink2, " ",
  "--bfile ", paste0(output_file, "_0.0.2"), " ",
  "--exclude ", paste0(output_file, "_0.0.2_variants_with_flipped_allele_to_be_removed.txt"), " ",
  "--make-bed ",
  "--out ", paste0(output_file, "_0.0.3")
)
cat(cmd, "\n")
system(cmd)
# 4. Sample lists: cogdx 1 = NCI, 2-6 = any CI; msex 1 = male, 0 = female
meta_data <- readxl::read_excel(metadata_xlsx)
id_pair <- function(specimen) data.frame(V1 = specimen, V2 = specimen)
fwrite(id_pair(filter(meta_data, cogdx == 1)$specimenID),
       "NCI_ids.txt",    sep = "\t", quote = FALSE, col.names = FALSE)
fwrite(id_pair(filter(meta_data, cogdx %in% 2:6)$specimenID),
       "any_CI_ids.txt", sep = "\t", quote = FALSE, col.names = FALSE)
fwrite(id_pair(filter(meta_data, msex == 1)$specimenID),
       "male_ids.txt",   sep = "\t", quote = FALSE, col.names = FALSE)
fwrite(id_pair(filter(meta_data, msex == 0)$specimenID),
       "female_ids.txt", sep = "\t", quote = FALSE, col.names = FALSE)
