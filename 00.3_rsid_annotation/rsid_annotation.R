# rsID (dbSNP) and nearest protein-coding gene (GENCODE) annotation for chrX variants; needs plink2 + bcftools
rm(list = ls())
library(dplyr)
library(data.table)
library(stringr)
# CONFIG
out_dir      <- getwd()
dbsnp_full   <- "path/to/dbsnp/GCF_000001405.40.gz"   # RefSeq-named contigs
sn_cov_file  <- "path/to/single_nucleus_covariates.tsv"   # used only for its sample columns
wgs_fam      <- "path/to/wgs_genotypes.fam"
wgs_bfile    <- "path/to/wgs_genotypes"
gencode_gtf  <- "path/to/gencode.basic.annotation.gtf.gz"
out_file     <- "chrX_QC"
dbsnp_chrX_acc <- "NC_000023.11"   # chrX RefSeq accession (GRCh38)
# 0. Extract dbSNP chrX and rename contig to X
run <- function(cmd, args) {
  message("\n>>> ", cmd, " ", paste(args, collapse = " "))
  status <- system2(cmd, args = args)
  if (status != 0L) {
    stop(sprintf("Command failed (exit %d): %s %s",
                 status, cmd, paste(args, collapse = " ")), call. = FALSE)
  }
}
chrX_vcf   <- file.path(out_dir, "dbsnp_chrX.vcf.gz")
rename_txt <- file.path(out_dir, "refseq_to_X.txt")
xname_vcf  <- file.path(out_dir, "dbsnp_chrX.Xname.vcf.gz")
dbsnp_chrX <- xname_vcf
run("bcftools", c("view",
                  "-r", dbsnp_chrX_acc,
                  shQuote(dbsnp_full),
                  "-Oz",
                  "-o", shQuote(chrX_vcf)))
run("bcftools", c("index", "-t", shQuote(chrX_vcf)))
writeLines(paste(dbsnp_chrX_acc, "X", sep = "\t"), rename_txt)
run("bcftools", c("annotate",
                  "--rename-chrs", shQuote(rename_txt),
                  shQuote(chrX_vcf),
                  "-Oz",
                  "-o", shQuote(xname_vcf)))
run("bcftools", c("index", "-t", shQuote(xname_vcf)))
message("\nDone. Output: ", xname_vcf)
# 1. Sample IDs shared by snRNA and WGS
sn_cov_df <- fread(sn_cov_file)
sn_id <- data.frame(V1 = colnames(sn_cov_df)[-1])
fam <- fread(wgs_fam)
overlapping_id <- semi_join(fam, sn_id)
overlapping_id <- overlapping_id[, c("V1", "V2")]
fwrite(overlapping_id, "extracted_overlapping_ids.txt", sep = "\t", quote = FALSE, col.names = FALSE)
# 2. Extract chrX; re-ID variants as chr:bp:REF:ALT; drop duplicates and flipped alleles
output_file <- file.path(out_dir, out_file)
cmd <- paste0(
  "plink2 ",
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
df_bim <- fread(paste0(output_file, "_0.0.1.bim"))
num_variant <- as.data.frame(table(df_bim$V2))
num_variant_duplicate <- filter(num_variant, Freq > 1)
write.table(as.character(num_variant_duplicate$Var1),
            paste0(output_file, "_0.0.1_duplicates_to_be_removed.txt"),
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
cmd <- paste0(
  "plink2 ",
  "--bfile ", paste0(output_file, "_0.0.1"), " ",
  "--exclude ", paste0(output_file, "_0.0.1_duplicates_to_be_removed.txt"), " ",
  "--make-bed ",
  "--out ", paste0(output_file, "_0.0.2")
)
cat(cmd, "\n")
system(cmd)
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
  "plink2 ",
  "--bfile ", paste0(output_file, "_0.0.2"), " ",
  "--exclude ", paste0(output_file, "_0.0.2_variants_with_flipped_allele_to_be_removed.txt"), " ",
  "--make-bed ",
  "--out ", paste0(output_file, "_0.0.3")
)
cat(cmd, "\n")
system(cmd)
# 3. bfile -> VCF -> rsID annotation -> bfile
bfile_003  <- paste0(output_file, "_0.0.3")
vcf_prefix <- paste0(output_file, "_0.0.3")
cmd <- paste0(
  "plink2 ",
  "--bfile ", bfile_003, " ",
  "--recode vcf bgz ",
  "--keep-allele-order ",
  "--out ", vcf_prefix
)
cat(cmd, "\n")
system(cmd)
cmd <- paste0("bcftools index -t ", vcf_prefix, ".vcf.gz")
cat(cmd, "\n")
system(cmd)
vcf_rsid <- paste0(vcf_prefix, ".rsid.vcf.gz")
cmd <- paste0(
  "bcftools annotate ",
  "-a ", dbsnp_chrX, " ",
  "-c ID ",
  vcf_prefix, ".vcf.gz ",
  "-Oz -o ", vcf_rsid
)
cat(cmd, "\n")
system(cmd)
cmd <- paste0("bcftools index -t ", vcf_rsid)
cat(cmd, "\n")
system(cmd)
bfile_rsid <- paste0(vcf_prefix, ".rsid")
cmd <- paste0(
  "plink2 ",
  "--vcf ", vcf_rsid, " ",
  "--make-bed ",
  "--keep-allele-order ",
  "--double-id ",
  "--out ", sub("\\.vcf\\.gz$", "", vcf_rsid)
)
cat(cmd, "\n")
system(cmd)
# 4. Nearest protein-coding gene (GENCODE)
df_bim_old <- fread(paste0(output_file, "_0.0.3.bim"))
df_bim_new <- fread(paste0(output_file, "_0.0.3.rsid.bim"))
df_bim_old$ID <- df_bim_old$V2
df_bim_new$ID <- paste0("chrX:", df_bim_new$V4, ":", df_bim_new$V6, ":", df_bim_new$V5)
d <- semi_join(df_bim_new, df_bim_old, by = 'ID')   # check nrow
d <- semi_join(df_bim_old, df_bim_new, by = 'ID')   # check nrow
ref_gencode <- fread(gencode_gtf)
setnames(
  ref_gencode,
  names(ref_gencode)[1:9],
  c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
)
ref_gencode <- filter(ref_gencode, feature == 'gene')
ref_gencode <- filter(ref_gencode, seqname == 'chrX')
ref_gencode <- filter(ref_gencode, grepl("protein_coding", attribute))
ref_gencode <- ref_gencode %>%
  mutate(
    gene_id = str_extract(attribute, 'gene_id "[^"]+"') %>%
      str_replace('gene_id "', '') %>%
      str_replace('"', ''),
    gene_name = str_extract(attribute, 'gene_name "[^"]+"') %>%
      str_replace('gene_name "', '') %>%
      str_replace('"', '')
  )
ref_gencode <- select(ref_gencode, start, end, gene_id, gene_name)
colnames(ref_gencode) <- c("start", "end", "gene_id", "gene_symbol")
# genic: variant inside a gene; intergenic: nearest gene assigned
df_bim_new$V4     <- as.integer(df_bim_new$V4)
ref_gencode$start <- as.integer(ref_gencode$start)
ref_gencode$end   <- as.integer(ref_gencode$end)
ord_gene <- order(ref_gencode$start)
ref_gene <- ref_gencode[ord_gene, ]
gene_start <- ref_gene$start
gene_end   <- ref_gene$end
gene_id    <- ref_gene$gene_id
gene_sym   <- ref_gene$gene_symbol
n_gene     <- length(gene_start)
pos <- df_bim_new$V4
n   <- length(pos)
df_bim_new$gene_id     <- NA_character_
df_bim_new$gene_symbol <- NA_character_
df_bim_new$region      <- NA_character_
idx <- findInterval(pos, gene_start)
valid_idx <- idx >= 1 & idx <= n_gene
genic <- rep(FALSE, n)
genic[valid_idx] <- pos[valid_idx] <= gene_end[idx[valid_idx]]
df_bim_new$gene_id[genic]     <- gene_id[idx[genic]]
df_bim_new$gene_symbol[genic] <- gene_sym[idx[genic]]
df_bim_new$region[genic]      <- "genic"
inter_idx <- which(!genic)
idx_i     <- idx[inter_idx]
pos_i     <- pos[inter_idx]
left  <- idx_i
right <- idx_i + 1
dist_left <- rep(Inf, length(pos_i))
ok_left   <- left >= 1 & left <= n_gene
dist_left[ok_left] <-
  pmax(0L,
       gene_start[left[ok_left]] - pos_i[ok_left],
       pos_i[ok_left] - gene_end[left[ok_left]])
dist_right <- rep(Inf, length(pos_i))
ok_right   <- right >= 1 & right <= n_gene
dist_right[ok_right] <-
  pmax(0L,
       gene_start[right[ok_right]] - pos_i[ok_right],
       pos_i[ok_right] - gene_end[right[ok_right]])
use_left <- dist_left <= dist_right
nearest  <- ifelse(use_left, left, right)
df_bim_new$gene_id[inter_idx]     <- gene_id[nearest]
df_bim_new$gene_symbol[inter_idx] <- gene_sym[nearest]
df_bim_new$region[inter_idx]      <- "intergenic"
df_rsid_annotation <- select(df_bim_new, ID, V2, gene_id, gene_symbol, region)
colnames(df_rsid_annotation) <- c("ID", "rsid", "gene_id", "gene_symbol", "region")
saveRDS(df_rsid_annotation, "rsid_gene_annotation_full.rds")
df_rsid_annotation <- dplyr::select(df_rsid_annotation, ID, rsid)
saveRDS(df_rsid_annotation, "rsid_gene_annotation.rds")