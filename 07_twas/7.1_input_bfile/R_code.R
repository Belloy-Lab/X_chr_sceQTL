# Copy the QC'd chrX bfile as TWAS input
bfile_dir <- "../../00.1_X_chr_bfile_QC/chrX_QC_EU.0.G"
out_dir <- "X_bfile23"
cmd <- paste0(
  "plink2 ",
  "--bfile ", bfile_dir, " ",
  "--make-bed ",
  "--out ", out_dir
)

cat(cmd, "\n")
system(cmd)