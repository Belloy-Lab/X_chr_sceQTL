# Filter f3 (PAR): remove variants in gnomAD low-complexity regions (LCR)
library(data.table)
# CONFIG
gnomad_prefix <- "path/to/gnomad.genomes.chr"   # expects <prefix><chr>_LCR.tsv
out_file      <- "chrX_QC"
ch            <- "X"
x1 <- fread(paste(gnomad_prefix, ch, "_LCR.tsv", sep = ""))
x2 <- fread(paste(out_file, "_EU.0.E.2.bim", sep = ""))
colnames(x2) <- c("CHR", "posID", "GENPOS", "BP", "ALT", "REF")
x3 <- x2[which(x2$BP %in% x1[which(x1$LCR == 1), BP]), ]
print(dim(x3)[1])
fwrite(as.data.frame(x3[, "posID"]),
       paste(out_file, "_EU.0.E.Gnomad3_LCR__remove.txt", sep = ""),
       row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
