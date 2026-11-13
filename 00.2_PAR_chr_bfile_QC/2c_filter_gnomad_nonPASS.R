# Filter f3 (PAR): remove variants that are non-PASS in gnomAD
rm(list = ls())
library(dplyr)
library(data.table)
# CONFIG
gnomad_prefix <- "path/to/gnomad.genomes.chr"   # expects <prefix><chr>_LCR.tsv
out_file      <- "chrX_QC"
ch            <- "X"
x1 <- fread(paste(gnomad_prefix, ch, "_LCR.tsv", sep = ""))
x2 <- fread(paste(out_file, "_EU.0.E.2.bim", sep = ""))
colnames(x2) <- c("CHR", "posID", "GENPOS", "BP", "ALT", "REF")
x1 <- x1[which(x1$FILTER != "PASS"), ]; x1 <- x1[which(x1$BP %in% x2$BP), ]
x1 <- x1[, p1 := paste(BP, REF, ALT, sep = ":")]
x2 <- x2[, p1 := paste(BP, REF, ALT, sep = ":")]; x2 <- x2[, p2 := paste(BP, ALT, REF, sep = ":")]
# match either allele orientation
mm1 <- inner_join(x = x1, y = x2, by = "p1")
mm2 <- inner_join(x = x1, y = x2, by = c("p1" = "p2"))
mm <- rbind(mm1[, c("posID", "ID", "FILTER")], mm2[, c("posID", "ID", "FILTER")])
mm <- mm[which(!duplicated(mm$posID)), ]
print(dim(mm)[1])
fwrite(as.data.frame(mm[, "posID"]),
       paste(out_file, "_EU.0.E.Gnomad3_nonPASS__remove.txt", sep = ""),
       row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
