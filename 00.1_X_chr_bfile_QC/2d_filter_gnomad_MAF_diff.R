# Filter f4: remove variants whose control MAF differs from gnomAD NFE AF by >= 10%
rm(list = ls())
library(dplyr)
library(data.table)
# CONFIG
gnomad_prefix <- "path/to/gnomad.genomes.chr"   # expects <prefix><chr>_AF_AC_AN.tsv
out_file      <- "chrX_QC"
ch            <- "X"
maf_diff_cut  <- 10                              # |MAF diff| in %
x1 <- fread(paste(gnomad_prefix, ch, "_AF_AC_AN.tsv", sep = ""))
x2 <- fread(paste0(out_file, "_EU.0.E.2.bim.CN.afreq"))
colnames(x2)[1:2] <- c("CHR", "posID")
tt <- strsplit(x2$posID, ':'); x2$BP <- unlist(tt)[seq(2, 4 * length(tt), 4)]
x1 <- x1[which(x1$POS %in% x2$BP), ]
x1 <- x1[, p1 := paste(POS, REF, ALT, sep = ":")]; x1 <- x1[, AF_nfe := as.numeric(AF_nfe)]
x2 <- x2[, p1 := paste(BP, REF, ALT, sep = ":")]; x2 <- x2[, p2 := paste(BP, ALT, REF, sep = ":")]
mm1 <- inner_join(x = x1, y = x2, by = "p1")
mm1$MAF_DIFF <- mm1$ALT_FREQS - mm1$AF_nfe; mm1$MAF_RATIO <- mm1$AF_nfe / mm1$ALT_FREQS
# allele-swapped match: use 1 - gnomAD AF
mm2 <- inner_join(x = x1, y = x2, by = c("p1" = "p2"))
mm2$AF_nfe <- (1 - mm2$AF_nfe)
mm2$MAF_DIFF <- mm2$ALT_FREQS - mm2$AF_nfe; mm2$MAF_RATIO <- mm2$AF_nfe / mm2$ALT_FREQS
# duplicate posID: keep smallest |MAF diff|
mm <- rbind(mm1[, c("ID", "posID", "ALT_FREQS", "AF_nfe", "MAF_DIFF", "MAF_RATIO")],
            mm2[, c("ID", "posID", "ALT_FREQS", "AF_nfe", "MAF_DIFF", "MAF_RATIO")])
mm$abs_MAF_DIFF <- abs(mm$MAF_DIFF)
mm <- arrange(mm, mm$abs_MAF_DIFF)
mm <- mm[which(!duplicated(mm$posID)), ]
mm$ALT_FREQS <- mm$ALT_FREQS * 100; mm$AF_nfe <- mm$AF_nfe * 100; mm$MAF_DIFF <- mm$MAF_DIFF * 100
tr <- mm[which(abs(mm$MAF_DIFF) >= maf_diff_cut), "posID"]
dim(tr)
fwrite(tr,
       paste(out_file, "_EU.0.E.Gnomad3_10%_MAF_diff__remove.txt", sep = ""),
       row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
