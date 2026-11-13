# Filter f5: Fisher test of allele frequency by sex in controls; remove P < 1e-5
rm(list = ls())
library(data.table)
library(plinkbinr)
plink_pathway <- get_plink_exe()   # plink1.9
# CONFIG
subj_list_CN <- "NCI_ids.txt"
out_file     <- "chrX_QC"
# 1. Phenotype = sex; Fisher test in controls
command <- paste(plink_pathway,
                 "--bfile", paste0(out_file, "_EU.0.E.2"),
                 "--keep-allele-order",
                 "--make-bed --out", paste0(out_file, "_EU.0.E.2.sex_as_pheno"),
                 sep = " ")
system(command)
fam <- fread(paste0(out_file, "_EU.0.E.2.sex_as_pheno.fam"))
fam$V6 <- fam$V5
fwrite(fam, paste0(out_file, "_EU.0.E.2.sex_as_pheno.fam"), sep = "\t", col.names = FALSE)
command <- paste(plink_pathway,
                 "--bfile", paste(out_file, "_EU.0.E.2.sex_as_pheno", sep = ""),
                 "--keep", subj_list_CN,
                 "--fisher --keep-allele-order",
                 "--out", paste(out_file, "_EU.0.E.2.sex_as_pheno.sex_fisher", sep = ""),
                 sep = " ")
system(command)
# 2. Write P-value table and remove-list
df  <- fread(paste(out_file, "_EU.0.E.2.sex_as_pheno.sex_fisher.assoc.fisher", sep = ""))
dft <- df[, c("SNP", "P")]; colnames(dft) <- c("posID", "fish.test.P")
fwrite(dft, paste(out_file, "_EU.0.E_cn_fish_test_sex_P.txt", sep = ""),
       col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
dim(((dft[which(dft$fish.test.P < 1e-5), "posID"])))
fwrite(((dft[which(dft$fish.test.P < 1e-5), "posID"])),
       paste(out_file, "_EU.0.E_cn_fish_test_sex__remove.txt", sep = ""),
       col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")