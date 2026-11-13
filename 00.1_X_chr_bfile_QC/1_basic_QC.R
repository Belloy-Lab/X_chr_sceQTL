# Basic chrX variant/sample QC: geno/mac, differential missingness, HWE, split PAR, set-hh-missing, mind
rm(list = ls())
library(data.table)
library(plinkbinr)
plink_pathway <- get_plink_exe()   # plink1.9; plink2 called by name
# CONFIG
subj_list_CN <- "NCI_ids.txt"
bfile_dir    <- "chrX_QC_0.0.3"
out_file     <- "chrX_QC"
# 1. geno 0.05 / mac 1; set fam phenotype (control = 1, case = 2)
command <- paste("plink2",
                 "--bfile", bfile_dir,
                 "--allow-no-sex --geno 0.05 --mac 1",
                 "--make-bed --out", paste0(out_file, "_EU.0.A"),
                 sep = " ")
system(command)
fam <- fread(paste(out_file, "_EU.0.A.fam", sep = ""))
control_ids <- fread("NCI_ids.txt", header = FALSE)
any_CI_ids  <- fread("any_CI_ids.txt", header = FALSE)
fam$V6 <- -9
fam$V6[fam$V2 %in% control_ids$V2] <- 1
fam$V6[fam$V2 %in% any_CI_ids$V2]  <- 2
fwrite(fam, paste(out_file, "_EU.0.A.fam", sep = ""), sep = "\t", col.names = FALSE)
# 2. geno 0.05 within each sex and dx group
command <- paste("plink2",
                 "--bfile", paste0(out_file, "_EU.0.A"),
                 "--keep any_CI_ids.txt",
                 "--geno 0.05",
                 "--make-bed --out", paste0(out_file, "_EU.0.A.cases"),
                 sep = " ")
system(command)
command <- paste("plink2",
                 "--bfile", paste0(out_file, "_EU.0.A"),
                 "--keep NCI_ids.txt",
                 "--geno 0.05",
                 "--make-bed --out", paste0(out_file, "_EU.0.A.controls"),
                 sep = " ")
system(command)
command <- paste("plink2",
                 "--bfile", paste0(out_file, "_EU.0.A"),
                 "--keep-males --geno 0.05",
                 "--make-bed --out", paste0(out_file, "_EU.0.A.males"),
                 sep = " ")
system(command)
command <- paste("plink2",
                 "--bfile", paste0(out_file, "_EU.0.A"),
                 "--keep-females --geno 0.05",
                 "--make-bed --out", paste0(out_file, "_EU.0.A.females"),
                 sep = " ")
system(command)
vars1 <- fread(paste0(out_file, "_EU.0.A.males.bim"))
vars2 <- fread(paste0(out_file, "_EU.0.A.females.bim"))
vars  <- intersect(unlist(vars1[, 2]), unlist(vars2[, 2]))
fwrite(as.data.frame(vars), paste0(out_file, "_EU.0.A.sex_specific_5%geno_vars_posID"), col.names = FALSE)
vars1 <- fread(paste(out_file, "_EU.0.A.cases.bim", sep = ""))
vars2 <- fread(paste(out_file, "_EU.0.A.controls.bim", sep = ""))
vars  <- intersect(unlist(vars1[, 2]), unlist(vars2[, 2]))
fwrite(as.data.frame(vars), paste0(out_file, "_EU.0.A.dx_specific_5%geno_vars_posID"), col.names = FALSE)
# 3. Differential missingness (dx: P < 1e-5; sex: P < 1e-10)
command <- paste(plink_pathway,
                 "--bfile", paste(out_file, "_EU.0.A", sep = ""),
                 "--extract", paste0(out_file, "_EU.0.A.dx_specific_5%geno_vars_posID"),
                 "--test-missing --keep-allele-order",
                 "--out", paste(out_file, "_EU.0.A.dx", sep = ""),
                 sep = " ")
system(command)
mis <- fread(paste(out_file, "_EU.0.A.dx.missing", sep = ""), header = TRUE)
mis <- mis[which(mis$P < 1e-5), "SNP"]
fwrite(mis, paste0(out_file, "_EU.0.A.dx_specific_diff_P1e-5_vars_posID"), col.names = FALSE)
command <- paste(plink_pathway,
                 "--bfile", paste0(out_file, "_EU.0.A"),
                 "--keep-allele-order",
                 "--make-bed --out", paste0(out_file, "_EU.0.A.for_sex_missingness"),
                 sep = " ")
system(command)
fam <- fread(paste0(out_file, "_EU.0.A.for_sex_missingness.fam"))
fam$V6 <- fam$V5
fwrite(fam, paste0(out_file, "_EU.0.A.for_sex_missingness.fam"), sep = "\t", col.names = FALSE)
command <- paste(plink_pathway,
                 "--bfile", paste0(out_file, "_EU.0.A.for_sex_missingness"),
                 "--extract", paste0(out_file, "_EU.0.A.sex_specific_5%geno_vars_posID"),
                 "--test-missing --keep-allele-order",
                 "--out", paste0(out_file, "_EU.0.A.sex"),
                 sep = " ")
system(command)
mis <- fread(paste0(out_file, "_EU.0.A.sex.missing"), header = TRUE)
mis <- mis[which(mis$P < 1e-10), "SNP"]
fwrite(mis, paste0(out_file, "_EU.0.A.sex_specific_diff_P1e-10_vars_posID"), col.names = FALSE)
# 4. Apply filters
command <- paste("plink2",
                 "--bfile", paste0(out_file, "_EU.0.A"),
                 "--extract", paste0(out_file, "_EU.0.A.dx_specific_5%geno_vars_posID"),
                 "--make-bed --out", paste0(out_file, "_EU.0.B.1"),
                 sep = " ")
system(command)
command <- paste("plink2",
                 "--bfile", paste0(out_file, "_EU.0.B.1"),
                 "--exclude", paste0(out_file, "_EU.0.A.dx_specific_diff_P1e-5_vars_posID"),
                 "--make-bed --out", paste0(out_file, "_EU.0.B.2"),
                 sep = " ")
system(command)
command <- paste("plink2",
                 "--bfile", paste0(out_file, "_EU.0.B.2"),
                 "--extract", paste0(out_file, "_EU.0.A.sex_specific_5%geno_vars_posID"),
                 "--make-bed --out", paste0(out_file, "_EU.0.B.3"),
                 sep = " ")
system(command)
command <- paste("plink2",
                 "--bfile", paste0(out_file, "_EU.0.B.3"),
                 "--exclude", paste0(out_file, "_EU.0.A.sex_specific_diff_P1e-10_vars_posID"),
                 "--make-bed --out", paste0(out_file, "_EU.0.B.4"),
                 sep = " ")
system(command)
count_vars <- function(prefix) {
  cat(prefix, ":",
      nrow(fread(paste0(prefix, ".bim"), header = FALSE)),
      "variants\n")
}
count_vars(paste0(out_file, "_EU.0.A"))
count_vars(paste0(out_file, "_EU.0.B.1"))
count_vars(paste0(out_file, "_EU.0.B.2"))
count_vars(paste0(out_file, "_EU.0.B.3"))
count_vars(paste0(out_file, "_EU.0.B.4"))
# 5. HWE in controls (P < 1e-5)
command <- paste("plink2",
                 "--bfile", paste0(out_file, "_EU.0.B.4"),
                 "--keep", paste0(subj_list_CN),
                 "--hardy",
                 "--out", paste0(out_file, "_EU.0.C"),
                 sep = " ")
system(command)
tt <- fread(paste(out_file, "_EU.0.C.hardy.x", sep = ""))
fwrite(tt[which(tt$P < 0.00001), "ID"],
       paste0(out_file, "_EU.0.C.hardy.x.var_list"), col.names = FALSE)
command <- paste("plink2",
                 "--bfile", paste0(out_file, "_EU.0.B.4"),
                 "--exclude", paste0(out_file, "_EU.0.C.hardy.x.var_list"),
                 "--make-bed --out", paste0(out_file, "_EU.0.C"),
                 sep = " ")
system(command)
# 6. Split PAR, keep chrX, set male het calls missing
command <- paste("plink2",
                 "--bfile", paste0(out_file, "_EU.0.C"),
                 "--split-par hg38",
                 "--make-bed --out", paste0(out_file, "_EU.0.D.1"),
                 sep = " ")
system(command)
command <- paste("plink2",
                 "--bfile", paste0(out_file, "_EU.0.D.1"),
                 "--chr X",
                 "--make-bed --out", paste0(out_file, "_EU.0.D.2"),
                 sep = " ")
system(command)
command <- paste("plink2",
                 "--bfile", paste0(out_file, "_EU.0.D.2"),
                 "--set-hh-missing",
                 "--make-bed --out", paste0(out_file, "_EU.0.E.1"),
                 sep = " ")
system(command)
# 7. mind 0.02; control allele frequencies
command <- paste("plink2",
                 "--bfile", paste0(out_file, "_EU.0.E.1"),
                 "--mind 0.02",
                 "--make-bed --out", paste0(out_file, "_EU.0.E.2"),
                 sep = " ")
system(command)
command <- paste("plink2",
                 "--bfile", paste0(out_file, "_EU.0.E.2"),
                 "--keep", paste0(subj_list_CN),
                 "--freq",
                 "--out", paste(out_file, "_EU.0.E.2.bim.CN", sep = ""),
                 sep = " ")
system(command)
