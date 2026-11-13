# Format XWAS summary statistics (AD, MSA, PD, LBD x rXCI/eXCI/male/female) for FUSION: SNP A1 A2 Z, IDs aligned to the bfile
rm(list = ls())
library(dplyr)
library(data.table)

bfile_non_PAR <- fread("../7.1_input_bfile/X_bfile23.bim")
bfile_id <- data.frame(ID=c(bfile_non_PAR$V2))
rm(bfile_non_PAR)

# match on either allele orientation
prepare_twas_sumstats <- function(df_GWAS,
                                  bfile_id,
                                  EA_col,
                                  NEA_col,
                                  BP_col,
                                  BETA_col,
                                  SE_col,
                                  chr_prefix = "chrX") {
  
  df <- copy(df_GWAS)
  
  df$Z <- df[[BETA_col]] / df[[SE_col]]
  
  df$id1 <- paste0(chr_prefix, ":", df[[BP_col]], ":", df[[NEA_col]], ":", df[[EA_col]])
  df$id2 <- paste0(chr_prefix, ":", df[[BP_col]], ":", df[[EA_col]], ":", df[[NEA_col]])
  
  df$match <- ifelse(df$id1 %in% bfile_id$ID, "id1",
                     ifelse(df$id2 %in% bfile_id$ID, "id2", NA))
  
  df <- df[!is.na(df$match), ]
  
  df$SNP <- ifelse(df$match == "id1", df$id1, df$id2)
  
  df_out <- df %>%
    transmute(
      SNP = SNP,
      A1  = .data[[EA_col]],
      A2  = .data[[NEA_col]],
      Z   = Z
    ) %>%
    na.omit()
  
  return(df_out)
}



# AD
# rxci
phenotype_name <- "AD_rXCI"
GWAS_path <- "path/to/AD/XWAS_rXCI.txt"

df_GWAS <- fread(GWAS_path)
savename <- paste0(phenotype_name,".txt")
df_TWAS <- prepare_twas_sumstats(
  df_GWAS  = df_GWAS,
  bfile_id = bfile_id,
  EA_col   = "ALLELE1",
  NEA_col  = "ALLELE0",
  BP_col   = "BP",
  BETA_col = "BETA",
  SE_col   = "SE"
)
fwrite(df_TWAS,savename,sep = '\t',quote = F)


# exci
phenotype_name <- "AD_eXCI"
GWAS_path <- "path/to/AD/XWAS_eXCI.txt"

df_GWAS <- fread(GWAS_path)
savename <- paste0(phenotype_name,".txt")
df_TWAS <- prepare_twas_sumstats(
  df_GWAS  = df_GWAS,
  bfile_id = bfile_id,
  EA_col   = "ALLELE1",
  NEA_col  = "ALLELE0",
  BP_col   = "BP",
  BETA_col = "BETA",
  SE_col   = "SE"
)
fwrite(df_TWAS,savename,sep = '\t',quote = F)


# male
phenotype_name <- "AD_male"
GWAS_path <- "path/to/AD/XWAS_male.txt"

df_GWAS <- fread(GWAS_path)
savename <- paste0(phenotype_name,".txt")
df_TWAS <- prepare_twas_sumstats(
  df_GWAS  = df_GWAS,
  bfile_id = bfile_id,
  EA_col   = "ALLELE1",
  NEA_col  = "ALLELE0",
  BP_col   = "BP",
  BETA_col = "BETA",
  SE_col   = "SE"
)
fwrite(df_TWAS,savename,sep = '\t',quote = F)


# female
phenotype_name <- "AD_female"
GWAS_path <- "path/to/AD/XWAS_female.txt"

df_GWAS <- fread(GWAS_path)
savename <- paste0(phenotype_name,".txt")
df_TWAS <- prepare_twas_sumstats(
  df_GWAS  = df_GWAS,
  bfile_id = bfile_id,
  EA_col   = "ALLELE1",
  NEA_col  = "ALLELE0",
  BP_col   = "BP",
  BETA_col = "BETA",
  SE_col   = "SE"
)
fwrite(df_TWAS,savename,sep = '\t',quote = F)










# MSA
# rxci
phenotype_name <- "MSA_rXCI"
GWAS_path <- "path/to/MSA/XWAS_rXCI.txt"

df_GWAS <- fread(GWAS_path)
savename <- paste0(phenotype_name,".txt")
df_TWAS <- prepare_twas_sumstats(
  df_GWAS  = df_GWAS,
  bfile_id = bfile_id,
  EA_col   = "effect_allele",
  NEA_col  = "other_allele",
  BP_col   = "base_pair_location",
  BETA_col = "beta",
  SE_col   = "standard_error"
)
fwrite(df_TWAS,savename,sep = '\t',quote = F)


# exci
phenotype_name <- "MSA_eXCI"
GWAS_path <- "path/to/MSA/XWAS_eXCI.txt"

df_GWAS <- fread(GWAS_path)
savename <- paste0(phenotype_name,".txt")
df_TWAS <- prepare_twas_sumstats(
  df_GWAS  = df_GWAS,
  bfile_id = bfile_id,
  EA_col   = "ALLELE1",
  NEA_col  = "ALLELE0",
  BP_col   = "BP",
  BETA_col = "BETA",
  SE_col   = "SE"
)
fwrite(df_TWAS,savename,sep = '\t',quote = F)


# male
phenotype_name <- "MSA_male"
GWAS_path <- "path/to/MSA/XWAS_male.txt"

df_GWAS <- fread(GWAS_path)
savename <- paste0(phenotype_name,".txt")
df_TWAS <- prepare_twas_sumstats(
  df_GWAS  = df_GWAS,
  bfile_id = bfile_id,
  EA_col   = "effect_allele",
  NEA_col  = "other_allele",
  BP_col   = "base_pair_location",
  BETA_col = "beta",
  SE_col   = "standard_error"
)
fwrite(df_TWAS,savename,sep = '\t',quote = F)


# female
phenotype_name <- "MSA_female"
GWAS_path <- "path/to/MSA/XWAS_female.txt"

df_GWAS <- fread(GWAS_path)
savename <- paste0(phenotype_name,".txt")
df_TWAS <- prepare_twas_sumstats(
  df_GWAS  = df_GWAS,
  bfile_id = bfile_id,
  EA_col   = "effect_allele",
  NEA_col  = "other_allele",
  BP_col   = "base_pair_location",
  BETA_col = "beta",
  SE_col   = "standard_error"
)
fwrite(df_TWAS,savename,sep = '\t',quote = F)






# PD
# rxci
phenotype_name <- "PD_rXCI"
GWAS_path <- "path/to/PD/XWAS_rXCI.txt"

df_GWAS <- fread(GWAS_path)
savename <- paste0(phenotype_name,".txt")
df_TWAS <- prepare_twas_sumstats(
  df_GWAS  = df_GWAS,
  bfile_id = bfile_id,
  EA_col   = "ALLELE1",
  NEA_col  = "ALLELE0",
  BP_col   = "BP",
  BETA_col = "BETA",
  SE_col   = "SE"
)
fwrite(df_TWAS,savename,sep = '\t',quote = F)


# exci
phenotype_name <- "PD_eXCI"
GWAS_path <- "path/to/PD/XWAS_eXCI.txt"

df_GWAS <- fread(GWAS_path)
savename <- paste0(phenotype_name,".txt")
df_TWAS <- prepare_twas_sumstats(
  df_GWAS  = df_GWAS,
  bfile_id = bfile_id,
  EA_col   = "ALLELE1",
  NEA_col  = "ALLELE0",
  BP_col   = "BP",
  BETA_col = "BETA",
  SE_col   = "SE"
)
fwrite(df_TWAS,savename,sep = '\t',quote = F)


# male
phenotype_name <- "PD_male"
GWAS_path <- "path/to/PD/XWAS_male.txt"

df_GWAS <- fread(GWAS_path)
savename <- paste0(phenotype_name,".txt")
df_TWAS <- prepare_twas_sumstats(
  df_GWAS  = df_GWAS,
  bfile_id = bfile_id,
  EA_col   = "ALLELE1",
  NEA_col  = "ALLELE0",
  BP_col   = "BP",
  BETA_col = "BETA",
  SE_col   = "SE"
)
fwrite(df_TWAS,savename,sep = '\t',quote = F)


# female
phenotype_name <- "PD_female"
GWAS_path <- "path/to/PD/XWAS_female.txt"

df_GWAS <- fread(GWAS_path)
savename <- paste0(phenotype_name,".txt")
df_TWAS <- prepare_twas_sumstats(
  df_GWAS  = df_GWAS,
  bfile_id = bfile_id,
  EA_col   = "ALLELE1",
  NEA_col  = "ALLELE0",
  BP_col   = "BP",
  BETA_col = "BETA",
  SE_col   = "SE"
)
fwrite(df_TWAS,savename,sep = '\t',quote = F)










# LBD
# rxci
phenotype_name <- "LBD_rXCI"
GWAS_path <- "path/to/LBD/XWAS_rXCI.txt"

df_GWAS <- fread(GWAS_path)
savename <- paste0(phenotype_name,".txt")
df_TWAS <- prepare_twas_sumstats(
  df_GWAS  = df_GWAS,
  bfile_id = bfile_id,
  EA_col   = "effect_allele",
  NEA_col  = "other_allele",
  BP_col   = "base_pair_location",
  BETA_col = "beta",
  SE_col   = "standard_error"
)
fwrite(df_TWAS,savename,sep = '\t',quote = F)


# exci
phenotype_name <- "LBD_eXCI"
GWAS_path <- "path/to/LBD/XWAS_eXCI.txt"

df_GWAS <- fread(GWAS_path)
savename <- paste0(phenotype_name,".txt")
df_TWAS <- prepare_twas_sumstats(
  df_GWAS  = df_GWAS,
  bfile_id = bfile_id,
  EA_col   = "ALLELE1",
  NEA_col  = "ALLELE0",
  BP_col   = "BP",
  BETA_col = "BETA",
  SE_col   = "SE"
)
fwrite(df_TWAS,savename,sep = '\t',quote = F)


# male
phenotype_name <- "LBD_male"
GWAS_path <- "path/to/LBD/XWAS_male.txt"

df_GWAS <- fread(GWAS_path)
savename <- paste0(phenotype_name,".txt")
df_TWAS <- prepare_twas_sumstats(
  df_GWAS  = df_GWAS,
  bfile_id = bfile_id,
  EA_col   = "effect_allele",
  NEA_col  = "other_allele",
  BP_col   = "base_pair_location",
  BETA_col = "beta",
  SE_col   = "standard_error"
)
fwrite(df_TWAS,savename,sep = '\t',quote = F)


# female
phenotype_name <- "LBD_female"
GWAS_path <- "path/to/LBD/XWAS_female.txt"

df_GWAS <- fread(GWAS_path)
savename <- paste0(phenotype_name,".txt")
df_TWAS <- prepare_twas_sumstats(
  df_GWAS  = df_GWAS,
  bfile_id = bfile_id,
  EA_col   = "effect_allele",
  NEA_col  = "other_allele",
  BP_col   = "base_pair_location",
  BETA_col = "beta",
  SE_col   = "standard_error"
)
fwrite(df_TWAS,savename,sep = '\t',quote = F)