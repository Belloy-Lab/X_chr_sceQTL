# Per sex x cell type x non-PAR gene: phenotype file (FID, IID, PHENO)
library(dplyr)
library(data.table)

sex_groups <- c("both", "female", "male")
cell_types <- c("Ast","End","Exc","Inh","Mic","Oli","OPC")

for (sex in sex_groups) {
  
  if (!dir.exists(sex)) dir.create(sex)
  
  for (ct in cell_types) {
    
    path <- file.path(sex, ct)
    
    if (!dir.exists(path)) {
      dir.create(path)
    }
  }
}



df_tss <- readxl::read_excel("../../03_merge_eQTL_res_X/tss_gene_X_chr.xlsx")
df_tss <- select(df_tss,gene,PAR_type,cis_low,cis_high)
df_tss_non_PAR <- filter(df_tss, PAR_type=='non_PAR')

for (ct in cell_types) {
  for (sex_status in sex_groups) {
    df_ct_gene_exp_extract_valid_id <- fread(paste0("../../01.2_cov_pheno/",sex_status,"_",ct,".pheno_chrX.tsv"))
    df_ct_gene_exp_extract_valid_id <- select(df_ct_gene_exp_extract_valid_id,
                                              FID,IID,any_of(df_tss_non_PAR$gene))
    ct_input_gene <- setdiff(colnames(df_ct_gene_exp_extract_valid_id), c("FID", "IID"))
    
    for (gene_name in ct_input_gene) {
      
      df_pheno <- select(df_ct_gene_exp_extract_valid_id,
                         FID,IID,all_of(gene_name))
      colnames(df_pheno)[3] <- "PHENO"
      
      cat("[Generating phenotype file] Cell type:",ct,"| Gene:",gene_name,"| sex_status:",sex_status,"| samplesize:",nrow(df_pheno),"\n")
      
      out_dir <- paste0(sex_status,"/",ct,"/",gene_name,"_",sex_status,"_",ct,".pheno")
      
      fwrite(df_pheno,out_dir,sep = '\t',quote = F)
      
      
    }
  }
}