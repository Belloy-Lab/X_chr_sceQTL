# Per sex x cell type x non-PAR gene: extract cis-window bfile
library(dplyr)
library(data.table)


bfile_dir <- "../7.1_input_bfile/X_bfile23"

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
df_tss <- dplyr::select(df_tss,gene,PAR_type,cis_low,cis_high)

for (ct in cell_types) {
  df_ct_gene_exp <- fread(paste0("../../01.2_cov_pheno/both_",ct,".pheno_chrX.tsv"))
  df_ct_gene_exp$FID <- NULL
  df_ct_gene_exp$IID <- NULL
  ct_input_gene <- data.frame(gene=colnames(df_ct_gene_exp))
  
  ct_input_gene <- left_join(ct_input_gene,df_tss,by='gene')
  ct_input_gene <- dplyr::filter(ct_input_gene,PAR_type=='non_PAR')
  
  cat("Generating gene expression weight files in",ct,"for",nrow(ct_input_gene),"non-PAR genes","\n")
  
  for (sex_status in sex_groups) {
    sample_id_path <- paste0("../../01.1_valid_sample_id/valid_sample_id_",sex_status,"_",ct,".txt")
    
    for (i in 1:nrow(ct_input_gene)) {
      
      gene_name <- ct_input_gene$gene[i]
      cis_low <- ct_input_gene$cis_low[i]
      cis_high <- ct_input_gene$cis_high[i]
      
      cat("Cell type:",ct,"| Gene:",gene_name,"| cis region:",cis_low,"-",cis_high,"\n")
      
      out_dir <- paste0(sex_status,"/",ct,"/",gene_name,"_",sex_status,"_",ct)
      
      cmd <- paste0(
        "plink2 ",
        "--bfile ", bfile_dir, " ",
        "--chr X ",
        "--from-bp ", cis_low, " ",
        "--to-bp ", cis_high, " ",
        "--keep ", sample_id_path, " ",
        "--make-bed ",
        "--out ", out_dir
      )
      
      cat(cmd, "\n")
      system(cmd)
      
      bim_file <- paste0(out_dir, ".bim")
      
      if (file.exists(bim_file)) {
        if (nrow(fread(bim_file)) == 0) {
          file.remove(paste0(out_dir,".bed"))
          file.remove(paste0(out_dir,".bim"))
          file.remove(paste0(out_dir,".fam"))
          next
        }
      }
    }
  }
}