# Per sex x cell type: covariate file (msex dropped for single-sex settings)
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



for (ct in cell_types) {
  for (sex_status in sex_groups) {
    df_ct_cov_extract_valid_id <- fread(paste0("../../01.2_cov_pheno/",sex_status,"_",ct,".cov.txt"))
    
    if (sex_status != "both") {
      df_ct_cov_extract_valid_id[, msex := NULL]
    }
    
    out_dir <- paste0(sex_status,"/",ct,"/",sex_status,"_",ct,".covar")
    fwrite(df_ct_cov_extract_valid_id,out_dir,sep = '\t',quote = F)
    cat("[Generating covariate file] Cell type:",ct,"| sex_status:",sex_status,"| samplesize:",nrow(df_ct_cov_extract_valid_id),"\n")
    
  }
}