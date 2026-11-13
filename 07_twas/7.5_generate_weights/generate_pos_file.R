# FUSION .pos file per sex x cell type from the .wgt.RDat files; P0/P1 = gene start/end
library(data.table)
# CONFIG
base_dir   <- "weights"
sex_groups <- c("both", "female", "male")
ref_rds    <- "../../03_merge_eQTL_res_X/TSS_all_chr.rds"
# 1. Gene coordinates
ref <- as.data.table(readRDS(ref_rds))
gene_pos <- unique(ref[, .(
  gene   = gene_symbol,
  P0     = pmin(gene_start, gene_end),
  P1     = pmax(gene_start, gene_end)
)])
gene_pos <- gene_pos[!duplicated(gene)]
# 2. .pos table per sex x cell type
for (sex_i in sex_groups) {
  
  sex_dir <- file.path(base_dir, sex_i)
  
  cell_types <- dir(sex_dir)
  
  for (cell_i in cell_types) {
    
    cat("Processing:", sex_i, cell_i, "\n")
    
    cell_dir <- file.path(sex_dir, cell_i)
    
    rdat_files <- list.files(cell_dir, pattern = "\\.wgt\\.RDat$")
    
    if (length(rdat_files) == 0) {
      cat("  No weight files found\n")
      next
    }
    
    gene_ids <- sub("\\.wgt\\.RDat$", "", rdat_files)
    
    pos_table <- data.table(
      WGT   = rdat_files,
      ID    = gene_ids,
      CHR   = 23,   # FUSION: 23 = X
      PANEL = cell_i
    )
    pos_table <- gene_pos[pos_table, on = c(gene = "ID")]
    setnames(pos_table, "gene", "ID")
    pos_table <- pos_table[, .(WGT, ID, CHR, P0, P1, PANEL)]
    
    # drop genes without coordinates
    missing <- pos_table[is.na(P0) | is.na(P1)]
    if (nrow(missing) > 0) {
      cat("  WARNING: no coordinates for", nrow(missing), "gene(s), dropped: ",
          paste(missing$ID, collapse = ", "), "\n")
      pos_table <- pos_table[!(is.na(P0) | is.na(P1))]
    }
    
    out_file <- file.path(cell_dir, paste0(cell_i, ".pos"))
    fwrite(pos_table, file = out_file, sep = "\t", quote = FALSE)
    
    cat("  Saved:", out_file, "\n")
    cat("  Gene number:", nrow(pos_table), "\n\n")
  }
}