# Per cell type Venn: female-biased eGenes vs female-upregulated DEGs vs known eXCI genes
rm(list = ls())
library(tidyverse)
library(ggvenn)
# CONFIG
indep_x_rds   <- "../05_independent_signal_X/df_eqtl_wide_after_second_clump.rds"
input_x_rds   <- "../03_merge_eQTL_res_X/df_eqtl_cis_all.rds"
input_par_rds <- "../03_merge_eQTL_res_PAR/df_eqtl_cis_all.rds"
deg_xlsx      <- "../06.3_sex_DEG/DEG_female_vs_male_X_chr_gene_non_adj_disease.xlsx"
exci_xlsx     <- "df_eXCI_gene_reported.xlsx"

cell_types      <- c("Exc","Inh","Ast","Oli","OPC","Mic","End")
score_threshold <- 1
# female-biased eGenes: sex-het FDR threshold
sex_het_fdr     <- 0.1
# TRUE: known eXCI restricted to this cell type's input genes (plus genes in the two female sets)
known_restrict_to_input <- TRUE
# 1. Per-cell-type input genes (X + PAR)
df_input_X   <- readRDS(input_x_rds)
df_input_PAR <- readRDS(input_par_rds)
df_input <- bind_rows(
  distinct(df_input_X,   cell_type, gene),
  distinct(df_input_PAR, cell_type, gene)
) %>% distinct(cell_type, gene)
# 2. Female-upregulated DEGs
df_female_up_xchr_gene <- readxl::read_excel(deg_xlsx)
df_female_up_xchr_gene$gene <- df_female_up_xchr_gene$updated_Gene_symbol
df_female_up_xchr_gene <- filter(df_female_up_xchr_gene, logFC > 0, adj.P.Val < 0.05)
# 3. Female-biased eGenes (sex-het FDR < sex_het_fdr)
df_eGene_female_biased <- readRDS(indep_x_rds)
df_eGene_female_biased <- filter(df_eGene_female_biased,
                                 P_het_FDR < sex_het_fdr,
                                 sex_het_label_nominal %in% c("Female-dominant (opposite)","Female-biased (concordant)"))
df_eGene_female_biased <- distinct(df_eGene_female_biased, cell_type, gene)
# 4. Known eXCI genes
df_exci <- readxl::read_excel(exci_xlsx)
df_exci <- filter(df_exci, escape_score >= score_threshold)
df_exci <- dplyr::select(df_exci, gene = gene_name)
df_exci <- distinct(df_exci, gene)
gene_known_all <- unique(df_exci$gene)
cat("number of known eXCI genes (score >=", score_threshold, "):", length(gene_known_all), "\n")
# 5. Venn per cell type; overlap CSVs
result_list <- list()
for (ct in cell_types) {
  
  cat("generate venn plot in cell_type:", ct, "\n")
  
  save_name <- paste0("venn3_female_biased_femaleUp_knownEXCI_", ct, ".pdf")
  
  gene_female_biased <- filter(df_eGene_female_biased, cell_type == ct)
  gene_female_biased <- gene_female_biased$gene
  gene_female_up_deg <- filter(df_female_up_xchr_gene, cell_type == ct)
  gene_female_up_deg <- gene_female_up_deg$updated_Gene_symbol
  if (known_restrict_to_input) {
    input_ct   <- df_input$gene[df_input$cell_type == ct]
    universe   <- union(input_ct, c(gene_female_biased, gene_female_up_deg))
    gene_known <- intersect(gene_known_all, universe)
  } else {
    gene_known <- gene_known_all
  }
  cat("   known eXCI in this cell type:", length(gene_known), "\n")
  
  venn_list <- list(
    `female-biased eGenes`     = gene_female_biased,
    `female-upregulated genes` = gene_female_up_deg,
    `known eXCI`               = gene_known
  )
  
  # rows: union of the two female sets; in_known is an annotation column
  all_genes <- unique(c(gene_female_biased, gene_female_up_deg))
  
  df_ct <- tibble(
    cell_type        = ct,
    gene             = all_genes,
    in_female_biased = all_genes %in% gene_female_biased,
    in_female_up     = all_genes %in% gene_female_up_deg,
    in_known         = all_genes %in% gene_known
  ) %>%
    mutate(
      n_sets = in_female_biased + in_female_up + in_known,
      intersection = dplyr::case_when(
        in_female_biased &  in_female_up &  in_known ~ "female_biased & female_up & known",
        in_female_biased &  in_female_up & !in_known ~ "female_biased & female_up",
        in_female_biased & !in_female_up &  in_known ~ "female_biased & known",
        !in_female_biased &  in_female_up &  in_known ~ "female_up & known",
        in_female_biased & !in_female_up & !in_known ~ "female_biased only",
        !in_female_biased &  in_female_up & !in_known ~ "female_up only",
        !in_female_biased & !in_female_up &  in_known ~ "known only",
        TRUE ~ "none"
      )
    )
  
  if (nrow(df_ct) > 0) {
    result_list[[ct]] <- df_ct
  }
  
  p <- ggvenn(
    venn_list,
    fill_color   = c("#33a02c", "red", "#6a3d9a"),
    stroke_size  = 0.8,
    stroke_color = "white",
    set_name_size = 5,
    text_size    = 8,
    show_percentage = FALSE
  ) +
    theme(
      plot.margin = margin(0, 0, 0, 0)
    )
  
  ggsave(
    save_name,
    plot = p,
    width = 6,
    height = 6,
    dpi = 300
  )
}
result_df <- bind_rows(result_list)
# read by 06.5_gene_cards
write.csv(
  result_df,
  "female_biased_female_up_adjP_overlap_genes.csv",
  row.names = FALSE
)
