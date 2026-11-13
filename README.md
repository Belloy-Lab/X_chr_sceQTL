# Brain Cell-Type-Specific X-Chromosome Gene Regulation Reveals Sex-Biased Candidate Genes in Neurodegeneration

## Background

Sex differences in neurodegenerative diseases may partly reflect X-chromosome dosage and cell-type-specific regulation, yet X-linked gene regulation in the human brain remains poorly resolved. We mapped X-chromosome cis-expression quantitative trait loci (cis-eQTLs) in the dorsolateral prefrontal cortex of 423 ROSMAP donors across seven brain cell types, integrating:

- An X-chromosome-tailored genotype quality-control pipeline covering both the pseudoautosomal (PAR) and non-pseudoautosomal (non-PAR) regions
- Dosage-aware cis-eQTL mapping: three settings for PAR variants (combined-sex, female-only, male-only) and four settings for non-PAR variants (random X-chromosome inactivation [rXCI], escape from X-chromosome inactivation [eXCI], female-only, male-only)
- Characterization of sex-heterogeneous eQTL effects and sex-biased gene expression
- Integration with X-chromosome-wide association studies (XWAS) of Alzheimer's disease (AD), Parkinson's disease (PD), Lewy body dementia (LBD), and multiple system atrophy (MSA) through cell-type-specific transcriptome-wide association (cTWAS), summary-data-based Mendelian randomization (SMR), and colocalization (COLOC)

This repository provides the code for each of these analyses.

## Repository structure

```
├── 00.1_X_chr_bfile_QC/                        # Non-PAR genotype QC (PLINK)
├── 00.2_PAR_chr_bfile_QC/                      # PAR genotype QC (PLINK)
├── 00.3_rsid_annotation/                       # rsID (dbSNP) and nearest-gene (GENCODE) annotation
├── 01.1_valid_sample_id/                       # Sample lists per cell type and sex stratum
├── 01.2_cov_pheno/                             # Expression normalization (TMM/voom, ComBat, inverse-normal) and covariates (ePCs, gPCs)
├── 02_eQTL_main/                               # cis-eQTL mapping with PLINK 2 (--xchr-model 1/2)
├── 03_merge_eQTL_res_PAR/                      # cis-window extraction and merging, PAR
├── 03_merge_eQTL_res_X/                        # cis-window extraction and merging, non-PAR
├── 04_significant_eQTL_PAR/                    # Two-step FDR and LD clumping, PAR
├── 04_significant_eQTL_X/                      # Two-step FDR and LD clumping, non-PAR
├── 05_independent_signal_PAR/                  # Independent signals across settings, PAR
├── 05_independent_signal_X/                    # Independent signals, sex-heterogeneity test, non-PAR
├── 06.1_manhattan_plot/                        # Manhattan plot of cis-eQTLs
├── 06.2_visualization_distribution_X_PAR_eGENE/ # eGene distribution across cell types
├── 06.3_sex_DEG/                               # Sex-biased differential expression (limma-voom)
├── 06.4_venn_plot/                             # Overlap of female-biased genes with known eXCI genes
├── 06.5_gene_cards_AD/                         # GeneCards AD relevance score analysis
├── 06.5_gene_cards_PD/                         # GeneCards PD relevance score analysis
├── 06.6_LocusZoom_sex_het_example/             # LocusZoom plots for sex-heterogeneous eQTL examples
├── 06.7_ggdensity/                             # MAF distribution of independent lead eQTL variants
├── 07_twas/                                    # cTWAS: FUSION weights, XWAS harmonization, association tests
├── 08_smr/                                     # SMR: BESD construction, XWAS harmonization, SMR/HEIDI
├── 09_coloc/                                   # COLOC (coloc.abf) between cis-eQTL and XWAS signals
├── 10_combine_twas_smr_coloc/                  # Integration of cTWAS, SMR, and COLOC for gene prioritization; LocusZoom plots of prioritized genes
├── README.md
└── LICENSE
```

Scripts are numbered in execution order. Paths to controlled-access input data are placeholders (`path/to/...`) and must be edited in the `CONFIG` block at the top of each script. Job-submission scripts (`submit_all.sh`) are written for an LSF cluster and need to be adapted to other schedulers.

## Data availability

### eQTL analyses

Cell-type-specific X-chromosome cis-eQTL summary statistics generated in this study are available through the GeneXplore browser and the Zenodo repository below. Individual-level genotype and single-nucleus RNA-seq data are controlled-access and are not hosted in this repository.

| Datasource | Link |
|---|---|
| Cell-type-specific X-chromosome cis-eQTL summary statistics (browser) | https://genexplore.wustl.edu/ |
| Cell-type-specific X-chromosome cis-eQTL summary statistics (Zenodo) | https://doi.org/10.5281/zenodo.23989152 |
| ROSMAP whole-genome sequencing (AD Knowledge Portal, Synapse; upon application) | https://www.synapse.org/Synapse:syn10901595 |
| ROSMAP single-nucleus RNA-seq and pseudobulk expression (Fujita et al. 2024; upon application) | https://www.synapse.org/Synapse:syn52335807 |

### XWAS summary statistics

| Trait | Study | Link |
|---|---|---|
| Alzheimer's disease | Belloy et al. 2024; Cook et al. 2026 | https://github.com/Belloy-Lab/AD_XWAS_Pleiotropy |
| Parkinson's disease | Le Guen et al. 2021 | https://www.ebi.ac.uk/gwas/studies/ (GCST90104085–GCST90104087) |
| Lewy body dementia | Bayram et al. 2024 | https://www.ebi.ac.uk/gwas/studies/ (GCST90502821–GCST90502823) |
| Multiple system atrophy | Ray et al. 2026 | https://www.ebi.ac.uk/gwas/studies/ (GCST90558035–GCST90558037) |

### Reference resources

| Resource | Link |
|---|---|
| gnomAD v3 (variant QC) | https://gnomad.broadinstitute.org/ |
| GENCODE v48 basic annotation | https://www.gencodegenes.org/human/ |
| dbSNP (rsID annotation) | https://www.ncbi.nlm.nih.gov/snp/ |
| GeneCards (AD and PD relevance scores) | https://www.genecards.org/ |
| Known eXCI genes (curated from 12 studies; Supplementary Table S6) | `06.4_venn_plot/df_eXCI_gene_reported.xlsx` |

### Software

| Tool | Link |
|---|---|
| PLINK 1.9 / 2.0 | https://www.cog-genomics.org/plink/ |
| FUSION | http://gusevlab.org/projects/fusion/ |
| SMR v1.3.1 (binary not included; place `smr` in `08_smr/8.1_data_clean_eQTL/` and `08_smr/8.3_smr_analysis/`) | https://yanglab.westlake.edu.cn/software/smr/ |
| coloc (R) | https://chr1swallace.github.io/coloc/ |
| ieugwasr / plinkbinr (R; LD clumping) | https://mrcieu.github.io/ieugwasr/ |
| locuszoomr (R; locus plots) | https://github.com/myles-lewis/locuszoomr |

## Citation

If you use these scripts, please cite our paper:

Zeng Y, Cook N, Fujita M, Johnson M, Wingo A, Wingo T, Pottier C, Cruchaga C, Menon V, Belloy ME. Brain Cell-Type-Specific X-Chromosome Gene Regulation Reveals Sex-Biased Candidate Genes in Neurodegeneration. *[preprint DOI to be added]* (2026)

## Contact

Michael E. Belloy – belloy@wustl.edu

## License (MIT)

Copyright (c) 2026 Michael Belloy

Permission is hereby granted, free of charge, to any person obtaining a copy of these softwares, scripts, pipelines, and associated documentation files (the "Code"), to deal in the Code without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Code, and to permit persons to whom the Code is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Code.

THE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE CODE OR THE USE OR OTHER DEALINGS IN THE CODE.
