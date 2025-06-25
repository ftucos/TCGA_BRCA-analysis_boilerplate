library(msigdbr)
library(GSVA)
library(data.table)
library(R.utils) # required for reading gzipped files
library(tidyverse)
library(glue)

# Load Data
exp <- fread("data/processed/deduplicated_and_filtered-data_mrna_seq_v2_rsem.tsv.gz") %>%
  column_to_rownames("external_gene_name") %>%
  as.matrix()

metadata <- readRDS("data/processed/TCGA_BRCA-metadata.rds")

msigdb <- msigdbr(species = "Homo sapiens")

signatures <- c("MALTA_CURATED_STEMNESS_MARKERS", "RAMALHO_STEMNESS_UP", "WONG_ADULT_TISSUE_STEM_MODULE", "WONG_EMBRYONIC_STEM_CELL_CORE", "LIM_MAMMARY_STEM_CELL_UP", "MIKKELSEN_PLURIPOTENT_STATE_UP", "ENGELMANN_CANCER_PROGENITORS_UP")

gsNames2signaturesList <- function(signature_name) {
  msigdb %>% filter(gs_name == signature_name) %>% pull("gene_symbol") %>% unique()
}

signatures_list <- signatures %>%
  set_names() %>%
  map(gsNames2signaturesList)

# check how many genes genes are present in gene names of the expression dataset
check_missing_genes <- function(signature_name, expression_matrix) {
  genes <- signatures_list[[signature_name]]
  n_genes <- length(genes)
  genes_in_data <- genes[genes %in% rownames(expression_matrix)]
  n_genes_in_data <- length(genes_in_data)
  missing_genes <- genes[!genes %in% rownames(expression_matrix)]
  pct_genes_in_data <- round(100*n_genes_in_data/n_genes,1) 
  print(glue("\n\n{signature_name}: {pct_genes_in_data}% ({n_genes_in_data}/{n_genes}) genes in expression dataset\n"))
  #if(n_genes_in_data < n_genes) {print("missing genes:", print(missing_genes, collapse = "; "))}
  if(n_genes_in_data < n_genes) {cat(paste0("- missing_genes: ", paste0(missing_genes, collapse = "; "), "\n"))}
}

walk(names(signatures_list), ~check_missing_genes(., expression_matrix = exp))

param <- gsvaParam(exp, geneSets = signatures_list)

GSVA.results <- gsva(param)

GSVA.result.df <- GSVA.results %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SAMPLE_ID")

write_tsv(GSVA.result.df, "data/interim/GSVA/stemness_scores-GSVA.tsv")
