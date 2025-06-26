library(msigdbr)
library(GSVA)
library(singscore)
library(GSEABase) # required for GeneSetCollection()
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

signatures <- c("HALLMARK_PI3K_AKT_MTOR_SIGNALING")

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

# Compute GSVA -------------------
param <- gsvaParam(exp, geneSets = signatures_list)

GSVA.results <- gsva(param)

GSVA.result.df <- GSVA.results %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SAMPLE_ID")

write_tsv(GSVA.result.df, "data/interim/GSVA/GSVA_scores.tsv")

# Compute singscore --------------

# filter genes with low expression as recommended by singscore
# https://www.bioconductor.org/packages/release/workflows/vignettes/SingscoreAMLMutations/inst/doc/workflow_transcriptional_mut_sig.html
# Check TCGA gene expression scale
median(colSums(exp)/1e6)
# remove genes with less than 10 counts in at least 50% of samples (recondised in case of rare molecular subtypes)
min_counts <- 10
patients_prop <- 0.5
genes_to_keep <- rowSums(exp >= min_counts)/ncol(exp) >= patients_prop
cat(glue("Keeping {sum(genes_to_keep)}/{nrow(exp)} genes with at least {min_counts} counts in at least {patients_prop*100}% of samples\n"))
exp_filtered <- exp[genes_to_keep, ]
log_exp <- log2(exp_filtered + 1)

# convert signatures list to GeneSetCollection required by singscore
signatures_collection <- signatures_list %>%
  names() %>%
  map(~GeneSet(signatures_list[[.]], setName = .)) %>%
  GeneSetCollection()

rankData <- rankGenes(log_exp)
singscore.result <- multiScore(rankData, upSet = signatures_collection)

singscore.result.df <- singscore.result$Scores %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SAMPLE_ID")

write_tsv(singscore.result.df, "data/interim/GSVA/singscore_scores.tsv")
