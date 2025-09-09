library(data.table)
library(tidyverse)
library(msigdbr)
library(R.utils) # required for reading gzipped files
library(glue)
library(ggstatsplot)
library(ggbeeswarm)
library(ggpubr) 
library(ggh4x) # to force ggplot panel sizes
library(survival)
library(survminer)
library(broom) # for tidy coxph output

# Load Load Data
log_transform_exp <- TRUE
exp <- fread("data/processed/deduplicated_and_filtered-data_mrna_seq_v2_rsem.tsv.gz")

if(log_transform_exp) {
  exp <- exp %>%
    column_to_rownames("external_gene_name") %>%
    as.matrix() %>%
    {log2(. + 1)} %>%
    as.data.frame() %>%
    rownames_to_column("external_gene_name")
}

metadata <- readRDS("data/processed/TCGA_BRCA-metadata.rds")

scores <- fread("data/interim/GSVA/singscore_scores.tsv")

metadata <- metadata %>%
  left_join(scores)

# list genes of interest to be added to metadata ------------
exp_of_interest <- c("CDK12", "CDK1")
CNA_of_interest <- c("PTEN")
mut_of_interes <- c("AKT1", "PIK3CA", "PTEN")

# merge expression of interest
if(length(exp_of_interest) > 0) {
  gene_expression <- exp %>%
    filter(external_gene_name %in% exp_of_interest) %>%
    column_to_rownames("external_gene_name") %>%
    t() %>%
    as.data.frame() %>%
    rename_all(~paste0(., "_EXP")) %>%
    rownames_to_column("SAMPLE_ID")
  
  metadata <- metadata %>%
    left_join(gene_expression)
}

if(length(CNA_of_interest) > 0) {
  gene_CNA <- fread("data/raw_data/brca_tcga_pan_can_atlas_2018/data_cna.txt") %>%
    select(-Entrez_Gene_Id) %>%
    filter(Hugo_Symbol %in% CNA_of_interest) %>%
    column_to_rownames("Hugo_Symbol") %>%
    as.matrix() %>% t() %>% as.data.frame() %>%
    rename_all(~paste0(., "_CNA")) %>%
    mutate_all(as.character) %>%
    mutate_all(~recode(.,
                       "-2" = "Deep/Homozygous Deletion",
                       "-1" = "Shallow Deletion",
                       "0" = "Diploid",
                       "1" = "Gain",
                       "2" = "Amplification"
    ) %>% factor(levels = c("Deep/Homozygous Deletion", "Shallow Deletion", "Diploid", "Gain", "Amplification"))) %>%
    rownames_to_column("SAMPLE_ID")
  
  metadata <- metadata %>%
    left_join(gene_CNA)
}

if(length(mut_of_interes) > 0) {
  # Regex string for PIK3CA variants of interest according to OncoKB (2026-06-01)
  PIK3CA_variants <- "(C420R|E542K|E545[ADGKQ]|G1049R|H1047[LRY]|M1043[IV]|N345K|Q546[EKPR]|R88Q)"
  # Regex string for PIK3CA variants of interest according to ESCAT 2019
  # PIK3CA_variants <- "H1047R|H1047L|E542K|E545K|E545A"
  gene_mutations_raw <- fread("data/raw_data/brca_tcga_pan_can_atlas_2018/data_mutations.txt")
  gene_mutations <- gene_mutations_raw %>%
    filter(Hugo_Symbol %in% mut_of_interes) %>%
    filter(Hugo_Symbol == "AKT1" & str_detect(HGVSp_Short, "E17K") | # excluded  p.*2* (1), p.L52R (1), p.V416= (1), p.V4L (1)
           Hugo_Symbol == "PIK3CA" & str_detect(HGVSp_Short, PIK3CA_variants) |
           Hugo_Symbol == "PTEN" & str_detect(HGVSp_Short, "R130*|Q171*|G129E")) %>%
    select(SAMPLE_ID = Tumor_Sample_Barcode)
  
  }

# used to plot all the scores at one with faceting variable
metadata.scores_long <-  metadata %>%
  full_join(gene_expression) %>%
  full_join(scores %>% gather(key="Signature", value = "GSVA", -1))

# plot --------------
ggplot(metadata.scores_long, aes(x=CDK12_EXP, y=GSVA, group=Signature)) +
  geom_point(size = 0.3) +
  stat_smooth(method="lm", linewidth =0.7) +
  stat_cor(label.x = 11, label.y=0.1, size = 8/.pt)+
  theme_bw(8) +
  facet_wrap(~Signature) +
  force_panelsizes(rows = unit(4.5, "cm"), cols = unit(5.5, "cm"))

ggsave("results/TCGA_BRCA/TCGA_BRCA-ABCC1_expression_and_stemness_signatures.pdf", width = 20, height = 20, unit = "cm" )

ggplot(metadata.scores_long %>% filter(Pam50 != "NC"), aes(x=CDK12_EXP, y=GSVA, group=Signature)) +
  geom_point(size = 0.3) +
  stat_smooth(method="lm", linewidth =0.7) +
  stat_cor(label.x = 11, label.y=0.1, size = 8/.pt)+
  theme_bw(8) +
  facet_grid(Pam50~Signature) +
  force_panelsizes(rows = unit(4.5, "cm"), cols = unit(5.5, "cm"))

ggsave("results/TCGA_BRCA/TCGA_BRCA-ABCC1_expression_and_stemness_signatures-by_subtype.pdf", width = 60, height = 30, unit = "cm" )

# explore non-linear relationship with survival --------------------------------
cox_res <- coxph(Surv(OS_years, OS) ~ pspline(CDK12_EXP), data = metadata)
cox_res %>% termplot()

# explore multiple survival cutoff and thresholds ------------------------------

source("R/extract_cox_stats.R")
endpoint <- c("OS", "DSS", "DFS")
variable <- c("CDK12_EXP", "CDK1_EXP")
quantiles <- list(0.1, 0.15, 0.25, 0.33, 0.5, 0.66, 0.75, 0.85, 0.9, "best")



# grid of model parameters to test
test_grid <- expand_grid(endpoint = endpoint,
                         variable = variable,
                         selected_quantile = quantiles)

# specify x as data for creating a partial function
map_funciton <- partial(extract_cox_stats, data = metadata)
res <- pmap_dfr(test_grid, map_funciton) %>%
  ungroup() %>%
  # highlight best quantile for each variable/endpoint combination
  mutate(best = ifelse(p.value == min(p.value), TRUE, FALSE), .by = c(variable, endpoint))

# plot thresholds results
res %>%
  # fill missing values for eventual vest quantiles to plot them as missing
  complete(variable, endpoint, quantile) %>%
  # add pvalue label
  mutate(signif = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01 ~ "**",
    p.value < 0.05 ~ "*",
    TRUE ~ "")) %>%
ggplot(aes(y=as.character(quantile), x=endpoint, fill = HR)) +
  geom_tile() +
  geom_text(aes(label = signif)) +
  facet_wrap(~variable) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey90") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  ylab("Quantile") +
  theme_bw() +
  theme(panel.grid = element_blank())
