library(data.table)
library(tidyverse)
library(msigdbr)
library(R.utils) # required for reading gzipped files
library(glue)
library(ggstatsplot)
library(ggbeeswarm)
library(ggpubr)
library(ggh4x)

# Load Load Data
exp <- fread("data/processed/deduplicated_and_filtered-data_mrna_seq_v2_rsem.tsv.gz")

metadata <- readRDS("data/processed/TCGA_BRCA-metadata.rds")

# scores <- fread("data/interim/GSVA/stemness_scores-GSVA.tsv")
metadata <- metadata %>%
  left_join(scores)

# list genes of interest to be added to metadata ------------
exp_of_interest <- c("CDK12", "CDK1")
CNA_of_interest <- c("PTEN")
mut_of_interes <- c("AKT1", "PIK3CA", "PTEN")

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

  gene_mutations <- fread("data/raw_data/brca_tcga_pan_can_atlas_2018/data_mutations.txt")
  
  }






# used to plot all the scores at one with faceting variable
metadata.scores_long <-  metadata %>%
  full_join(gene_expression) %>%
  full_join(scores %>% gather(key="Signature", value = "GSVA", -1))

# plot --------------
ggplot(metadata.scores_long, aes(x=log2(ABCC1_EXP+1), y=GSVA, group=Signature)) +
  geom_point(size = 0.3) +
  stat_smooth(method="lm", linewidth =0.7) +
  stat_cor(label.x = 6.5, label.y=-0.6, size = 8/.pt)+
  theme_bw(8) +
  facet_wrap(~Signature) +
  force_panelsizes(rows = unit(4.5, "cm"), cols = unit(5.5, "cm"))

ggsave("results/TCGA_BRCA/TCGA_BRCA-ABCC1_expression_and_stemness_signatures.pdf", width = 20, height = 20, unit = "cm" )

ggplot(metadata.scores_long %>% filter(Pam50 != "NC"), aes(x=log2(ABCC1_EXP+1), y=GSVA, group=Signature)) +
  geom_point(size = 0.3) +
  stat_smooth(method="lm", linewidth =0.7) +
  stat_cor(label.x = 6.5, label.y=-0.6, size = 8/.pt)+
  theme_bw(8) +
  facet_grid(Pam50~Signature) +
  force_panelsizes(rows = unit(4.5, "cm"), cols = unit(5.5, "cm"))

ggsave("results/TCGA_BRCA/TCGA_BRCA-ABCC1_expression_and_stemness_signatures-by_subtype.pdf", width = 60, height = 30, unit = "cm" )


# ABCC1 expression and SRCIN1 CNA ---------
ggstatsplot::ggbetweenstats(metadata.1 %>% mutate(LOG2_ABCC1_EXP = log2(ABCC1_EXP+1)),
                            x=SRCIN1_CNA, y=LOG2_ABCC1_EXP,
                            type = "parametric", bf.message = F, centrality.plotting = F)

metadata.1 %>% filter(Pam50 != "") %>% mutate(Pam50 = str_remove(Pam50, "BRCA_")) %>%
  mutate(SRCIN1_CNA = factor(SRCIN1_CNA, levels = c("Diploid", "Deep/Homozygous Deletion", "Shallow Deletion", "Gain", "Amplification"))) %>%
  mutate(Pam50 = factor(Pam50, levels = c("LumA", "LumB", "Her2", "Basal", "Normal"))) %>%
  lm(ABCC1_EXP ~ SRCIN1_CNA + Pam50, data = .) %>%
    summary()
