library(genefu)
library(data.table)
library(tidyverse)
library(skimr)
library(readxl)
library(ggbeeswarm)

set.seed(17)

# Preprocess RNA seq data -----------------------------
exp_raw <- fread("data/raw_data/brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt")

dup_symbol <- exp_raw %>% group_by(Hugo_Symbol) %>% summarize(n = n()) %>% filter(n > 1)
if (nrow(dup_symbol) > 0) {
  cat("Duplicated gene symbols found:\n")
  dup_symbol %>% 
    mutate(symbol = paste0(Hugo_Symbol, " (", n, ")")) %>%
    pull(symbol) %>%
    cat(sep = ", ")
} else {
  cat("No duplicated gene symbols found.\n")
}

# remove duplicated genes, keeping the one with the highest mean expression
exp <- exp_raw %>%
  # add a column with mean expression for each gene
  cbind(.,
        select(., -Hugo_Symbol, -Entrez_Gene_Id) %>% 
          rowMeans() %>%
          as.data.frame()  %>%
          rename("mean_counts" = ".")
  ) %>%
  # replace with Entrez gene Id, genes with missing symbol
  mutate(Hugo_Symbol = ifelse(Hugo_Symbol == "", Entrez_Gene_Id, Hugo_Symbol)) %>%
  select(-Entrez_Gene_Id) %>%
  # remove eventual remaining genes with missing Symbol/Entrez_Gene_Id
  filter(Hugo_Symbol != "") %>%
  # for each gene, keep only the row with the highest median
  group_by(Hugo_Symbol) %>%
  slice_max(mean_counts, n = 1, with_ties = FALSE) %>%
  # drop the helper column
  select(-mean_counts) %>%
  column_to_rownames("Hugo_Symbol") %>%
  as.matrix()

dup_entrez <- exp_raw %>% group_by(Entrez_Gene_Id) %>% summarize(n = n()) %>% filter(n > 1)
if (nrow(dup_entrez) > 0) {
  cat("Duplicated Entrez ID found:\n")
  dup_entrez %>% 
    mutate(Entrez_Gene_Id = paste0(Entrez_Gene_Id, " (", n, ")")) %>%
    pull(Entrez_Gene_Id) %>%
    cat(sep = ", ")
} else {
  cat("No duplicated Entrez ID found.\n")
}


# remove duplicated genes, keeping the one with the highest mean expression
exp_entrez <- exp_raw %>%
  # add a column with mean expression for each gene
  cbind(.,
        select(., -Hugo_Symbol, -Entrez_Gene_Id) %>% 
          rowMeans() %>%
          as.data.frame()  %>%
          rename("mean_counts" = ".")
  ) %>%
  select(-Hugo_Symbol) %>%
  # remove eventual remaining genes with missing Symbol/Entrez_Gene_Id
  filter(Entrez_Gene_Id != "") %>%
  # for each gene, keep only the row with the highest median
  group_by(Entrez_Gene_Id) %>%
  slice_max(mean_counts, n = 1, with_ties = FALSE) %>%
  # drop the helper column
  select(-mean_counts) %>%
  column_to_rownames("Entrez_Gene_Id") %>%
  as.matrix()


# Preprocess clinical ------------------
patient <- fread("data/raw_data/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt", skip = 4, na.strings = c("")) %>%
  filter(PATIENT_ID %in% (colnames(exp) %>% str_remove("-01$")))

sample <- fread("data/raw_data/brca_tcga_pan_can_atlas_2018/data_clinical_sample.txt", skip = 4, na.strings = c("")) %>%
  filter(SAMPLE_ID %in% colnames(exp))

patient_legacy <- fread("data/raw_data/brca_tcga/data_clinical_patient.txt", skip = 4, na.strings = c("", "[Not Available]")) %>%
  filter(PATIENT_ID %in% (colnames(exp) %>% str_remove("-01$")))

sample_legacy <- fread("data/raw_data/brca_tcga/data_clinical_sample.txt", skip = 4, na.strings = c("", "[Not Available]")) %>%
  filter(SAMPLE_ID %in% colnames(exp))

# curated histotypes # from A Thennavan - 2021 --------------------
# re-annotation of papillary neoplasm
table_S1 <- read_excel("data/raw_data/TCGA_third_party_data/Thennavan_et_al-2021-Cell_Genomics/Supplementary Tables S1 S4.xlsx", sheet = "TableS1", skip = 1) %>%
  select(CLID, "2019 Re-annotation")
# re-annotation of metaplastic carcinomas
table_S4 <- read_excel("data/raw_data/TCGA_third_party_data/Thennavan_et_al-2021-Cell_Genomics/Supplementary Tables S1 S4.xlsx", sheet = "TableS4", skip = 1) %>%
  select(CLID, "2019 Re-annotation")

histology_curated <- read_excel("data/raw_data/TCGA_third_party_data/Thennavan_et_al-2021-Cell_Genomics/Supplementary Data S1.xlsx", na = c("-", "na", "[Not Applicable]", "N/A", "[Not Evaluated]")) %>%
  left_join(rbind(table_S1, table_S4)) %>%
  mutate(Histotype_Thennavan_2021 = case_when(
    str_detect(`2019 Re-annotation`, "Invasive ductal carcinoma") ~ "Invasive ductal carcinoma",
    str_detect(`2019 Re-annotation`, "Encapsulated papillary carcinoma") ~ "Encapsulated papillary carcinoma",
    str_detect(`2019 Re-annotation`, "Invasive solid papillary carcinoma") ~ "Invasive solid papillary carcinoma",
    str_detect(`2019 Re-annotation`, "Solid papillary carcinoma") ~ "Solid papillary carcinoma",
    `2019 Re-annotation` == "DCIS, cribriform and micropapillary; possibly focus of IDC" ~ "DCIS",
    TRUE ~ `2016 Histology Annotations`
  )) %>%
  # fix erroneous TCGA-BH-A0WA-01A classification
  mutate(Histotype_Thennavan_2021 = ifelse(CLID == "TCGA-BH-A0WA-01A", "Invasive lobular carcinoma", Histotype_Thennavan_2021))
  
# curated clinical informations 
# same survival times and outcomes, used only for menopausal state
# from Liu et Al 2018
clinical_curated <- read_excel("data/raw_data/TCGA_third_party_data/Liu_et_al-2018-Cell/mmc1.xlsx",
                               guess_max = 2000, na = c("[Not Available]", "[Unknown]", "-", "[Not Applicable]", "[Not Evaluated]")) %>%
  filter(type == "BRCA") %>%
  mutate(menopause_status = recode(menopause_status,
                                   "Indeterminate (neither Pre or Postmenopausal)" = "Indeterminate/Unknown",
                                   "Peri (6-12 months since last menstrual period)" = "Peri",
                                   "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)" = "Post",
                                   "Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)" = "Pre",
                                   .missing = "Indeterminate/Unknown") %>%
           factor(levels = c( "Post", "Peri", "Pre", "Indeterminate/Unknown"))) %>%
  select(PATIENT_ID = bcr_patient_barcode, menopause_status)

# from 1082 patients, removed 1 phylloides, and 2 solid/papillary carcinoma, 12 males --> 1067
metadata <- full_join(
  patient,
  sample %>% select(
    SAMPLE_ID, PATIENT_ID, CANCER_TYPE_DETAILED, TUMOR_TYPE,
    ANEUPLOIDY_SCORE, TMB_NONSYNONYMOUS, TBL_SCORE, MSI_SCORE_MANTIS, MSI_SENSOR_SCORE)
  ) %>%
  # reorganize some columns and rename others
  mutate(
         pT = str_extract(PATH_T_STAGE, "T[0-4X]") %>% factor(levels = c("T1", "T2", "T3", "T4", "TX")),
         pN = str_extract(PATH_N_STAGE, "N[0-3X]") %>% factor(levels = c("N0", "N1", "N2", "N3", "NX")),
         pM = str_extract(PATH_M_STAGE, "M[0-1X]") %>% factor(levels = c("M0", "M1", "MX")), # AJCC now recommends reporting MX as cM0,
         AJCC_Stage = str_extract(AJCC_PATHOLOGIC_TUMOR_STAGE, "(?<=STAGE )[IVX]+") %>% factor(levels = c("I", "II", "III", "IV", "X")),
         OS        = str_extract(OS_STATUS, "^[0-1]") %>% as.numeric(),
         OS_years  = OS_MONTHS/12,
         DSS       = str_extract(DSS_STATUS, "^[0-1]") %>% as.numeric(),
         DSS_years = DSS_MONTHS/12,
         DFS       = str_extract(DFS_STATUS, "^[0-1]") %>% as.numeric(),
         DFS_years = DFS_MONTHS/12,
         PFS       = str_extract(PFS_STATUS, "^[0-1]") %>% as.numeric(),
         PFS_years = PFS_MONTHS/12,
         # more convenient scale to interpret HR in multivariable analysis
         Age_10_years_increase = AGE/10,
         ) %>%
  # add histological subtype according to Firehose Legacy release (same as GDC release)
  left_join(sample_legacy %>%
              select(SAMPLE_ID, CANCER_TYPE_DETAILED), by = "SAMPLE_ID", suffix = c("_PANCANCERATLAS", "_LEGACY")) %>%
  mutate(
    # cases manually reviewed (original report and pthology slide when available)
    Histotype = case_when(
      SAMPLE_ID %in% c("TCGA-A7-A5ZV-01") ~ "Breast Invasive Ductal Carcinoma", # actually annotated as "Breast Invasive Carcinoma No Special Type", no AR expression by RPPA and RNAseq not supporting apocrine differentation reported by Thennavan_2021
      SAMPLE_ID %in% c("TCGA-A1-A0SG-01", "TCGA-EW-A2FV-01") ~ "Breast Invasive Micropapillary Carcinoma",
      SAMPLE_ID %in% c("TCGA-A1-A0SK-01") ~ "Neuroendocrine Carcinoma",
      SAMPLE_ID %in% c("TCGA-D8-A1XS-01") ~ "Breast Invasive Ductal Carcinoma",
      SAMPLE_ID %in% c("TCGA-AO-A03U-01") ~ "Breast Secretory Carcinoma",
      CANCER_TYPE_DETAILED_PANCANCERATLAS == CANCER_TYPE_DETAILED_LEGACY ~ CANCER_TYPE_DETAILED_LEGACY,
      CANCER_TYPE_DETAILED_LEGACY == "Paget Disease of the Nipple" ~ CANCER_TYPE_DETAILED_PANCANCERATLAS,
      CANCER_TYPE_DETAILED_LEGACY %in% c("Solid Papillary Carcinoma of the Breast", "Adenoid Cystic Breast Cancer", "Basal Cell Carcinoma", "Metaplastic Breast Cancer", "Malignant Phyllodes Tumor of the Breast" , "Breast Mixed Ductal and Lobular Carcinoma", "Breast Mixed Ductal and Lobular Carcinoma") ~ CANCER_TYPE_DETAILED_LEGACY,
      CANCER_TYPE_DETAILED_PANCANCERATLAS == "Breast Invasive Lobular Carcinoma" & CANCER_TYPE_DETAILED_LEGACY == "Breast Invasive Ductal Carcinoma" ~ CANCER_TYPE_DETAILED_PANCANCERATLAS,
      CANCER_TYPE_DETAILED_LEGACY == "Breast Invasive Lobular Carcinoma" ~ CANCER_TYPE_DETAILED_LEGACY,
      CANCER_TYPE_DETAILED_LEGACY == "Breast Invasive Mixed Mucinous Carcinoma" ~ CANCER_TYPE_DETAILED_LEGACY,
      CANCER_TYPE_DETAILED_LEGACY == "Breast Invasive Ductal Carcinoma" & CANCER_TYPE_DETAILED_PANCANCERATLAS == "Breast Invasive Carcinoma (NOS)" ~ CANCER_TYPE_DETAILED_LEGACY,
      is.na(CANCER_TYPE_DETAILED_LEGACY) ~ CANCER_TYPE_DETAILED_PANCANCERATLAS,
      CANCER_TYPE_DETAILED_LEGACY == "" & CANCER_TYPE_DETAILED_PANCANCERATLAS == "" ~ CANCER_TYPE_DETAILED_LEGACY,
      CANCER_TYPE_DETAILED_LEGACY == "Invasive Breast Carcinoma" & CANCER_TYPE_DETAILED_PANCANCERATLAS == "Breast Invasive Ductal Carcinoma" ~ CANCER_TYPE_DETAILED_PANCANCERATLAS,
      CANCER_TYPE_DETAILED_LEGACY == "Basal Cell Carcinoma" & CANCER_TYPE_DETAILED_PANCANCERATLAS == "Breast Invasive Ductal Carcinoma" ~ CANCER_TYPE_DETAILED_PANCANCERATLAS,
      TRUE ~ CANCER_TYPE_DETAILED_LEGACY
    )) %>%
  # add histotype revision according to Thennavan et al. 2021
  left_join(histology_curated %>%
              mutate(SAMPLE_ID = CLID %>% str_remove("[A-Z]$")) %>%
              select(SAMPLE_ID, Histotype_Thennavan_2021)) %>%
  # add curated menopausal state
  left_join(clinical_curated) %>%
  # for patients missing in the annotation
  mutate(menopause_status = replace_na(menopause_status, "Indeterminate/Unknown")) %>%
  select("PATIENT_ID", "SAMPLE_ID", "Age" = "AGE", "Age_10_years_increase", "Sex" = "SEX", "menopause_status", "Neoadjuvant_therapy" = "HISTORY_NEOADJUVANT_TRTYN",
         "Pam50" = "SUBTYPE",
         "Histotype", "Histotype_PanCancerAtlas" = "CANCER_TYPE_DETAILED_PANCANCERATLAS", "Histotype_Legacy" = "CANCER_TYPE_DETAILED_LEGACY", "Histotype_Thennavan_2021",
         "Aneuploidy_score" = "ANEUPLOIDY_SCORE", "TMB" = "TMB_NONSYNONYMOUS", "TBL" = "TBL_SCORE", "MSI_MANTIS" = "MSI_SCORE_MANTIS", "MSI_Sensor" = "MSI_SENSOR_SCORE",
         "OS", "OS_years", "DSS", "DSS_years", "DFS", "DFS_years", "PFS", "PFS_years"
         )

compare_histology_annotations <- xtabs(data = metadata, ~ Histotype_PanCancerAtlas + Histotype_Legacy + Histotype + Histotype_Thennavan_2021) %>%
  as.data.frame() %>%
  filter(Freq > 0)

metadata.1 <- metadata %>%
  # Remove non epithelial tumors
  filter(!Histotype %in% c("Neuroendocrine Carcinoma", "Malignant Phyllodes Tumor of the Breast")) 

# Impute menopausal state ------------------------------
ggplot(metadata.1 %>% filter(Sex == "Female")) +
  geom_quasirandom( aes(x=Age, y = menopause_status, color = menopause_status),size = 0.2) +
  geom_vline(xintercept = 45, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 55, linetype = "dashed", color = "black") +
  theme_bw() + theme(legend.position = "none")

metadata.2 <- metadata.1 %>%
  # impute menopausal state
  mutate(menopause_status_imputed = case_when(
    Sex == "Male" ~ "Indeterminate",
    menopause_status != "Indeterminate/Unknown" ~ menopause_status,
    Age < 45 ~ "Pre",
    Age >= 45 & Age < 55 ~ "Peri",
    Age >= 55 ~ "Post",
  ) %>% factor(levels = c("Post", "Peri", "Pre", "Indeterminate/Unknown"))) %>%
  relocate(menopause_status_imputed, .after = menopause_status)

# Recompute Pam50 and SCMOD2 ---------------

# recompute SCMOD2 and Pam50 classification
annot <- data.frame(
  probe        = rownames(exp_entrez),    # required but can be any unique id
  #Gene.Symbol  = rownames(exp),    # <-- key column for mapping
  EntrezGene.ID  = rownames(exp_entrez),    # <-- key column for mapping
  stringsAsFactors = FALSE,
  row.names = rownames(exp_entrez)
)

data(scmod2.robust)
data(pam50.robust)

scmod2_classification <- molecular.subtyping(
  sbt.model = "scmod2",          # the pretrained model
  data      = t(log2(exp_entrez + 1)),     # your expression matrix
  annot     = annot,         
  do.mapping = T
)

scmod2_class <- scmod2_classification[["subtype"]] %>% 
  as.data.frame() %>%
  rownames_to_column("SAMPLE_ID") %>%
  rename("SCMOD2" =  ".") %>%
  mutate(SCMOD2 = factor(SCMOD2, levels = c("ER+/HER2- Low Prolif", "ER+/HER2- High Prolif", "HER2+", "ER-/HER2-")))

pam50_classification <- molecular.subtyping(
  sbt.model = "pam50",          # the pretrained model
  data      = t(log2(exp_entrez + 1)),     # your expression matrix
  annot     = annot,         
  do.mapping = T
)

pam50_class <- pam50_classification[["subtype"]] %>% 
  as.data.frame() %>%
  rownames_to_column("SAMPLE_ID") %>%
  rename("Pam50" =  ".") %>%
  mutate(Pam50 = factor(Pam50, levels = c("LumA", "LumB", "Her2", "Basal", "Normal")))


# update metadata with new classifications --------------------
metadata.3 <- metadata.2 %>%
  left_join(pam50_class, by = "SAMPLE_ID", suffix = c("", "_NEW")) %>%
  left_join(scmod2_class, by = "SAMPLE_ID")

# some differences in Her2 and LumB vs A classification
table(metadata.3$Pam50, metadata.3$Pam50_NEW)


metadata.4 <- metadata.3 %>% 
  select(-Pam50, Pam50 = Pam50_NEW) %>%
  # reorder columns
  relocate(Pam50, SCMOD2, .before = Aneuploidy_score)

# Export --------------
write_tsv(metadata.4, "data/processed/TCGA_BRCA-metadata.tsv")
saveRDS(metadata.4, "data/processed/TCGA_BRCA-metadata.rds")

# Export counts
exp_df <- exp[,metadata.4$SAMPLE_ID] %>%
  as.data.frame() %>%
  rownames_to_column("external_gene_name") 

exp_df_zscore <- exp[,metadata.4$SAMPLE_ID] %>%
  {log2(. + 1)} %>%
  t() %>% scale() %>% t() %>%
  as.data.frame() %>%
  rownames_to_column("external_gene_name") 

fwrite(exp_df, "data/processed/deduplicated_and_filtered-data_mrna_seq_v2_rsem.tsv.gz", compress = "gzip", sep = "\t", quote = FALSE)
fwrite(exp_df_zscore, "data/processed/deduplicated_and_filtered-data_mrna_seq_v2_rsem_log_zscores.tsv.gz", compress = "gzip", sep = "\t", quote = FALSE)
