# 07_stratified_neoplasm_type.R
# Stratified analysis by neoplasm type (benign vs malignant).
# Uses explicit ICD-9/10 code lists.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ---- load paths ------------------------------------------------------------------
source(file.path("Scripts", "00_setup.R"))

normalize_icd <- function(x) {
  x <- as.character(x)
  x <- toupper(x)
  x <- gsub("\\.", "", x)
  x <- gsub("\\s+", "", x)
  x
}

# ---- load data -------------------------------------------------------------------
all_subjects <- fread(out_universe, showProgress = FALSE)
diagnoses    <- fread(out_diagnoses,
                      select = c("subject_id","icd_code","icd_version"),
                      showProgress = FALSE)

diagnoses <- diagnoses %>%
  mutate(icd_code = normalize_icd(icd_code))

# Exposure tables
antidep_df <- fread(file.path("Output", "antidepressant_analysis.tsv"))
statin_df  <- fread(file.path("Output", "statin_analysis.tsv"))

antidep_ids <- antidep_df %>% filter(exposed) %>% distinct(subject_id)
statin_ids  <- statin_df  %>% filter(exposed) %>% distinct(subject_id)

# ---- explicit ICD lists for benign / malignant -----------------------------------

neoplasm_icd9_benign_list <- c(
  "2250","2251","2252","2253","2254","2258","2259","2273","2274"
)

neoplasm_icd10_benign_list <- c(
  "D320","D321","D329","D330","D331","D332","D333",
  "D334","D337","D339","D352","D354"
)

neoplasm_icd9_malignant_list <- c(
  "1910","1911","1912","1913","1914","1915","1916",
  "1917","1918","1919","1920","1921","1922","1923",
  "1928","1929","1943","1944","1983","1984"
)

neoplasm_icd10_malignant_list <- c(
  "C700","C701","C709","C710","C711","C712","C713",
  "C714","C715","C716","C717","C718","C719","C720",
  "C721","C7220","C7221","C7222","C7230","C7231",
  "C7232","C7240","C7241","C7242","C7250","C7259",
  "C729","C753","C751","C7931","C7932","C7940","C7949"
)

# ---- identify benign / malignant patients ----------------------------------------

diag9  <- diagnoses %>% filter(icd_version %in% c(9L,"9"))
diag10 <- diagnoses %>% filter(icd_version %in% c(10L,"10"))

# benign
cancer_benign <- bind_rows(
  diag9  %>% filter(icd_code %in% neoplasm_icd9_benign_list),
  diag10 %>% filter(icd_code %in% neoplasm_icd10_benign_list)
) %>% distinct(subject_id)

# malignant
cancer_malignant <- bind_rows(
  diag9  %>% filter(icd_code %in% neoplasm_icd9_malignant_list),
  diag10 %>% filter(icd_code %in% neoplasm_icd10_malignant_list)
) %>% distinct(subject_id)

# ---- function to print 2x2 + fisher ----------------------------------------------
run_tables <- function(exposed_ids, cancer_ids, label) {
  
  df <- all_subjects %>%
    mutate(
      exposed = subject_id %in% exposed_ids$subject_id,
      cancer  = subject_id %in% cancer_ids$subject_id
    )
  
  tab <- table(df$cancer, df$exposed)
  
  cat("\n", label, " — Table 1 (Cancer x Exposure):\n", sep = "")
  print(tab)
  
  ft <- fisher.test(tab)
  
  cat("\nOdds ratio: ", round(unname(ft$estimate), 3), "\n",
      "95% CI: ", round(ft$conf.int[1], 3), " - ", round(ft$conf.int[2], 3), "\n",
      "p-value: ", signif(ft$p.value, 3), "\n\n", sep = "")
}

# ---- run analyses ----------------------------------------------------------------

# Benign
run_tables(antidep_ids, cancer_benign,
           "Antidepressants-NeoplasmType-Benign")

run_tables(statin_ids,  cancer_benign,
           "Statins-NeoplasmType-Benign")

# Malignant
run_tables(antidep_ids, cancer_malignant,
           "Antidepressants-NeoplasmType-Malignant")

run_tables(statin_ids,  cancer_malignant,
           "Statins-NeoplasmType-Malignant")
