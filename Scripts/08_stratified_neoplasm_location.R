# 08_stratified_neoplasm_location.R
# Stratified analysis by neoplasm location for malignant CNS tumors:
# Cerebrum vs Cerebellum vs Brainstem/Spinal cord.
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
all_subjects <- fread(out_universe, showProgress = FALSE)         # subject_id
diagnoses    <- fread(out_diagnoses,
                      select = c("subject_id","icd_code","icd_version"),
                      showProgress = FALSE)

diagnoses <- diagnoses %>%
  mutate(icd_code = normalize_icd(icd_code))

# Exposure tables (already have exposed flag)
antidep_df <- fread(file.path("Output", "antidepressant_analysis.tsv"))
statin_df  <- fread(file.path("Output", "statin_analysis.tsv"))

antidep_ids <- antidep_df %>% filter(exposed) %>% distinct(subject_id)
statin_ids  <- statin_df  %>% filter(exposed) %>% distinct(subject_id)

# ---- explicit ICD lists for location (malignant CNS) -----------------------------

neoplasm_icd9_cerebrum_list <- c("1910","1911","1912","1913","1914")
neoplasm_icd10_cerebrum_list <- c("C710","C711","C712","C713","C714")

neoplasm_icd9_cerebellum_list <- c("1916")
neoplasm_icd10_cerebellum_list <- c("C716")

neoplasm_icd9_stem_list <- c("1917","1922")
neoplasm_icd10_stem_list <- c("C717","C720")

# ---- identify subjects by location ----------------------------------------------

diag9  <- diagnoses %>% filter(icd_version %in% c(9L,"9"))
diag10 <- diagnoses %>% filter(icd_version %in% c(10L,"10"))

# Cerebrum
cns_cerebrum <- bind_rows(
  diag9  %>% filter(icd_code %in% neoplasm_icd9_cerebrum_list),
  diag10 %>% filter(icd_code %in% neoplasm_icd10_cerebrum_list)
) %>% distinct(subject_id)

# Cerebellum
cns_cerebellum <- bind_rows(
  diag9  %>% filter(icd_code %in% neoplasm_icd9_cerebellum_list),
  diag10 %>% filter(icd_code %in% neoplasm_icd10_cerebellum_list)
) %>% distinct(subject_id)

# Brainstem / Spinal cord
cns_stem <- bind_rows(
  diag9  %>% filter(icd_code %in% neoplasm_icd9_stem_list),
  diag10 %>% filter(icd_code %in% neoplasm_icd10_stem_list)
) %>% distinct(subject_id)

# ---- helper to print 2x2 + Fisher test ------------------------------------------

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

# ---- run analyses for each location ---------------------------------------------

# Cerebrum
run_tables(antidep_ids, cns_cerebrum,
           "Antidepressants-NeoplasmLocation-Cerebrum")
run_tables(statin_ids,  cns_cerebrum,
           "Statins-NeoplasmLocation-Cerebrum")

# Cerebellum
run_tables(antidep_ids, cns_cerebellum,
           "Antidepressants-NeoplasmLocation-Cerebellum")
run_tables(statin_ids,  cns_cerebellum,
           "Statins-NeoplasmLocation-Cerebellum")

# Brainstem / Spinal cord
run_tables(antidep_ids, cns_stem,
           "Antidepressants-NeoplasmLocation-Brainstem_SpinalCord")
run_tables(statin_ids,  cns_stem,
           "Statins-NeoplasmLocation-Brainstem_SpinalCord")
