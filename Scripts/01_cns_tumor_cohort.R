# 01_cns_tumor_cohort.R
# Load core MIMIC-IV tables, apply adult filter, import ICD definitions, and construct CNS tumor cohort.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(stringr)
  library(readxl)
})

# paths and output filenames
source(file.path("Scripts", "00_setup.R"))

# ---- switches --------------------------------------------------------------------
ADULTS_ONLY <- TRUE

# ---- helpers ---------------------------------------------------------------------
normalize_icd <- function(x) {
  x <- as.character(x)
  x <- toupper(x)
  x <- gsub("\\.", "", x)
  x <- gsub("\\s+", "", x)
  x
}

clean_icd_vec <- function(v) {
  v <- as.character(v)
  v <- v[!is.na(v) & nzchar(v)]
  v <- toupper(gsub("\\.|\\s+", "", v))
  unique(v)
}

# ---- load MIMIC (minimal columns) -------------------------------------------------
diagnoses <- fread(
  path_diagnoses,
  select = c("subject_id","hadm_id","icd_code","icd_version"),
  showProgress = FALSE
)
prescriptions <- fread(
  path_prescriptions,
  select = c("subject_id","hadm_id","drug","ndc","gsn"),
  showProgress = FALSE
)
patients <- fread(
  path_patients,
  select = c("subject_id","anchor_age"),
  showProgress = FALSE
)

# ICD to joinable form
diagnoses <- diagnoses %>% mutate(icd_code = normalize_icd(icd_code))

# ---- adult restriction -------------------------------------------------------------
if (ADULTS_ONLY) {
  adult_ids <- patients %>%
    mutate(anchor_age = suppressWarnings(as.numeric(anchor_age))) %>%
    filter(!is.na(anchor_age), anchor_age >= 18) %>%
    distinct(subject_id)
  
  diagnoses     <- diagnoses     %>% semi_join(adult_ids, by = "subject_id")
  prescriptions <- prescriptions %>% semi_join(adult_ids, by = "subject_id")
}

# ---- ICD definitions (Excel, multiple sheets) -------------------------------------
# Only CNS tumor sheets are used now
icd9_cns <- read_excel(path_icd_defs, sheet = "ICD9_CNS") %>%
  transmute(icd_code = normalize_icd(CODE)) %>%
  filter(!is.na(icd_code) & nzchar(icd_code)) %>%
  distinct()

icd10_cns <- read_excel(path_icd_defs, sheet = "ICD10_CNS") %>%
  transmute(icd_code = normalize_icd(CODE)) %>%
  filter(!is.na(icd_code) & nzchar(icd_code)) %>%
  distinct()

# ---- cohort: CNS tumors ------------------------------------------------------------
diag9  <- diagnoses %>% filter(icd_version %in% c(9L, "9"))
diag10 <- diagnoses %>% filter(icd_version %in% c(10L, "10"))

cns_icd9  <- diag9  %>% semi_join(icd9_cns,  by = "icd_code") %>% distinct(subject_id)
cns_icd10 <- diag10 %>% semi_join(icd10_cns, by = "icd_code") %>% distinct(subject_id)
cns_ids   <- bind_rows(cns_icd9, cns_icd10) %>% distinct(subject_id)

# ---- universe ---------------------------------------------------------------------
all_subjects_universe <- bind_rows(
  diagnoses     %>% distinct(subject_id),
  prescriptions %>% distinct(subject_id)
) %>% distinct(subject_id)

# ---- writes (plain .tsv) ----------------------------------------------------------
fwrite(diagnoses,     out_diagnoses,     sep = "\t")
fwrite(prescriptions, out_prescriptions, sep = "\t")
fwrite(patients,      out_patients,      sep = "\t")

fwrite(cns_ids,       out_cns_cohort,    sep = "\t")
fwrite(all_subjects_universe, out_universe, sep = "\t")

# ---- console summary --------------------------------------------------------------
cat(
  "\n[extract_patients_icd]\n",
  "diagnoses_filtered: ", nrow(diagnoses), "\n",
  "prescriptions_filtered: ", nrow(prescriptions), "\n",
  "patients_filtered: ", nrow(patients), "\n",
  "cns_tumor_cohort: ", nrow(cns_ids), "\n",
  "all_subjects_universe: ", nrow(all_subjects_universe), "\n\n",
  sep = ""
)
