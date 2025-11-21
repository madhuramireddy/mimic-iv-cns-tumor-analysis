# 00_setup.R
# Project setup: create main subdirectories and define key paths for inputs/outputs.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(stringr)
  library(readxl)
})

# ---- Project directories ----------------------------------------------------------
project_dir <- getwd()

# Core subdirectories under the project root
subdirs <- c("Data", "Output", "Scripts")
for (subdir in subdirs) {
  dir_path <- file.path(project_dir, subdir)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Created directory: ", dir_path)
  } else {
    message("Directory already exists: ", dir_path)
  }
}

# ---- Input file paths -------------------------------------------------------------
# MIMIC data
dir_mimic <- file.path(project_dir, "Data", "mimic-iv-3-2.1", "hosp")
path_diagnoses     <- file.path(dir_mimic, "diagnoses_icd.csv.gz")
path_prescriptions <- file.path(dir_mimic, "prescriptions.csv.gz")
path_patients      <- file.path(dir_mimic, "patients.csv.gz")

# Reference lists (your renamed files)
path_icd_defs <- file.path(project_dir, "Data", "Depression_Cancer_ICD_Codes.xlsx")
path_antidep  <- file.path(project_dir, "Data", "Antidepressants.csv")
path_statins  <- file.path(project_dir, "Data", "Statins.csv")

# ---- Output naming (plain .tsv, no compression) -----------------------------------
out_diagnoses     <- file.path(project_dir, "Output", "diagnoses_filtered.tsv")
out_prescriptions <- file.path(project_dir, "Output", "prescriptions_filtered.tsv")
out_patients      <- file.path(project_dir, "Output", "patients_filtered.tsv")

out_cns_cohort <- file.path(project_dir, "Output", "cns_tumor_cohort.tsv")
out_universe   <- file.path(project_dir, "Output", "all_subjects_universe.tsv")

message("Subdirectory setup complete. Paths initialized.")
