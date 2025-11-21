# 06_stratified_sex.R
# Sex-stratified analysis:
# Prints full 2x2 tables (Cancer x Exposure) + Fisher's test for antidepressants and statins, by sex.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ---- load paths from setup -------------------------------------------------------
source(file.path("Scripts", "00_setup.R"))

# ---- load data -------------------------------------------------------------------
# Enriched exposure tables from earlier scripts
antidep_df <- fread(file.path("Output", "antidepressant_analysis.tsv"))
statin_df  <- fread(file.path("Output", "statin_analysis.tsv"))

# Full patients table from MIMIC to get gender
patients_full <- fread(
  path_patients,
  select = c("subject_id", "gender"),
  showProgress = FALSE
)

# ---- build sex variable ----------------------------------------------------------
sex_map <- patients_full %>%
  mutate(
    sex = dplyr::case_when(
      gender %in% c("M", "m") ~ "Male",
      gender %in% c("F", "f") ~ "Female",
      TRUE ~ NA_character_
    )
  ) %>%
  select(subject_id, sex)

# ---- merge sex into exposure tables ---------------------------------------------
antidep_sex <- antidep_df %>%
  left_join(sex_map, by = "subject_id")

statin_sex <- statin_df %>%
  left_join(sex_map, by = "subject_id")

# ---- helper to print 2x2 tables + Fisher results --------------------------------
run_tables <- function(df, sex_label, exposure_label) {
  df_sub <- df %>%
    filter(sex == sex_label)
  
  if (nrow(df_sub) == 0) {
    cat("\n", exposure_label, "bySex", sex_label,
        ": no subjects in this group.\n", sep = "")
    return(invisible(NULL))
  }
  
  tab <- table(df_sub$cns_tumor, df_sub$exposed)
  
  cat("\n", exposure_label, "bySex", sex_label,
      " — Table 1 (Cancer x Exposure):\n", sep = "")
  print(tab)
  
  ft <- fisher.test(tab)
  
  cat("\nOdds ratio: ", round(unname(ft$estimate), 3), "\n",
      "95% CI: ", round(ft$conf.int[1], 3), " - ", round(ft$conf.int[2], 3), "\n",
      "p-value: ", signif(ft$p.value, 3), "\n\n", sep = "")
}

# ---- run for each sex level ------------------------------------------------------
sex_levels <- c("Male", "Female")
# Alternatively, derive from data:
# sex_levels <- sort(unique(na.omit(sex_map$sex)))

for (sx in sex_levels) {
  run_tables(antidep_sex, sx, "Antidepressants")
  run_tables(statin_sex,  sx, "Statins")
}
