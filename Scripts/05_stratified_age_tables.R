# 05_stratified_age_tables.R
# Age-stratified analysis:
# Prints full 2x2 tables + Fisher's test (no tidy summaries).

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ---- load paths from setup -------------------------------------------------------
source(file.path("Scripts", "00_setup.R"))

# ---- load data -------------------------------------------------------------------
antidep_df <- fread(file.path("Output", "antidepressant_analysis.tsv"))
statin_df  <- fread(file.path("Output", "statin_analysis.tsv"))
patients   <- fread(out_patients)    # subject_id + anchor_age

# ---- define age groups -----------------------------------------------------------
patients_age <- patients %>%
  mutate(anchor_age = suppressWarnings(as.numeric(anchor_age))) %>%
  filter(!is.na(anchor_age)) %>%
  mutate(
    age_group = case_when(
      anchor_age < 19 ~ "<19",
      anchor_age >= 19 & anchor_age <= 40 ~ "19-40",
      anchor_age > 40 ~ ">40",
      TRUE ~ NA_character_
    )
  ) %>%
  select(subject_id, age_group)

# ---- merge age into exposure tables ---------------------------------------------
antidep_age <- antidep_df %>% left_join(patients_age, by = "subject_id")
statin_age  <- statin_df  %>% left_join(patients_age, by = "subject_id")

# ---- table function ----------------------------------------------------
run_tables <- function(df, age_label, exposure_label) {
  
  df_sub <- df %>% filter(age_group == age_label)
  
  # Build 2x2 table
  tab <- table(df_sub$cns_tumor, df_sub$exposed)
  
  cat("\n", exposure_label, "byAge", age_label, " — Table 1 (Cancer x Exposure):\n", sep = "")
  print(tab)
  
  # Fisher test
  ft <- fisher.test(tab)
  
  cat("\nOdds ratio: ", round(unname(ft$estimate), 3), "\n",
      "95% CI: ", round(ft$conf.int[1], 3), " - ", round(ft$conf.int[2], 3), "\n",
      "p-value: ", signif(ft$p.value, 3), "\n\n", sep = "")
}

# ---- run all age levels ----------------------------------------------------------
age_levels <- c("<19", "19-40", ">40")

for (ag in age_levels) {
  run_tables(antidep_age, ag, "Antidepressants")
  run_tables(statin_age,  ag, "Statins")
}
