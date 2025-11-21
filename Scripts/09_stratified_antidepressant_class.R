# 09_stratified_antidepressant_class.R
# Stratified analysis by antidepressant class (SSRI, SNRI, TCA).
# Prints 2x2 tables (Cancer x Exposure) + Fisher's test for each class.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ---- load paths ------------------------------------------------------------------
source(file.path("Scripts", "00_setup.R"))

# helper: normalize text for matching
clean_lower <- function(x) {
  tolower(trimws(as.character(x)))
}

# regex OR helper (same idea as in your earlier script)
rx_or <- function(v) {
  v <- unique(v[!is.na(v) & nzchar(v)])
  if (length(v) == 0) return("a^")  # match nothing if empty
  pattern_parts <- stringr::str_replace_all(v, "([\\W])", "\\\\\\1")
  paste0("(", paste0(pattern_parts, collapse = "|"), ")")
}

# ---- load data -------------------------------------------------------------------
all_subjects <- fread(out_universe,      showProgress = FALSE)  # subject_id
cns_cohort   <- fread(out_cns_cohort,    showProgress = FALSE)  # subject_id
prescriptions <- fread(out_prescriptions,
                       select = c("subject_id","drug"),
                       showProgress = FALSE)

prescriptions <- prescriptions %>%
  mutate(drug_lc = clean_lower(drug))

# Antidepressant list with class info
ad_raw <- fread(path_antidep, header = TRUE)

# Expect columns: "Brand Name", "Generic", "Class"
# If names differ slightly, adjust here.
if (!all(c("Generic","Class") %in% names(ad_raw))) {
  stop("Antidepressants file must contain at least 'Generic' and 'Class' columns.")
}

ad_raw <- ad_raw %>%
  mutate(
    generic_lc = clean_lower(Generic),
    class      = trimws(as.character(Class))
  ) %>%
  filter(!is.na(generic_lc), nzchar(generic_lc),
         !is.na(class),      nzchar(class))

# ---- function to get exposed subjects for a given AD class -----------------------

get_exposed_ids_for_class <- function(class_name) {
  generics <- ad_raw %>%
    filter(class == class_name) %>%
    pull(generic_lc) %>%
    unique()
  
  if (length(generics) == 0) {
    warning("No generics found for class: ", class_name)
    return(data.frame(subject_id = integer(0)))
  }
  
  pattern <- rx_or(generics)
  
  exposed <- prescriptions %>%
    filter(!is.na(drug_lc) & str_detect(drug_lc, pattern)) %>%
    distinct(subject_id)
  
  exposed
}

# ---- helper to print 2x2 + Fisher -----------------------------------------------

run_tables <- function(exposed_ids, cancer_ids, label) {
  df <- all_subjects %>%
    mutate(
      exposed = subject_id %in% exposed_ids$subject_id,
      cancer  = subject_id %in% cancer_ids$subject_id
    )
  
  tab <- table(df$cancer, df$exposed)  # rows: cancer (FALSE/TRUE), cols: exposed (FALSE/TRUE)
  
  cat("\n", label, " — Table 1 (Cancer x Exposure):\n", sep = "")
  print(tab)
  
  ft <- fisher.test(tab)
  
  cat("\nOdds ratio: ", round(unname(ft$estimate), 3), "\n",
      "95% CI: ", round(ft$conf.int[1], 3), " - ", round(ft$conf.int[2], 3), "\n",
      "p-value: ", signif(ft$p.value, 3), "\n\n", sep = "")
}

# ---- run for each antidepressant class ------------------------------------------

classes <- c("SSRI", "SNRI", "TCA")

for (cl in classes) {
  exposed_ids <- get_exposed_ids_for_class(cl)
  label <- paste0("Antidepressants-Class-", cl)
  run_tables(exposed_ids, cns_cohort, label)
}
