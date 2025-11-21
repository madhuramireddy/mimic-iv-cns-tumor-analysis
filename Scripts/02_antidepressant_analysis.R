# 02_antidepressant_analysis.R
# Build antidepressant exposure flags across the full universe, merge with CNS tumor cohort, and compute association statistics.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ---- paths from setup -------------------------------------------------------------
source(file.path("Scripts", "00_setup.R"))

# ---- inputs -----------------------------------------------------------------------
prescriptions <- fread(out_prescriptions, showProgress = FALSE)
all_subjects  <- fread(out_universe,      showProgress = FALSE)
cns_cohort    <- fread(out_cns_cohort,    showProgress = FALSE)

# Antidepressant list
ad_raw <- fread(path_antidep, header = TRUE)
ad_names <- if ("Generic" %in% names(ad_raw)) ad_raw$Generic else ad_raw[[2]]
ad_names <- tolower(trimws(ad_names))
ad_names <- unique(ad_names[!is.na(ad_names) & nzchar(ad_names) & ad_names != "generic"])

# ---- exposure flagging ------------------------------------------------------------
pattern <- paste0("(", paste0(str_replace_all(ad_names, "([\\W])", "\\\\\\1"),
                              collapse = "|"), ")")

antidep_ids <- prescriptions %>%
  mutate(drug_lc = tolower(drug)) %>%
  filter(!is.na(drug_lc) & str_detect(drug_lc, pattern)) %>%
  distinct(subject_id)

# ---- merge into universe ----------------------------------------------------------
df <- all_subjects %>%
  mutate(
    exposed   = subject_id %in% antidep_ids$subject_id,
    cns_tumor = subject_id %in% cns_cohort$subject_id
  )

# Write enriched table
out_antidep_analysis <- file.path("Output", "antidepressant_analysis.tsv")
fwrite(df, out_antidep_analysis, sep = "\t")

# ---- 2x2 table + Fisher test ------------------------------------------------------
tab <- table(df$cns_tumor, df$exposed)
ft  <- fisher.test(tab)

cat("\n[antidepressant_analysis]\n")
print(tab)
cat("\nOdds ratio: ", round(unname(ft$estimate), 3),
    "\n95% CI: ", round(ft$conf.int[1], 3), " - ", round(ft$conf.int[2], 3),
    "\np-value: ", signif(ft$p.value, 3), "\n", sep = "")
