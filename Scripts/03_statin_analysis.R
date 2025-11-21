# 03_statin_analysis.R
# Build statin exposure flags across the full universe,
# merge with CNS tumor cohort, and compute association statistics.

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

# Statin list
st_raw <- fread(path_statins, header = TRUE)

# Collect name-like columns and clean them
name_cols <- intersect(c("synonym","inn","name","V3","V5","V6"), names(st_raw))
st_names <- st_raw %>%
  select(all_of(name_cols)) %>%
  as.list() %>%
  unlist(use.names = FALSE)

st_names <- tolower(trimws(st_names))
st_names <- unique(st_names[!is.na(st_names) & nzchar(st_names)])
st_names <- st_names[!grepl("^level_|^code_|^c10", st_names)]  # drop metadata tokens

# ---- exposure flagging ------------------------------------------------------------
pattern <- paste0("(", paste0(str_replace_all(st_names, "([\\W])", "\\\\\\1"),
                              collapse = "|"), ")")

statin_ids <- prescriptions %>%
  mutate(drug_lc = tolower(drug)) %>%
  filter(!is.na(drug_lc) & str_detect(drug_lc, pattern)) %>%
  distinct(subject_id)

# ---- merge into universe ----------------------------------------------------------
df <- all_subjects %>%
  mutate(
    exposed   = subject_id %in% statin_ids$subject_id,
    cns_tumor = subject_id %in% cns_cohort$subject_id
  )

# Write enriched table
out_statin_analysis <- file.path("Output", "statin_analysis.tsv")
fwrite(df, out_statin_analysis, sep = "\t")

# ---- 2x2 table + Fisher test ------------------------------------------------------
tab <- table(df$cns_tumor, df$exposed)
ft  <- fisher.test(tab)

cat("\n[statin_analysis]\n")
print(tab)
cat("\nOdds ratio: ", round(unname(ft$estimate), 3),
    "\n95% CI: ", round(ft$conf.int[1], 3), " - ", round(ft$conf.int[2], 3),
    "\np-value: ", signif(ft$p.value, 3), "\n", sep = "")
