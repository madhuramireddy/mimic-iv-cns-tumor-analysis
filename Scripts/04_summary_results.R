# 04_summary_results.R
# Combine results from antidepressant and statin analyses into one summary table.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ---- paths from setup -------------------------------------------------------------
source(file.path("Scripts", "00_setup.R"))

path_antidep_analysis <- file.path("Output", "antidepressant_analysis.tsv")
path_statin_analysis  <- file.path("Output", "statin_analysis.tsv")

# ---- helper to compute stats ------------------------------------------------------
compute_stats <- function(df, exposure_label) {
  tab <- table(df$cns_tumor, df$exposed)
  ft  <- fisher.test(tab)
  
  tibble::tibble(
    exposure   = exposure_label,
    n_total    = nrow(df),
    n_exposed  = sum(df$exposed),
    n_cns      = sum(df$cns_tumor),
    n_overlap  = sum(df$exposed & df$cns_tumor),
    odds_ratio = as.numeric(ft$estimate),
    ci_lower   = ft$conf.int[1],
    ci_upper   = ft$conf.int[2],
    p_value    = ft$p.value
  ) %>%
    mutate(
      odds_ratio = round(odds_ratio, 3),
      ci_lower   = round(ci_lower,   3),
      ci_upper   = round(ci_upper,   3),
      p_value    = signif(p_value,   3)
    )
}

# ---- read analysis outputs --------------------------------------------------------
df_antidep <- fread(path_antidep_analysis)
df_statin  <- fread(path_statin_analysis)

# ---- compute results --------------------------------------------------------------
res_antidep <- compute_stats(df_antidep, "Antidepressants")
res_statin  <- compute_stats(df_statin,  "Statins")

summary_tbl <- bind_rows(res_antidep, res_statin)

# ---- write summary ----------------------------------------------------------------
out_summary <- file.path("Output", "exposure_summary.tsv")
fwrite(summary_tbl, out_summary, sep = "\t")

# ---- console summary --------------------------------------------------------------
cat("\n[exposure_summary]\n")
print(summary_tbl)
cat("\nResults written to: ", out_summary, "\n")
