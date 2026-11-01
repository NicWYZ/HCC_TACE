# Propensity score matching for the three pairwise sequencing comparisons.
#
# Input:  data/processed/df_{tace_sys,tace_same,same_sys}.csv (from notebooks/01)
# Output: data/processed/matched_*.csv and love plots in results/figures/
#
# Run from the repository root:  Rscript R/psm_matching.R

suppressPackageStartupMessages({
  library(MatchIt)
  library(cobalt)
  library(readr)
  library(ggplot2)
})

processed_dir <- file.path("data", "processed")
figures_dir <- file.path("results", "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

ps_formula <- sequence_group ~ AGE + SEX + RACE + INSURANCE_STATUS +
  Charlson_Deyo_Score + HISTOLOGY + COMBINED_STAGE_GROUP +
  YEAR_OF_DIAGNOSIS + Surgery + Radiation_Therapy

# Comparison -> matching ratio (controls per treated patient), caliper 0.2 SD
comparisons <- list(
  tace_sys  = 4,  # TACE first (0) vs systemic first (1)
  tace_same = 3,  # TACE first (0) vs simultaneous (1)
  same_sys  = 1   # simultaneous (0) vs systemic first (1)
)

for (name in names(comparisons)) {
  df <- read_csv(file.path(processed_dir, paste0("df_", name, ".csv")), show_col_types = FALSE)

  m_out <- matchit(ps_formula, data = df, method = "nearest",
                   ratio = comparisons[[name]], caliper = 0.2)

  cat("\n====", name, "====\n")
  print(summary(m_out, standardize = TRUE))

  p <- love.plot(m_out, binary = "std", thresholds = c(m = .1),
                 var.order = "unadjusted", abs = TRUE)
  ggsave(file.path(figures_dir, paste0("love_", name, ".pdf")), p, width = 8, height = 6)

  write_csv(match.data(m_out), file.path(processed_dir, paste0("matched_", name, ".csv")))
}
