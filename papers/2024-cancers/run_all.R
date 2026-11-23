# Run all analyses for paper 2024-cancers (file papers/2024-cancers/run_all.R)
# Master script to reproduce all figures for Cumming et al. 2024 (Cancers).
# Run from the project root (where MODEL.Rproj lives, so renv activates):
#
#   Rscript papers/2024-cancers/run_all.R
#
# Prerequisites:
#   1. R >= 4.2
#   2. renv::restore()               # restore the package library
#   3. git lfs pull                  # download the derived .rds files
#      (or regenerate them: see "Regenerating the .rds files" in
#      papers/2024-cancers/README.md)
#
# By default, 01_bootstrap_params.R is SKIPPED because its output
# (data-derived/resampled_pars_fin.rds) is shipped through Git LFS.
# To regenerate it:
#   Rscript papers/2024-cancers/run_all.R --bootstrap

library(here)

paper_dir <- here("papers", "2024-cancers")
args <- commandArgs(trailingOnly = TRUE)
run_bootstrap <- "--bootstrap" %in% args

if (run_bootstrap) {
  message("==> Running 01_bootstrap_params.R ...")
  source(file.path(paper_dir, "analysis", "01_bootstrap_params.R"))
} else {
  message(
    "==> Skipping 01_bootstrap_params.R ",
    "(using shipped data-derived/ artifacts)."
  )
  message("    Pass --bootstrap to regenerate from scratch.")
}

message("==> Running 02_sensitivity.R ...")
source(file.path(paper_dir, "analysis", "02_sensitivity.R"))

message("==> Running 03_progression_model.R ...")
source(file.path(paper_dir, "analysis", "03_progression_model.R"))

message("==> Running 04_clinical_scenarios.R ...")
source(file.path(paper_dir, "analysis", "04_clinical_scenarios.R"))

message("==> Done. Figures written to ", file.path(paper_dir, "figures"))
