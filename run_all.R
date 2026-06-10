# run_all.R
# Master script to reproduce all figures from the published HCC surveillance
# model.  Run from the project root (where MODEL.Rproj lives):
#
#   Rscript run_all.R
#
# Prerequisites:
#   1. R >= 4.2
#   2. renv::restore()               # restore the package library
#   3. git lfs pull                  # download data-derived/ RDS files
#
# By default, the expensive bootstrap estimation step (01_bootstrap_params.R)
# is SKIPPED because its outputs are already shipped in data-derived/.
# To regenerate those files from scratch (requires HPC-scale resources):
#   Rscript run_all.R --bootstrap

library(here)

args <- commandArgs(trailingOnly = TRUE)
run_bootstrap <- "--bootstrap" %in% args

if (run_bootstrap) {
  message(
    "==> Running bootstrap parameter estimation (this will take a long time)..."
  )
  source(here("analysis", "01_bootstrap_params.R"))
} else {
  message(
    "==> Skipping 01_bootstrap_params.R (using shipped data-derived/ artifacts)."
  )
  message("    Pass --bootstrap to regenerate from scratch.")
}

message("==> Running 02_sensitivity.R ...")
source(here("analysis", "02_sensitivity.R"))

message("==> Running 03_progression_model.R ...")
source(here("analysis", "03_progression_model.R"))

message("==> Running 04_clinical_scenarios.R ...")
source(here("analysis", "04_clinical_scenarios.R"))

message("==> Done. Figures written to figures/")
