# Shared test setup (testthat sources helper-*.R before the tests).
# Loads the model code exactly as the analysis scripts do, plus the packages
# that neg_log_lik() needs (hypogeom for dhypogeom, DescTools for dOrder).

suppressPackageStartupMessages({
  library(testthat)
  library(hypogeom)
  library(DescTools)
})

options(hcc.paper = "2024-cancers") # paths: paper_path() in R/setup.R
source(here::here("R", "setup.R"))
source(here::here("R", "model_core.R"))
source(here::here("R", "bootstrap.R"))
source(here::here("R", "scenarios.R"))
source(here::here("R", "plot_theme.R"))
source(here::here("R", "plot_sweeps.R"))
source(here::here("R", "plot_scenarios.R"))
source(here::here("R", "plot_calibration.R"))

# The 8 likelihood parameters at the CSV point estimates (= paper Tables 1-2 =
# the Nelder-Mead MLE; finding G7), in the order of NLL_PAR_NAMES.
csv_par <- unlist(params_def[NLL_PAR_NAMES])

resampled_pars_fin <- readRDS(
  paper_path("data-derived", "resampled_pars_fin.rds")
)

csv_main <- paper_path("data", "Parameters.csv")
csv_25 <- paper_path("data", "Parameters_25.csv")
read_params <- function(csv) read.csv(csv, row.names = 1)

golden_path <- function() {
  here::here("tests", "testthat", "fixtures", "golden.rds")
}

# The 15x15 matrix used inside neg_log_lik(), for the 8-vector `probs`
nll_matrix <- function(probs = csv_par) {
  likelihood_matrix(probs[1], probs[2], probs[6], probs[7])
}
