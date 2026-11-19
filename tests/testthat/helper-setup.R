# Shared test setup (testthat sources helper-*.R before the tests).
# Loads the model code exactly as the analysis scripts do, plus the packages
# that neg_log_lik() needs (hypogeom for dhypogeom, DescTools for dOrder).

suppressPackageStartupMessages({
  library(testthat)
  library(hypogeom)
  library(DescTools)
})

source(here::here("R", "setup.R"))
source(here::here("R", "model_core.R"))
source(here::here("R", "bootstrap.R"))

# Observed data and the hard-coded point estimates used by analysis/01
# (copied verbatim from analysis/01_bootstrap_params.R; neg_log_lik()'s
# default `meds` argument looks for these globals).
med_0 <- 38
med_a <- 25
med_b <- 10
med_c <- 7
med_d <- 4
meds <- c(med_0, med_a, med_b, med_c, med_d)
hypo_probs <- c(
  0.07982013,
  0.05540597,
  0.33922926,
  0.35681215,
  0.16023503,
  0.21878495,
  0.54017653,
  0.01015650
)

resampled_pars_fin <- readRDS(
  here::here("data-derived", "resampled_pars_fin.rds")
)

csv_main <- here::here("data", "Parameters.csv")
csv_25 <- here::here("data", "Parameters_25.csv")
read_params <- function(csv) read.csv(csv, row.names = 1)

golden_path <- function()
  here::here("tests", "testthat", "fixtures", "golden.rds")

# Evaluate the statements of a function body in order, inside `env`, until the
# assignment `<name> <- ...` has run, and return that value. Lets tests reach
# matrices that the code builds inline and does not expose.
extract_assigned <- function(fn, name, env) {
  for (stmt in as.list(body(fn))[-1]) {
    eval(stmt, env)
    if (
      is.call(stmt) &&
        identical(stmt[[1]], as.name("<-")) &&
        identical(as.character(stmt[[2]]), name)
    ) {
      return(get(name, envir = env))
    }
  }
  stop("assignment to '", name, "' not found in function body")
}

# The 15x15 matrix built inside neg_log_lik()
nll_matrix <- function(probs = hypo_probs) {
  env <- list2env(list(hypo_probs = probs), parent = globalenv())
  extract_assigned(neg_log_lik, "matrix_of_trans", env)
}
