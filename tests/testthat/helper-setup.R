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

# The 7x7 matrix built inside fromf4() (after its annual->monthly conversion)
fromf4_tmat <- function(
  n = 100,
  t = 120,
  p_HCC = p_HCC_def,
  p_0 = p_0_def,
  p_A = p_A_def,
  p_B = p_B_def,
  p_C = p_C_def,
  p_D = p_D_def
) {
  env <- list2env(
    list(
      n = n,
      t = t,
      p_HCC = p_HCC,
      p_0 = p_0,
      p_A = p_A,
      p_B = p_B,
      p_C = p_C,
      p_D = p_D
    ),
    parent = globalenv()
  )
  extract_assigned(fromf4, "tmat", env)
}

# The 7x7 matrix built inside us_theor_pa()
us_theor_pa_tmat <- function(
  p_HCC = p_HCC_def,
  p_0 = p_0_def,
  p_A = p_A_def,
  p_B = p_B_def,
  p_C = p_C_def,
  p_D = p_D_def
) {
  env <- list2env(
    list(
      n = 100,
      t = 120,
      m = 6,
      Pr_Att = 1,
      p_HCC = p_HCC,
      p_0 = p_0,
      p_A = p_A,
      p_B = p_B,
      p_C = p_C,
      p_D = p_D,
      Sens_0 = Sens_0_def,
      Sens_A = Sens_A_def,
      Sens_B = Sens_B_def,
      Sens_C = Sens_C_def,
      Sens_D = Sens_D_def,
      Sym_0 = Sym_0_def,
      Sym_A = Sym_A_def,
      Sym_B = Sym_B_def,
      Sym_C = Sym_C_def,
      Sym_D = Sym_D_def
    ),
    parent = globalenv()
  )
  extract_assigned(us_theor_pa, "tmat", env)
}

# The 15x15 matrix built inside neg_log_lik()
nll_matrix <- function(probs = hypo_probs) {
  env <- list2env(list(hypo_probs = probs), parent = globalenv())
  extract_assigned(neg_log_lik, "matrix_of_trans", env)
}
