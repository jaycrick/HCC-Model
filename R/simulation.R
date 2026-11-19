# R/simulation.R
# Stochastic (Monte Carlo) check of the deterministic cohort model.
#
# These functions replay, person by person (as binomial counts per state), the
# monthly sequence of operations in screen_cohort() / screen_by_adherence()
# (R/model_core.R): progress -> symptomatic detection -> ultrasound on
# months with i %% m == 0. Their expectation therefore equals the
# deterministic model's output, and tests/testthat/test-simulation.R checks
# that. Not used in the published analysis (Cumming et al. 2024, Cancers
# 16(15):2745), which uses the deterministic model.
# Like model_core.R, nothing here uses globals or sets a seed.

#' Stochastic analogue of screen_cohort().
#'
#' @inheritParams screen_cohort
#' @return Numeric vector of length 7 (named by `stage_list`): people
#'   diagnosed per stage.
simulate_cohort <- function(
  Pr_Att,
  ...,
  n = 1000,
  t = 120,
  m = 6,
  params = params_def
) {
  params <- modify_params(params, list(...))
  tmat <- transition_matrix(
    annual_to_monthly(params$p_HCC),
    unlist(params[STAGE_P_NAMES])
  )
  p_move <- tmat[cbind(1:6, 2:7)] # progression probability of states 1..6
  sym_prob <- c(
    0,
    unlist(params[paste0("Sym_", c("0", "A", "B", "C", "D"))]),
    1
  )
  us_prob <- c(
    0,
    unlist(params[paste0("Sens_", c("0", "A", "B", "C", "D"))]),
    1
  )
  attend <- c(0, rep(Pr_Att, 5), 1)

  vec <- c(n, rep(0, 6)) # undiagnosed people per state
  dx <- rep(0, 7)
  for (i in 1:t) {
    move <- rbinom(6, vec[1:6], p_move) # all draws use start-of-month counts
    vec <- vec - c(move, 0) + c(0, move)
    sym <- rbinom(7, vec, sym_prob)
    dx <- dx + sym
    vec <- vec - sym
    if (i %% m == 0) {
      go <- rbinom(7, vec, attend)
      us <- rbinom(7, go, us_prob)
      dx <- dx + us
      vec <- vec - us
    }
  }
  setNames(dx, stage_list)
}

#' Stochastic analogue of screen_by_adherence().
#'
#' People are split into non-, partially and completely adherent groups
#' (multinomial draw), each group is simulated with `simulate_cohort()`
#' (attendance 0, `p_Att`, 1) and the diagnoses are summed.
#'
#' @inheritParams screen_by_adherence
#' @return Numeric vector of length 7: people diagnosed per stage.
simulate_by_adherence <- function(
  ...,
  n = 1000,
  t = 120,
  m = 6,
  params = params_def
) {
  params <- modify_params(params, list(...))
  prob <- c(params$p_NA, 1 - params$p_NA - params$p_CA, params$p_CA)
  n_group <- as.vector(rmultinom(1, n, prob))
  att <- c(0, params$p_Att, 1)
  rowSums(mapply(
    function(ng, a) simulate_cohort(a, n = ng, t = t, m = m, params = params),
    n_group,
    att
  ))
}
