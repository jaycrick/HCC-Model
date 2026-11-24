# Stochastic simulation (file R/simulation.R)
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
# Ref: not in the publication; mirrors Cumming 2024, Sec. 2.1 and Fig. 1.

#' Stochastic analogue of screen_cohort().
#'
#' Ref: Cumming 2024, Fig. 1 (same monthly sequence as `screen_cohort()`;
#' binomial draws instead of expected counts).
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
  # Ref: Cumming 2024, Fig. 1 (arrow probabilities p_HCC, p_0 .. p_D)
  p_move <- tmat[cbind(1:6, 2:7)] # progression probability of states 1..6
  sym_prob <- c(0, unlist(params[STAGE_SYM_NAMES]), 1)
  us_prob <- c(0, unlist(params[STAGE_SENS_NAMES]), 1)
  # Ref: Cumming 2024, Fig. 1 (attend ultrasound?)
  attend <- c(0, rep(Pr_Att, 5), 1)

  vec <- c(n, rep(0, 6)) # undiagnosed people per state
  dx <- rep(0, 7)
  for (i in 1:t) {
    # Ref: Cumming 2024, Fig. 1 (progress to next BCLC stage)
    move <- rbinom(6, vec[1:6], p_move) # all draws use start-of-month counts
    vec <- vec - c(move, 0) + c(0, move)
    # Ref: Cumming 2024, Fig. 1 (symptoms or incidental diagnosis, Sym_X)
    sym <- rbinom(7, vec, sym_prob)
    dx <- dx + sym
    vec <- vec - sym
    if (i %% m == 0) {
      # Ref: Cumming 2024, Fig. 1 (attend ultrasound, then Sens_X)
      go <- rbinom(7, vec, attend)
      us <- rbinom(7, go, us_prob)
      dx <- dx + us
      vec <- vec - us
    }
  }
  setNames(dx, stage_list)
}

#' Stochastic analogue of screen_cohort_due_mat() (one cohort).
#'
#' Same slots and order of operations as `screen_cohort_due_mat()`: per clock
#' slot, binomial progression then symptomatic detection; then the due slot
#' attends (binomial with `q`) and attenders are detected (binomial with
#' Sens_X). See that function for the clock and `init` semantics (`q` is per
#' due month).
#'
#' Ref: plan 2026-11-24, Sec. 4 S3 (not in the publication).
#'
#' @inheritParams screen_cohort
#' @param q     Monthly attendance probability while due (scalar).
#' @param init  "due" or "gap".
#' @return Numeric vector of length 7 (named by `stage_list`): people
#'   diagnosed per stage.
simulate_cohort_due <- function(
  q,
  ...,
  n = 1000,
  t = 120,
  m = 6,
  params = params_def,
  init = c("due", "gap")
) {
  init <- match.arg(init)
  params <- modify_params(params, list(...))
  tmat <- transition_matrix(
    annual_to_monthly(params$p_HCC),
    unlist(params[STAGE_P_NAMES])
  )
  p_move <- tmat[cbind(1:6, 2:7)]
  sym_prob <- c(0, unlist(params[STAGE_SYM_NAMES]), 1)
  us_prob <- c(0, unlist(params[STAGE_SENS_NAMES]), 1)
  slots <- replicate(m, numeric(7), simplify = FALSE)
  slots[[if (init == "due") m else 1]] <- c(n, rep(0, 6))
  dx <- rep(0, 7)
  for (i in 1:t) {
    for (k in seq_len(m)) {
      vec <- slots[[k]]
      move <- rbinom(6, vec[1:6], p_move)
      vec <- vec - c(move, 0) + c(0, move)
      sym <- rbinom(7, vec, sym_prob)
      dx <- dx + sym
      slots[[k]] <- vec - sym
    }
    due <- slots[[m]]
    att <- rbinom(7, due, q)
    us <- rbinom(7, att, us_prob)
    dx <- dx + us
    new <- replicate(m, numeric(7), simplify = FALSE)
    new[[1]] <- att - us
    for (k in seq_len(m - 1)) {
      new[[k + 1]] <- slots[[k]]
    }
    new[[m]] <- new[[m]] + (due - att)
    slots <- new
  }
  setNames(dx, stage_list)
}

#' Stochastic analogue of screen_by_adherence().
#'
#' Ref: Cumming 2024, Sec. 2.1 (p. 4, three adherence cohorts) and Table 3.
#'
#' People are split into non-, partially and completely adherent groups
#' (multinomial draw), each group is simulated with `simulate_cohort()`
#' (attendance 0, `p_Att`, 1) and the diagnoses are summed.
#'
#' With `adherence = "due"` (paper 2, not in the publication) each group is
#' simulated with `simulate_cohort_due()` and monthly due-attendance
#' probabilities `Att_Low`, `Att_Med`, `Att_High` and class shares `Pi_*` (required;
#' error if missing; p_NA, p_CA, p_Att are not read); `Att_*` are per due month.
#' `init` is only used in that mode.
#'
#' @inheritParams screen_by_adherence
#' @return Numeric vector of length 7: people diagnosed per stage.
simulate_by_adherence <- function(
  ...,
  n = 1000,
  t = 120,
  m = 6,
  params = params_def,
  adherence = c("scheduled", "due"),
  init = c("due", "gap")
) {
  adherence <- match.arg(adherence)
  init <- match.arg(init)
  params <- modify_params(params, list(...))
  due <- if (adherence == "due") due_class_params(params)
  # Ref: Supplement Sec. 3.1 (p. 11, partially adherent = 1 - p_NA - p_CA)
  prob <- if (adherence == "due") {
    unlist(due$pi)
  } else {
    c(params$p_NA, 1 - params$p_NA - params$p_CA, params$p_CA)
  }
  n_group <- as.vector(rmultinom(1, n, prob))
  # Ref: Cumming 2024, Sec. 2.1 (attendance 0, p_Att, 1 per cohort)
  if (adherence == "due") {
    q <- unlist(due$q)
    return(rowSums(mapply(
      function(ng, a) {
        simulate_cohort_due(a, n = ng, t = t, m = m, params = params, init = init)
      },
      n_group,
      q
    )))
  }
  att <- c(0, params$p_Att, 1)
  rowSums(mapply(
    function(ng, a) simulate_cohort(a, n = ng, t = t, m = m, params = params),
    n_group,
    att
  ))
}
