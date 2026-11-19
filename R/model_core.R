# R/model_core.R
# Core Markov model functions: parameter handling, transition matrix,
# untreated progression, ultrasound detection, adherence model, life-years
# calculation, and sensitivity-analysis helpers.
# All functions are pure (no file I/O, no set.seed) -- they depend only on
# objects initialised in R/setup.R (`stage_list`, `params_def`).
#
# Model structure: monthly discrete-time Markov cohort over the seven states
# in `stage_list` (Cirrhosis -> BCLC 0 -> A -> B -> C -> D -> Mortality);
# Cumming et al. 2024, Cancers 16:2745, Sec. 2.1 and Fig. 1.
# p_HCC is an ANNUAL probability; p_0..p_D, Sens_*, Sym_* are MONTHLY or
# per-event probabilities (paper Tables 1-3).

# ----- Constants --------------------------------------------------------------

# Sym_B = SYM_B_RATIO * Sym_A: symptomatic detection of BCLC B is assumed eight
# times that of BCLC A. Source: Cumming et al. 2024, Supplement Sec. 2.1.3,
# footnote 4 (midpoint of the Sym_A = 0 and Sym_B = 0 MLEs).
SYM_B_RATIO <- 8

# Names of the five monthly stage-to-stage progression probabilities, in the
# order of the transition matrix rows BCLC 0 .. BCLC D.
STAGE_P_NAMES <- c("p_0", "p_A", "p_B", "p_C", "p_D")

# ----- Parameter handling -----------------------------------------------------

#' Override named parameters, erroring on any name that is not a parameter.
#'
#' This is what stops a typo such as `p_att` from silently returning the
#' baseline (finding C2).
#'
#' @param params     Named list of parameters (e.g. `params_def`).
#' @param overrides  Named list of replacement values.
#' @return `params` with the overrides applied.
modify_params <- function(params, overrides) {
  if (length(overrides) == 0) {
    return(params)
  }
  nms <- names(overrides)
  if (is.null(nms) || any(nms == "")) {
    stop("unnamed parameter override(s); all overrides must be named")
  }
  unknown <- setdiff(nms, names(params))
  if (length(unknown) > 0) {
    stop("unknown parameter(s): ", paste(unknown, collapse = ", "))
  }
  params[nms] <- overrides
  params
}

# ----- Probability conversion and transition matrix ---------------------------

#' Convert an annual probability to the equivalent monthly probability.
#'
#' Assumes a constant hazard: p_month = 1 - (1 - p_annual)^(1/12).
#' Fleurence & Hollenbeak 2007, PharmacoEconomics 25(1):3-6.
#'
#' @param p  Annual probability (scalar or vector).
#' @return Monthly probability, same shape as `p`.
annual_to_monthly <- function(p) {
  1 - (1 - p)^(1 / 12)
}

#' Build the 7x7 monthly transition matrix.
#'
#' Cirrhosis -> BCLC 0 with probability `p_HCC_month`; each BCLC stage
#' progresses to the next with its own probability; BCLC D -> Mortality;
#' Mortality is absorbing. Cumming et al. 2024 Cancers 16:2745, Sec. 2.1, Fig. 1.
#'
#' @param p_HCC_month  Monthly probability of developing HCC from cirrhosis.
#' @param p_stage      Monthly progression probabilities (p_0, p_A, p_B, p_C,
#'                     p_D), in that order.
#' @return 7x7 row-stochastic matrix with dimnames `stage_list`.
transition_matrix <- function(p_HCC_month, p_stage) {
  stopifnot(length(p_stage) == 5)
  p <- c(p_HCC_month, unname(p_stage)) # progression out of states 1..6
  tm <- matrix(0, 7, 7, dimnames = list(stage_list, stage_list))
  tm[cbind(1:6, 1:6)] <- 1 - p # stay in state
  tm[cbind(1:6, 2:7)] <- p # move to the next state
  tm[7, 7] <- 1 # Mortality absorbing
  tm
}

# ----- Base progression model -------------------------------------------------

#' Markov progression of an untreated cohort.
#'
#' @param ...     Parameter overrides by name (e.g. `p_HCC = 0.03`); unknown
#'                names are an error.
#' @param n       Starting cohort size (all in Cirrhosis at t = 0).
#' @param t       Number of months to run.
#' @param params  Named parameter list (default `params_def`).
#' @return Tibble with one column per `stage_list` state (people per stage).
progress_untreated <- function(..., n = 100, t = 120, params = params_def) {
  params <- modify_params(params, list(...))
  tmat <- transition_matrix(
    annual_to_monthly(params$p_HCC), # p_HCC is annual
    unlist(params[STAGE_P_NAMES]) # p_0..p_D are already monthly (paper Table 1)
  )
  stage <- c(n, rep(0, 6)) %*% (tmat %^% t) # Chapman-Kolmogorov, t steps
  as_tibble(stage, .name_repair = ~stage_list)
}


# ----- Ultrasound detection model (by attendance probability) -----------------

#' Theoretical ultrasound detection for a cohort with a given attendance prob.
#'
#' Each month the cohort moves one step through the transition matrix, then
#' symptomatic detection removes `Sym_*` of each state, and every `m` months
#' attending people are screened with sensitivity `Sens_*`.
#' (Cumming et al. 2024 Cancers 16:2745, Sec. 2.1-2.2.)
#'
#' @param Pr_Att  Probability of attending a scheduled ultrasound (0 = none,
#'                1 = complete adherence).
#' @param ...     Parameter overrides by name; unknown names are an error.
#' @param n       Cohort size.
#' @param t       Total follow-up in months.
#' @param m       Ultrasound interval in months.
#' @param params  Named parameter list (default `params_def`).
#' @return List with `dx` (people diagnosed per stage, length 7) and
#'   `undiagnosed` (people still undiagnosed in each state at month `t`,
#'   length 7); both named by `stage_list`, and `sum(dx) + sum(undiagnosed)`
#'   equals `n`.
screen_cohort <- function(
  Pr_Att,
  ...,
  n = 100,
  t = 120,
  m = 6,
  params = params_def
) {
  params <- modify_params(params, list(...))
  tmat <- transition_matrix(
    annual_to_monthly(params$p_HCC), # p_HCC is annual
    unlist(params[STAGE_P_NAMES])
  )
  vec <- c(n, rep(0, 6)) # undiagnosed people per state
  Dx <- rep(0, 7) # diagnosed people
  sym_prob <- c(
    0,
    unlist(params[paste0("Sym_", c("0", "A", "B", "C", "D"))]),
    1
  )
  US_probs <- c(
    0,
    unlist(params[paste0("Sens_", c("0", "A", "B", "C", "D"))]),
    1
  )
  attend <- c(0, rep(Pr_Att, 5), 1) # who goes to the ultrasound
  for (i in 1:t) {
    vec <- vec %*% tmat
    sym <- vec * sym_prob # detected through symptoms
    Dx <- Dx + sym
    vec <- vec - sym
    if (i %% m == 0) {
      go <- vec * attend # how many actually go to the ultrasound
      us <- go * US_probs # detected by ultrasound
      Dx <- Dx + us
      vec <- vec - us
    }
  }
  list(
    dx = setNames(as.numeric(Dx), stage_list),
    undiagnosed = setNames(as.numeric(vec), stage_list)
  )
}


# ----- Adherence-stratified detection model -----------------------------------

#' Combine detection across three adherence groups (non-, partial, complete).
#'
#' @param ...     Parameter overrides by name (e.g. `p_NA`, `p_CA`, `p_Att`,
#'                `Sens_0`); unknown names are an error.
#' @param n       Cohort size.
#' @param t       Follow-up months.
#' @param m       Ultrasound interval months.
#' @param params  Named parameter list (default `params_def`).
#' @return Tibble with one column per `stage_list` state (people diagnosed).
screen_by_adherence <- function(
  ...,
  n = 100,
  t = 120,
  m = 6,
  params = params_def
) {
  params <- modify_params(params, list(...))
  # people in each adherence group: non-, partially, completely adherent
  adh_group <- n * c(params$p_NA, 1 - params$p_NA - params$p_CA, params$p_CA)

  caught <- function(n_group, Pr_Att) {
    dx <- screen_cohort(Pr_Att, n = n_group, t = t, m = m, params = params)$dx
    as_tibble(as.list(dx))
  }

  map2_dfr(adh_group, c(0, params$p_Att, 1), caught) %>%
    summarise_all(sum)
}


# ----- Life-years added by surveillance and treatment -------------------------

#' Total life years added by the HCC surveillance programme.
#'
#' Sums the product of stage-specific diagnoses and stage-specific life-year
#' gains (L_*) across all adherence groups.
#'
#' @param ...     Parameter overrides by name (any of `names(params_def)`);
#'                unknown names are an error.
#' @param n       Cohort size.
#' @param t       Follow-up months.
#' @param m       Ultrasound interval months.
#' @param params  Named parameter list (default `params_def`).
#' @return Scalar life-years gained per cohort of n.
life_add <- function(..., n = 100, t = 120, m = 6, params = params_def) {
  params <- modify_params(params, list(...))
  stage_caught <- screen_by_adherence(n = n, t = t, m = m, params = params)
  L <- unname(unlist(params[paste0("L_", c("0", "A", "B", "C", "D"))]))
  sum(stage_caught * c(0, L, 0)) # life-years: sum_stage Dx_stage * L_stage
}


# ----- Sensitivity-analysis helpers -------------------------------------------

#' Build a one-row tibble of named parameter values for use with pmap.
dfs <- function(Parameter, value, ...) {
  new_df <- tibble(Parameter, value) %>%
    pivot_wider(names_from = Parameter, values_from = value)
  return(new_df)
}

#' Keep only the model-parameter columns of a draws tibble.
#'
#' Bootstrap tibbles carry scenario/result columns next to the parameter
#' columns; life_add() rejects unknown names (finding C2), so every
#' pmap(..., life_add) call selects the parameter columns explicitly.
#'
#' @param df  Tibble of parameter draws (one row per draw).
#' @return `df` restricted to columns named in `names(params_def)`.
draws <- function(df) select(df, any_of(names(params_def)))

#' Apply life_add across each row of a data-frame of parameter draws.
years <- function(df) pmap_dbl(df, life_add)

#' Sum of absolute values.
abssum <- function(x) sum(abs(x))

#' Build the parameter tibble used for univariate sensitivity tornado plots.
#'
#' @param csv  Full path to a Parameters CSV (e.g. here("data","Parameters.csv")).
#' @return Tibble with one row per parameter × bound combination plus computed
#'   life-year columns.
life_tibble <- function(csv) {
  lower <- load_params(csv, "LOWER")
  upper <- load_params(csv, "UPPER")
  param_tidy <- tibble(
    Parameter = names(lower), # parameters by name, in CSV order
    LOWER = unlist(lower, use.names = FALSE),
    UPPER = unlist(upper, use.names = FALSE)
  ) %>%
    pivot_longer(-Parameter, names_to = "Bound")
  param_tbl <- param_tidy %>%
    mutate(
      dfs_new = pmap(., dfs),
      Bound = as.factor(str_to_title(Bound)),
      Parameter = as.factor(Parameter),
      dfs_low = map(
        dfs_new,
        ~ if (names(.x) == "p_HCC") .x else tibble(.x, p_HCC = lower$p_HCC)
      ),
      dfs_high = map(
        dfs_new,
        ~ if (names(.x) == "p_HCC") .x else tibble(.x, p_HCC = upper$p_HCC)
      ),
      years_mid = map_dbl(dfs_new, years),
      years_low = map_dbl(dfs_low, years),
      years_high = map_dbl(dfs_high, years),
      years_mid_comp = years_mid - life_add(),
      years_low_comp = years_low - life_add(p_HCC = lower$p_HCC),
      years_high_comp = years_high - life_add(p_HCC = upper$p_HCC)
    )
  return(param_tbl)
}
