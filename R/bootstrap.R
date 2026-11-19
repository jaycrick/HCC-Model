# R/bootstrap.R
# Functions used in the bootstrap parameter estimation and uncertainty analysis.
# Requires: R/setup.R, R/model_core.R, and the hypogeom and DescTools packages.
# Ref: Cumming 2024, Sec. 2.2 and Supplement Sec. 2 (calibration, bootstrap).
# Paper notation (Supp. Table S1): m_X median months, n_X patients, k, D_0,
# D_A, D_alpha, s-bar. Here: `med`/`meds`, `n`, `n_early`, `Dx_0`, `Dx_A`,
# `Dx_sym_early`, `n_sym`. These are argument names of neg_log_lik(), used by
# tests and scripts, so they are not renamed.

# ----- Calibration data -------------------------------------------------------

# Names and order of the 8 free parameters of the likelihood: the five monthly
# stage-progression probabilities, ultrasound sensitivity in BCLC 0 and A, and
# symptomatic detection in BCLC A (Sym_B = SYM_B_RATIO * Sym_A).
NLL_PAR_NAMES <- c(STAGE_P_NAMES, "Sens_0", "Sens_A", "Sym_A")

#' Observed calibration data, defined once (Cumming et al. 2024, Cancers
#' 16:2745, and its Supplement, p. 3 and footnote 2).
#'
#' One row per BCLC stage, in the order 0, A, B, C, D:
#'   n      patients observed in the stage;
#'   med    observed median months spent in the stage;
#'   mlen   length of the order-statistic sample passed to dOrder(): n if n is
#'          odd, n - 1 if n is even (Supp. p.3 footnote 2: the median is then
#'          the middle value of an odd-length sample). So 12 -> 11,
#'          138 -> 137, 210 -> 209; 101 and 139 are unchanged. This is the
#'          intended rule, not an off-by-one (finding G6, closed);
#'   j      the median's rank in that sample, (mlen + 1) / 2.
#' Ref: Supplement Sec. 2.1.1 (p. 2-3: n_X = (12, 101, 139, 138, 210),
#' m_X = (38, 25, 10, 7, 4), (k, D_0, D_A) = (109, 36, 51), s-bar = 158).
#' `dx` holds the diagnosis counts: `Dx_0` of `n_early` early-stage cases were
#' found at BCLC 0 and `Dx_A` of the remaining `n_early - Dx_0` at BCLC A,
#' and `Dx_sym_early` of `n_sym` symptomatic cases were found at BCLC 0 or A.
#'
#' @format A list: `stage` (data frame: stage, n, med, mlen, j) and `dx`
#'   (named numeric: Dx_0, Dx_A, Dx_sym_early, n_early, n_sym).
calibration_data <- local({
  stage <- data.frame(
    stage = c("0", "A", "B", "C", "D"),
    n = c(12, 101, 139, 138, 210),
    med = c(38, 25, 10, 7, 4)
  )
  # Ref: Supplement p. 3, footnote 2 (even n: use n - 1 observations)
  stage$mlen <- ifelse(stage$n %% 2 == 1, stage$n, stage$n - 1)
  # Ref: Supplement p. 3, footnote 2 (median = middle order statistic)
  stage$j <- (stage$mlen + 1) / 2
  list(
    stage = stage,
    # Finding P3: paper Sec. 2.2.2 says "52 of n = 156"; the supplement
    # (Sec. 2.1.1, p. 3: s-bar = 158, D_alpha = 127 - 75 = 52) and the code
    # use 158, which reproduces the published MLE.
    dx = c(Dx_0 = 36, Dx_A = 51, Dx_sym_early = 52, n_early = 109, n_sym = 158)
  )
})

# ----- Transition matrix ------------------------------------------------------

#' The 15 x 15 monthly transition matrix of Supp. matrix S1 / Table S2.
#'
#' Ref: Supplement Table S2 (p. 5, states) and matrix (S1) (p. 16). Finding
#' P5 (note): Table S2 lists 14 states; the 15th, absorbing `past_A`, is the
#' code's "progressed past BCLC A undiagnosed" state.
#'
#' States (rows = from, columns = to), named in the dimnames:
#'   dx_0, dx_A         absorbing: diagnosed at BCLC 0 / A;
#'   u0_6 ... u0_1      undiagnosed BCLC 0 with 6 ... 1 months to the next
#'                      ultrasound (u0_1 is screened on the next step);
#'   uA_6 ... uA_1      the same for undiagnosed BCLC A;
#'   past_A             absorbing: progressed past BCLC A undiagnosed.
#' # uncertain: whether "k months to ultrasound" counts the current month; the
#' # code and Supp. Table S2 only fix the transitions, which are used here.
#'
#' @param p_0,p_A     Monthly progression probabilities from BCLC 0 and A.
#' @param sens_0,sens_A  Ultrasound sensitivity in BCLC 0 and A.
#' @return 15 x 15 numeric matrix with named dimensions.
likelihood_matrix <- function(p_0, p_A, sens_0, sens_A) {
  states <- c(
    "dx_0",
    "dx_A",
    paste0("u0_", 6:1),
    paste0("uA_", 6:1),
    "past_A"
  )
  tm <- matrix(0, 15, 15, dimnames = list(states, states))
  # Ref: Supplement matrix (S1), p. 16 (diagnosed states are absorbing)
  tm["dx_0", "dx_0"] <- 1
  tm["dx_A", "dx_A"] <- 1
  # Finding P5: 15th state, not in Supplement Table S2 (14 states)
  tm["past_A", "past_A"] <- 1
  # Ref: Supplement matrix (S1), p. 16 (entries 1 - p_0 and p_0)
  for (k in 6:2) {
    # a month passes: stay in stage 0 (one month closer to ultrasound) or
    # progress to stage A at the same point in the screening cycle
    tm[paste0("u0_", k), paste0("u0_", k - 1)] <- 1 - p_0
    tm[paste0("u0_", k), paste0("uA_", k - 1)] <- p_0
  }
  # screening month from BCLC 0: detected, or restart the 6-month cycle
  # Ref: Supplement matrix (S1), p. 16 ((1-p_0)Sens_0, p_0 Sens_A,
  # (1-p_0)(1-Sens_0), p_0(1-Sens_A))
  tm["u0_1", "dx_0"] <- (1 - p_0) * sens_0
  tm["u0_1", "dx_A"] <- p_0 * sens_A
  tm["u0_1", "u0_6"] <- (1 - p_0) * (1 - sens_0)
  tm["u0_1", "uA_6"] <- p_0 * (1 - sens_A)
  # Ref: Supplement matrix (S1), p. 16 (entries 1 - p_A and p_A)
  for (k in 6:2) {
    tm[paste0("uA_", k), paste0("uA_", k - 1)] <- 1 - p_A
    tm[paste0("uA_", k), "past_A"] <- p_A
  }
  # Ref: Supplement matrix (S1), p. 16 ((1-p_A)Sens_A, (1-p_A)(1-Sens_A), p_A)
  tm["uA_1", "dx_A"] <- (1 - p_A) * sens_A
  tm["uA_1", "uA_6"] <- (1 - p_A) * (1 - sens_A)
  tm["uA_1", "past_A"] <- p_A
  tm
}

#' Probability of ending each state after `t` months, starting from the
#' uniform position in the 6-month screening cycle (Supp. footnote 3).
#'
#' Ref: Supplement Sec. 2.1.3 (p. 4: initial distribution t^(0), limiting
#' distribution approximated by t^(100)) and footnote 3 (p. 5, convergence
#' check against the exact solution for Pr(d_0)).
#'
#' Pr(d_0) is the first element. The closed form is
#' s0 (1-(1-p0)^6) / (6 p0 (1-(1-s0)(1-p0)^6)); t = 100 months stands in for
#' infinity (error about 1.5e-6, tested against the closed form).
#'
#' @param tm  Matrix from likelihood_matrix().
#' @param sens_0  Ultrasound sensitivity in BCLC 0 (sets the start vector).
#' @param t   Number of months.
#' @return Numeric vector of length 15 (names as `tm`).
end_state_probs <- function(tm, sens_0, t = 100) {
  # Ref: Supplement p. 4 (initial distribution t^(0))
  start <- c(sens_0 / 6, 0, (1 - sens_0) / 6, rep(1 / 6, 5), rep(0, 7))
  # Ref: Supplement p. 4 (t^(100) = t^(0) T^100)
  drop(start %*% (tm %^% t))
}

# ----- Negative log-likelihood for parametric bootstrap ----------------------

#' Negative log-likelihood for joint parameter estimation.
#'
#' Combines hypergeometric-order-statistics likelihoods for median survival
#' times with binomial likelihoods for diagnosis proportions.
#'
#' Ref: Supplement Sec. 2.1.3 (p. 3-6): L = L_m0 ... L_mD x L_D0 x L_DA|D0 x
#' L_Dalpha; minimised by Nelder-Mead (p. 6).
#'
#' @param par  Length-8 numeric vector of candidate probabilities, in the
#'             order of NLL_PAR_NAMES.
#' @param meds  Observed median months per stage (0, A, B, C, D).
#' @param Dx_0,Dx_A,Dx_sym_early  Observed diagnosis counts.
#' @param n_early,n_sym  Their denominators (see calibration_data).
#' @return Scalar NLL, or NaN if any probability is outside [0,1].
neg_log_lik <- function(
  par,
  meds = calibration_data$stage$med,
  Dx_0 = calibration_data$dx[["Dx_0"]],
  Dx_A = calibration_data$dx[["Dx_A"]],
  Dx_sym_early = calibration_data$dx[["Dx_sym_early"]],
  n_early = calibration_data$dx[["n_early"]],
  n_sym = calibration_data$dx[["n_sym"]]
) {
  if (any(par < 0 | par > 1)) {
    return(NaN)
  }

  p_A <- par[2]
  p_B <- par[3]
  Sym_A <- par[8]
  Sym_B <- SYM_B_RATIO * Sym_A # Ref: Supplement Sec. 2.1.3, footnote 4

  tm <- likelihood_matrix(par[1], p_A, par[6], par[7])
  t100 <- end_state_probs(tm, par[6]) # t^(100), Supplement p. 4

  # Order-statistic part: the median of the time in stage k is the sum of the
  # times in stages k..D, each at least one month, so its distribution is
  # hypogeometric in par[k:5] after subtracting (6 - k) months.
  # Ref: Supplement Sec. 2.1.2-2.1.3 (p. 3: HypoGeom(p_D, p_C, ..., p_X),
  # median order statistic via DescTools::dOrder).
  cal <- calibration_data$stage
  ll_median <- 0
  for (k in 5:1) {
    ll_median <- ll_median +
      DescTools::dOrder(
        meds[k] - (6 - k),
        distn = "hypogeom",
        probs = par[k:5],
        mlen = cal$mlen[k],
        j = cal$j[k],
        log = TRUE
      )
  }

  # Share of symptomatic cases found at BCLC A, then at B given not A.
  # Supp. p.6 prints Pr(d_alpha) without the (1 - Sym) factors; that is a
  # typo. This form, Sym / (1 - (1-p)(1-Sym)), is the correct one: it gives
  # 0.3289 at the MLE against 52/158 = 0.329 observed.
  # Ref: Supplement p. 6, Pr(d_alpha) (finding P1: printed form is a typo)
  pr_sym_A <- Sym_A / (1 - (1 - p_A) * (1 - Sym_A))
  pr_sym_B <- Sym_B / (1 - (1 - p_B) * (1 - Sym_B))
  pr_d_alpha <- pr_sym_A + (1 - pr_sym_A) * pr_sym_B

  # Ref: Supplement p. 4-5 (L_D0 = Binom(D_0; k, t_1), L_DA|D0 =
  # Binom(D_A; k - D_0, t_2 / (1 - t_1))) and p. 6 (L_Dalpha).
  # Finding P5: t_1 uses the t^(100) approximation, not the exact closed form
  # (difference about 1.5e-6); kept because it produced the published MLE.
  -(ll_median +
    dbinom(Dx_0, n_early, t100[1], log = TRUE) +
    dbinom(
      Dx_A,
      n_early - Dx_0,
      t100[2] / (1 - t100[1]),
      log = TRUE
    ) +
    dbinom(Dx_sym_early, n_sym, pr_d_alpha, log = TRUE))
}

#' Maximum-likelihood point estimate by Nelder-Mead.
#' The published MLE is NM's stopping point at default tolerance (CHANGES P7).
#'
#' Ref: Supplement Sec. 2.1.3 (p. 6: "numerically solved using the
#' Nelder-Mead method") and Cumming 2024, Sec. 2.2.1.
#'
#' @param start  Starting values (length 8, see NLL_PAR_NAMES).
#' @return The optim() result; `$par` is named by NLL_PAR_NAMES.
fit_point_estimate <- function(start) {
  fit <- stats::optim(
    start,
    neg_log_lik,
    method = "Nelder-Mead",
    control = list(maxit = 10000)
  )
  fit$par <- setNames(fit$par, NLL_PAR_NAMES)
  fit
}

# ----- Bootstrap outputs (shipped RDS inputs) ---------------------------------

#' Bootstrapped calibration inputs (m_*, D_*) in long form, from the HPC
#' list-of-lists output.
#'
#' Ref: Supplement Sec. 2.1.4 (p. 7, bootstrap samples m_X*b, D_0*b, D_A*b,
#' D_alpha*b) and Fig. S1 (p. 8).
#'
#' @param path  Path to the shipped RDS (list of resampled inputs).
#' @return Long tibble with columns `name` and `value`.
bootstrap_inputs_long <- function(
  path = here::here("data-derived", "list_of_resampled_milton_output.rds")
) {
  nms <- c("m_0", "m_A", "m_B", "m_C", "m_D", "D_0", "D_A", "D_alpha")
  tibble(all = readRDS(path)) |>
    unnest_wider(all, names_sep = "Yes", names_repair = ~nms) |>
    suppressMessages() |>
    pivot_longer(all_of(nms))
}

#' The 10 000 bootstrap MLEs (8 columns, NLL_PAR_NAMES) from the HPC output.
#'
#' Ref: Supplement Sec. 2.1.5 (p. 7, bootstrap maximum likelihood estimates)
#' and Fig. S2 (p. 9); Cumming 2024, Sec. 2.2.1 (10,000 resamples).
#'
#' @param path  Path to the shipped RDS (HPC optimiser output).
#' @return Tibble with one row per resample, columns `NLL_PAR_NAMES`.
bootstrap_estimates <- function(
  path = here::here("data-derived", "miltonoutput.rds")
) {
  tibble(list_of_pars = readRDS(path)) |>
    unnest_wider(list_of_pars) |>
    select(par) |>
    unnest_wider(par, names_sep = "k", names_repair = ~NLL_PAR_NAMES) |>
    suppressMessages()
}

# Scenario definitions live in R/scenarios.R (run_scenarios()).

# ----- LaTeX label helper ----------------------------------------------------

#' Convert a variable name to its TeX expression label.
#'
#' @param value  Parameter name(s), e.g. "p_HCC".
#' @return The matching expression(s) from `TeX_labels`.
Tex_labeller <- function(value) {
  return(TeX_labels[value])
}
