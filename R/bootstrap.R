# Calibration and bootstrap (file R/bootstrap.R)
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

# Nelder-Mead starting values for the 8 likelihood parameters (NLL_PAR_NAMES).
NM_START <- c(0.08, 0.06, 0.34, 0.35, 0.16, 0.22, 0.54, 0.01)

# The published MLE (Cumming 2024, Table 3; the POINT column of
# data/Parameters.csv up to tag v1.1.0, before the G9 refit). Second starting
# point of the refit's point fit (tools/refit_calibration.R). It is a fixed
# constant, not read from Parameters.csv, because the refit rewrites that CSV:
# starting from the CSV made a rerun start from the previous run's output and
# moved 253 of the 10 000 bootstrap fits by up to 2e-5 (relative).
PUBLISHED_MLE <- c(
  p_0 = 0.07982013,
  p_A = 0.05540597,
  p_B = 0.34353183,
  p_C = 0.35197652,
  p_D = 0.16046211,
  Sens_0 = 0.20916773,
  Sens_A = 0.53267066,
  Sym_A = 0.01015362
)

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
    stage = STAGE_IDS,
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

#' Probability of ending in each state, starting from the uniform position in
#' the 6-month screening cycle (Supp. footnote 3).
#'
#' Ref: Supplement Sec. 2.1.3 (p. 4: initial distribution t^(0)) and footnote 3
#' (p. 5: absorbing chain, fundamental matrix N = (I - Q)^-1). Finding P5 is
#' applied (gate G8): with `t = NULL` the absorption probabilities are exact,
#' B = (I - Q)^-1 R, instead of the published t^(100) = t^(0) T^100
#' approximation (error about 1.4e-6 in Pr(d_0)).
#'
#' Transient states are 3:14 (u0_6 ... uA_1); absorbing states are 1, 2, 15
#' (dx_0, dx_A, past_A). Pr(d_0) is the first element. Its closed form is
#' s0 (1-(1-p0)^6) / (6 p0 (1-(1-s0)(1-p0)^6)).
#'
#' @param tm  Matrix from likelihood_matrix().
#' @param sens_0  Ultrasound sensitivity in BCLC 0 (sets the start vector).
#' @param t   `NULL` (exact, default) or a number of months, which gives the
#'   old power approximation (kept for the convergence test).
#' @return Numeric vector of length 15 (names as `tm`); transient entries are
#'   0 when `t` is `NULL`.
end_state_probs <- function(tm, sens_0, t = NULL) {
  # Ref: Supplement p. 4 (initial distribution t^(0))
  start <- c(sens_0 / 6, 0, (1 - sens_0) / 6, rep(1 / 6, 5), rep(0, 7))
  if (!is.null(t)) {
    # Ref: Supplement p. 4 (t^(t) = t^(0) T^t)
    return(drop(start %*% (tm %^% t)))
  }
  tr <- 3:14
  ab <- c(1, 2, 15)
  # Ref: Supplement footnote 3 (absorbing chain; B = (I - Q)^-1 R)
  B <- solve(diag(12) - tm[tr, tr], tm[tr, ab])
  out <- setNames(numeric(15), rownames(tm))
  out[ab] <- start[ab] + drop(start[tr] %*% B)
  out
}

# ----- Negative log-likelihood for parametric bootstrap ----------------------

#' Order-statistic part of the log-likelihood: the observed median months in
#' stages 0, A, B, C and D under the hypergeometric model.
#'
#' The median of the time in stage k is the sum of the times in stages k..D,
#' each at least one month, so its distribution is hypogeometric in par[k:5]
#' after subtracting (6 - k) months.
#' Ref: Supplement Sec. 2.1.2-2.1.3 (p. 3: HypoGeom(p_D, p_C, ..., p_X),
#' median order statistic via DescTools::dOrder).
#'
#' @param par   Likelihood parameter vector (see NLL_PAR_NAMES).
#' @param meds  Observed median months per stage (0, A, B, C, D).
#' @return Scalar log-likelihood.
loglik_medians <- function(par, meds) {
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
  ll_median
}

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
  absorb <- end_state_probs(tm, par[6]) # exact, Supplement fn 3 (P5, G8)

  ll_median <- loglik_medians(par, meds)

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
  # Finding P5 (applied, G8): t_1 and t_2 are the exact absorption
  # probabilities, not the published t^(100) approximation (about 1.4e-6).
  -(ll_median +
    dbinom(Dx_0, n_early, absorb[1], log = TRUE) +
    dbinom(
      Dx_A,
      n_early - Dx_0,
      absorb[2] / (1 - absorb[1]),
      log = TRUE
    ) +
    dbinom(Dx_sym_early, n_sym, pr_d_alpha, log = TRUE))
}

# Box for the polish step (probability scale). Reasons:
#  H1: a tiny p makes the hypogeometric FFT support huge (up to 8 s and many
#      GB per evaluation at p = 1e-4), so the five p_* stay above 1e-3;
#  H2: Sym_B = SYM_B_RATIO * Sym_A must be <= 1, so Sym_A < 1 / SYM_B_RATIO;
#  H4: fits can sit at p = 1, so the upper bound is 1 - 1e-6, not 1.
FIT_LB <- setNames(c(rep(1e-3, 5), 1e-6, 1e-6, 1e-6), NLL_PAR_NAMES)
FIT_UB <- setNames(
  c(rep(1 - 1e-6, 7), (1 / SYM_B_RATIO) * (1 - 1e-6)),
  NLL_PAR_NAMES
)
# Penalty for a non-finite NLL inside the polish only (H3: Inf at p_D near 1).
FIT_PENALTY <- 1e10

#' Polish a Nelder-Mead result by L-BFGS-B on the scaled-logit scale
#' z = qlogis(par / ub), with `par` clipped into the box first.
#'
#' Ref: Supplement Sec. 2.1.3 (p. 6) uses Nelder-Mead alone; the polish is the
#' P7 improvement (gate G9).
#'
#' @param par  NM result (length 8).
#' @param ... Data arguments passed to neg_log_lik().
#' @return The optim() result with `$par` on the probability scale, or NULL if
#'   the polish errors.
polish_fit <- function(par, ...) {
  par <- pmin(pmax(par, FIT_LB), FIT_UB)
  to_par <- function(z) FIT_UB * plogis(z)
  obj <- function(z) {
    v <- neg_log_lik(to_par(z), ...)
    if (is.finite(v)) v else FIT_PENALTY
  }
  tryCatch(
    {
      fit <- stats::optim(
        qlogis(par / FIT_UB),
        obj,
        method = "L-BFGS-B",
        lower = qlogis(FIT_LB / FIT_UB),
        upper = qlogis(1 - 1e-6),
        control = list(factr = 1e3, maxit = 1000)
      )
      fit$par <- to_par(fit$par)
      fit
    },
    error = function(e) NULL
  )
}

#' Stage 1 of one start: Nelder-Mead, or (polish-only) the clipped start.
#'
#' @param start  Numeric vector (NM + polish) or `list(par = <numeric>, nm =
#'   FALSE)` (polish only; `par` is clipped into the box).
#' @param ... Data arguments passed to neg_log_lik().
#' @return An optim()-like list (`par`, `value`, `convergence`, `message`).
start_stage <- function(start, ...) {
  if (!is.list(start)) {
    start <- list(par = start, nm = TRUE)
  }
  if (isFALSE(start$nm)) {
    par <- pmin(pmax(start$par, FIT_LB), FIT_UB)
    return(list(
      par = par,
      value = neg_log_lik(par, ...),
      convergence = 99L,
      message = "not optimised (nm = FALSE)",
      nm_run = FALSE
    ))
  }
  nm <- stats::optim(
    start$par,
    neg_log_lik,
    ...,
    method = "Nelder-Mead",
    control = list(maxit = 10000)
  )
  nm$nm_run <- TRUE
  nm
}

#' One start: Nelder-Mead (unless `nm = FALSE`), then the polish; keep the
#' better of the two.
#'
#' @param start  See start_stage().
#' @param ... Data arguments passed to neg_log_lik().
#' @return List: `fit` (chosen optim result), `nm_value` (NA if Nelder-Mead
#'   was skipped), `polish_value`, `polish_convergence`.
fit_one_start <- function(start, ...) {
  nm <- start_stage(start, ...)
  pol <- polish_fit(nm$par, ...)
  use_pol <- !is.null(pol) && !(pol$value >= nm$value)
  list(
    fit = if (use_pol) pol else nm,
    nm_value = if (nm$nm_run) nm$value else NA_real_,
    polish_value = if (is.null(pol)) NA_real_ else pol$value,
    polish_convergence = if (is.null(pol)) NA_integer_ else pol$convergence
  )
}

#' Which parameters sit on the box? One-sided: on or outside the lower or
#' upper limit, within a relative 1e-6 (the upper tolerance is 2e-6 so that the
#' polish limit, ub * plogis(qlogis(1 - 1e-6)), is flagged too).
#'
#' @param par  Named parameter vector (length 8).
#' @return Named logical vector.
fit_at_bound <- function(par) {
  par <- setNames(par, NLL_PAR_NAMES)
  par <= FIT_LB * (1 + 1e-6) | par >= FIT_UB * (1 - 2e-6)
}

#' Maximum-likelihood estimate: Nelder-Mead then L-BFGS-B polish, multi-start.
#'
#' Ref: Supplement Sec. 2.1.3 (p. 6: "numerically solved using the
#' Nelder-Mead method") and Cumming 2024, Sec. 2.2.1. Finding P7 is applied
#' (gate G9): each start is polished on the scaled-logit scale inside the box
#' FIT_LB..FIT_UB; the polish is kept only if it lowers the NLL.
#'
#' @param starts  List of starts. Each is a numeric vector (length 8, see
#'   NLL_PAR_NAMES; Nelder-Mead then polish) or `list(par = <numeric>, nm =
#'   FALSE)` (polish only, clipped into the box).
#' @param meds,Dx_0,Dx_A,Dx_sym_early,n_early,n_sym  Data, as in neg_log_lik().
#' @return List: `par` (named), `value`, `convergence` (of the chosen stage;
#'   0 = OK; 99 if the only stage was an unpolished polish-only start),
#'   `message`, `at_bound` (named logical, see fit_at_bound()) and
#'   `starts` (data frame: start, nm_value, polish_value, polish_convergence,
#'   chosen).
fit_mle <- function(
  starts,
  meds = calibration_data$stage$med,
  Dx_0 = calibration_data$dx[["Dx_0"]],
  Dx_A = calibration_data$dx[["Dx_A"]],
  Dx_sym_early = calibration_data$dx[["Dx_sym_early"]],
  n_early = calibration_data$dx[["n_early"]],
  n_sym = calibration_data$dx[["n_sym"]]
) {
  res <- lapply(
    starts,
    fit_one_start,
    meds = meds,
    Dx_0 = Dx_0,
    Dx_A = Dx_A,
    Dx_sym_early = Dx_sym_early,
    n_early = n_early,
    n_sym = n_sym
  )
  vals <- vapply(res, function(r) r$fit$value, numeric(1))
  best <- which.min(vals)
  fit <- res[[best]]$fit
  par <- setNames(fit$par, NLL_PAR_NAMES)
  list(
    par = par,
    value = fit$value,
    convergence = fit$convergence,
    message = if (is.null(fit$message)) "" else fit$message,
    at_bound = fit_at_bound(par),
    starts = data.frame(
      start = seq_along(res),
      nm_value = vapply(res, function(r) r$nm_value, numeric(1)),
      polish_value = vapply(res, function(r) r$polish_value, numeric(1)),
      polish_convergence = vapply(
        res,
        function(r) r$polish_convergence,
        numeric(1)
      ),
      chosen = seq_along(res) == best
    )
  )
}

#' Maximum-likelihood point estimate from one start (thin wrapper).
#' P7 is applied (gate G9): see fit_mle().
#'
#' @param start  Starting values (length 8, see NLL_PAR_NAMES).
#' @return The fit_mle() result for the single start.
fit_point_estimate <- function(start) {
  fit_mle(list(start))
}

#' Does a fit need the fallback? Yes if it did not converge, sits on the box,
#' or has no finite value (e.g. the polish failed).
#'
#' @param fit  A fit_mle() result.
#' @return Logical scalar.
needs_fallback <- function(fit) {
  !is.finite(fit$value) || fit$convergence != 0 || any(fit$at_bound)
}

#' Fit one bootstrap resample, with an optional fallback in the same worker.
#'
#' `best_start` codes the start that gave the chosen result: an index into the
#' primary starts list; or, if the fallback ran and was better, the length of
#' the primary list plus the index into the fallback list.
#'
#' @param x  numeric(8): m_0..m_D, D_0, D_A, D_alpha.
#' @param starts,fallback  Lists of starts (see fit_mle()); `fallback` may be
#'   NULL.
#' @return A fit_mle() result with extra elements `fallback` (logical) and
#'   `best_start`.
fit_resample <- function(x, starts, fallback = NULL) {
  fit_x <- function(s) {
    fit_mle(
      s,
      meds = x[1:5],
      Dx_0 = x[6],
      Dx_A = x[7],
      Dx_sym_early = x[8]
    )
  }
  fit <- fit_x(starts)
  fit$best_start <- which(fit$starts$chosen)[1]
  fit$fallback <- FALSE
  if (!is.null(fallback) && needs_fallback(fit)) {
    alt <- fit_x(fallback)
    fit$fallback <- TRUE
    if (alt$value < fit$value || !is.finite(fit$value)) {
      alt$best_start <- length(starts) + which(alt$starts$chosen)[1]
      alt$fallback <- TRUE
      fit <- alt
    }
  }
  fit
}

#' Number of forked workers parallel::mclapply() can use: `cores`, except on
#' Windows, where forking is unavailable and mclapply() refuses more than 1.
#'
#' @param cores  Requested number of workers.
#' @return Integer number of workers.
fork_cores <- function(cores) {
  if (.Platform$OS.type == "windows") 1L else as.integer(cores)
}

#' Refit the MLE for every bootstrap resample, in parallel and deterministic.
#'
#' Ref: Supplement Sec. 2.1.4-2.1.5 (p. 7: bootstrap samples, bootstrap MLEs).
#' Each element of `inputs` is numeric(8): m_0..m_D, D_0, D_A, D_alpha (the
#' order of bootstrap_inputs_long()), mapped to `meds = x[1:5]`, `Dx_0 = x[6]`,
#' `Dx_A = x[7]`, `Dx_sym_early = x[8]`; `n_early` and `n_sym` are unchanged
#' (calibration_data$dx). Resamples are forked one by one
#' (`mc.preschedule = FALSE`, their cost varies a lot); there is no RNG, and
#' the results are in the order of `ids`, whatever `cores` is.
#'
#' @param inputs  List read from list_of_resampled_milton_output.rds.
#' @param starts_fn  Function of the resample index b returning a list of
#'   starts.
#' @param cores  Forked workers; the result does not depend on it. Forking is
#'   unavailable on Windows, so there it always runs serially (fork_cores()).
#' @param fallback_fn  NULL, or a function of b returning a list of starts to
#'   try (and compare) if the primary fit needs the fallback.
#' @param ids  Resample indices to fit (default all).
#' @return List: `pars` (tibble, columns NLL_PAR_NAMES) and `diag` (tibble: b,
#'   value, convergence, best_start, fallback, and one `at_bound_<name>` per
#'   parameter). See fit_resample() for the `best_start` coding.
refit_bootstrap <- function(
  inputs,
  starts_fn,
  cores = max(1L, parallel::detectCores() - 1L),
  fallback_fn = NULL,
  ids = seq_along(inputs)
) {
  one <- function(b) {
    fb <- if (is.null(fallback_fn)) NULL else fallback_fn(b)
    fit_resample(inputs[[b]], starts_fn(b), fb)
  }
  res <- parallel::mclapply(
    ids,
    function(b) try(one(b), silent = TRUE),
    mc.cores = fork_cores(cores),
    mc.preschedule = FALSE
  )
  bad <- which(vapply(res, inherits, logical(1), what = "try-error"))
  if (length(bad) > 0) {
    stop(
      "refit_bootstrap failed for resample ",
      ids[bad[1]],
      ": ",
      as.character(res[[bad[1]]]),
      call. = FALSE
    )
  }
  pars <- do.call(rbind, lapply(res, function(r) r$par))
  at_bound <- do.call(rbind, lapply(res, function(r) r$at_bound))
  colnames(at_bound) <- paste0("at_bound_", NLL_PAR_NAMES)
  list(
    pars = tibble::as_tibble(pars),
    diag = tibble::tibble(
      b = ids,
      value = vapply(res, function(r) r$value, numeric(1)),
      convergence = vapply(res, function(r) r$convergence, numeric(1)),
      best_start = vapply(
        res,
        function(r) as.integer(r$best_start),
        integer(1)
      ),
      fallback = vapply(res, function(r) r$fallback, logical(1))
    ) |>
      dplyr::bind_cols(tibble::as_tibble(at_bound))
  )
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
  path = paper_path("data-derived", "list_of_resampled_milton_output.rds")
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
  path = paper_path("data-derived", "miltonoutput.rds")
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
  TeX_labels[value]
}
