# Model core (file R/model_core.R)
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
# Notation follows the paper (Tables 2-3, Fig. 1): p_*, Sens_*, Sym_*, L_*,
# p_NA, p_CA, p_Att. Short local names (vec, Dx, sym, us) have no paper symbol.

# ----- Constants --------------------------------------------------------------

# Sym_B = SYM_B_RATIO * Sym_A: symptomatic detection of BCLC B is assumed eight
# times that of BCLC A. Source: Cumming et al. 2024, Supplement Sec. 2.1.3,
# footnote 4 (midpoint of the Sym_A = 0 and Sym_B = 0 MLEs).
SYM_B_RATIO <- 8

# ----- Parameter handling -----------------------------------------------------

#' Override named parameters, erroring on any name that is not a parameter.
#'
#' This is what stops a typo such as `p_att` from silently returning the
#' baseline (finding C2).
#'
#' Ref: none (code safeguard, not part of the publication's model).
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
#' Ref: Fleurence & Hollenbeak 2007, PharmacoEconomics 25(1):3-6; applied to
#' p_HCC, the only annual probability (Cumming 2024, Table 2).
#'
#' @param p  Annual probability (scalar or vector).
#' @return Monthly probability, same shape as `p`.
annual_to_monthly <- function(p) {
  # Ref: Fleurence & Hollenbeak 2007 (rate-to-probability conversion)
  1 - (1 - p)^(1 / 12)
}

#' Build the 7x7 monthly transition matrix.
#'
#' Cirrhosis -> BCLC 0 with probability `p_HCC_month`; each BCLC stage
#' progresses to the next with its own probability; BCLC D -> Mortality;
#' Mortality is absorbing.
#'
#' Ref: Cumming 2024, Sec. 2.1 and Fig. 1 (base model of undiagnosed HCC;
#' arrows p_HCC, p_0, p_A, p_B, p_C, p_D).
#'
#' @param p_HCC_month  Monthly probability of developing HCC from cirrhosis.
#' @param p_stage      Monthly progression probabilities (p_0, p_A, p_B, p_C,
#'                     p_D), in that order.
#' @return 7x7 row-stochastic matrix with dimnames `stage_list`.
transition_matrix <- function(p_HCC_month, p_stage) {
  stopifnot(length(p_stage) == 5)
  # Ref: Cumming 2024, Fig. 1 (arrow probabilities out of each state)
  p <- c(p_HCC_month, unname(p_stage)) # progression out of states 1..6
  tm <- matrix(0, 7, 7, dimnames = list(stage_list, stage_list))
  # Ref: Cumming 2024, Fig. 1 (stay with probability 1 - p)
  tm[cbind(1:6, 1:6)] <- 1 - p # stay in state
  # Ref: Cumming 2024, Fig. 1 (move to the next state with probability p)
  tm[cbind(1:6, 2:7)] <- p # move to the next state
  tm[7, 7] <- 1 # Mortality absorbing; Ref: Cumming 2024, Fig. 1 (Death)
  tm
}

# ----- Base progression model -------------------------------------------------

#' Markov progression of an untreated cohort.
#'
#' Ref: Cumming 2024, Sec. 2.1 and Fig. 1 (base model of undiagnosed HCC).
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
    annual_to_monthly(params$p_HCC), # p_HCC is annual (Table 2)
    unlist(params[STAGE_P_NAMES]) # p_0..p_D are already monthly (Table 1)
  )
  # Ref: Cumming 2024, Sec. 2.1 (monthly steps; all start in Cirrhosis).
  # Chapman-Kolmogorov: t monthly steps.
  stage <- c(n, rep(0, 6)) %*% (tmat %^% t)
  as_tibble(stage, .name_repair = ~stage_list)
}


# ----- Ultrasound detection model (by attendance probability) -----------------

#' N x k matrix of the named parameters, recycling length-1 entries to N rows.
#'
#' Ref: none (helper for the vectorised engine).
#'
#' @param params  Named list; each element has length 1 or N.
#' @param nms     Names of the parameters to extract.
#' @param N       Number of rows.
#' @return N x length(nms) numeric matrix.
param_matrix <- function(params, nms, N) {
  matrix(unlist(lapply(params[nms], rep_len, N)), nrow = N)
}

#' Number of parallel cohorts implied by the (possibly vector) inputs.
#'
#' Ref: none (helper for the vectorised engine).
#'
#' @param params  Named parameter list (each element length 1 or N).
#' @param ...     Further vectors (`Pr_Att`, `n`) of length 1 or N.
#' @return Integer N (the longest input length).
n_cohorts <- function(params, ...) {
  max(lengths(params), lengths(list(...)))
}

#' Vectorised ultrasound detection model: N cohorts at once.
#'
#' Same model and same order of operations as the scalar loop (each month:
#' progress, symptomatic detection, and every `m` months attended ultrasound),
#' but on an N x 7 state matrix with element-wise column updates instead of
#' `vec %*% tmat` (the transition matrix is almost all zeros).
#'
#' Ref: Cumming 2024, Fig. 1 (treatment layer: symptoms or incidental
#' diagnosis with probability Sym_X, else 6-monthly ultrasound that is
#' attended with probability p_Att and detects with probability Sens_X) and
#' Sec. 2.1 (p. 4, scheduled surveillance).
#'
#' @param Pr_Att  Attendance probability, length 1 or N.
#' @param params  Named parameter list (like `params_def`); each element has
#'                length 1 or N.
#' @param n       Cohort size, length 1 or N.
#' @param t       Total follow-up in months.
#' @param m       Ultrasound interval in months.
#' @return List with `dx` and `undiagnosed`, both N x 7 matrices with
#'   colnames `stage_list`.
screen_cohort_mat <- function(Pr_Att, params, n, t, m) {
  N <- n_cohorts(params, Pr_Att, n)
  # Ref: Cumming 2024, Fig. 1 (arrow probabilities; p_HCC is annual)
  p <- cbind(
    annual_to_monthly(rep_len(params$p_HCC, N)),
    param_matrix(params, STAGE_P_NAMES, N)
  )
  stay <- 1 - p # same expression as transition_matrix()
  # Ref: Cumming 2024, Fig. 1 and Tables 2-3 (Sym_X; Cirrhosis 0, Death 1)
  sym_prob <- cbind(0, param_matrix(params, STAGE_SYM_NAMES, N), 1)
  # Ref: Cumming 2024, Fig. 1 and Tables 2-3 (Sens_X)
  US_probs <- cbind(0, param_matrix(params, STAGE_SENS_NAMES, N), 1)
  # Ref: Cumming 2024, Fig. 1 (attend ultrasound: p_NA, p_CA, p_Att)
  attend <- cbind(0, matrix(rep_len(Pr_Att, N), N, 5), 1)
  vec <- cbind(rep_len(n, N), matrix(0, N, 6)) # undiagnosed per state
  Dx <- matrix(0, N, 7)
  for (i in 1:t) {
    # Ref: Cumming 2024, Fig. 1 ("1 month", progress to next BCLC stage)
    vec <- progress_step(vec, p, stay)
    # Ref: Cumming 2024, Fig. 1 (symptoms or incidental diagnosis, Sym_X)
    sym <- vec * sym_prob
    Dx <- Dx + sym
    vec <- vec - sym
    # Ref: Cumming 2024, Fig. 1 (6 months since last scheduled ultrasound?)
    if (i %% m == 0) {
      # Ref: Cumming 2024, Fig. 1 (attend ultrasound? HCC detected? Sens_X)
      us <- (vec * attend) * US_probs
      Dx <- Dx + us
      vec <- vec - us
    }
  }
  colnames(Dx) <- colnames(vec) <- stage_list
  list(dx = Dx, undiagnosed = vec)
}

#' Vectorised "due-until-attend" detection model: N cohorts at once.
#'
#' Alternative to `screen_cohort_mat()` for the 3-class adherence model of
#' paper 2 (not in the publication). There is no fixed scan calendar. A person
#' is either not due (they attended a scan fewer than `m` months ago) or due,
#' and a due person attends with probability `q` EACH MONTH until they attend.
#' So `q` is a per-due-month probability, not a per-scheduled-scan one as
#' `Pr_Att` is in `screen_cohort_mat()`.
#'
#' State: a list of `m` N x 7 matrices (undiagnosed people per disease state)
#' indexed by clock slot. At the attendance step of a month, slot k < m holds
#' people who attended k months earlier (not due) and slot m holds the due
#' people. After a scan in month s the person is therefore next able to attend
#' in month s + m, as in `due_exposure()` (R/adherence_data.R).
#'
#' Each month, in the same order as `screen_cohort_mat()`: (1) progression on
#' every slot, (2) symptomatic detection (Sym_X) on every slot, (3) due people
#' attend with probability `q` and attenders are detected with Sens_X (the
#' Cirrhosis column has Sens 0 and the Mortality column 1, as in
#' `screen_cohort_mat()`); attended-but-undetected people go to slot 1, other
#' not-due people move one slot on and due non-attenders stay due.
#'
#' `init = "gap"`: at t = 0 everyone has just been scanned (slot 1 at month 1),
#' so the first possible scan is in month m. With q = 1 the scan months are
#' m, 2m, 3m, ..., exactly those of `screen_cohort_mat()` with Pr_Att = 1.
#' `init = "due"`: everyone is due in month 1. Assumption (off by one): the
#' data fit (`due_exposure()`) counts month 0 as the first due month, whereas
#' the first model step is month 1; a person with exposure T months in the data
#' corresponds to t = T model months, month i being data month i - 1.
#'
#' Ref: plan 2026-11-24, Sec. 4 S3 (not in the publication); the monthly
#' sequence follows Cumming 2024, Fig. 1.
#'
#' @param q       Monthly attendance probability while due, length 1 or N.
#' @param params  Named parameter list; each element has length 1 or N.
#' @param n       Cohort size, length 1 or N.
#' @param t       Total follow-up in months.
#' @param m       Months from an attended scan to the next possible one.
#' @param init    "due" (everyone due in month 1) or "gap" (everyone just
#'                scanned at month 0).
#' @return List with `dx` and `undiagnosed` (N x 7 matrices, colnames
#'   `stage_list`; `undiagnosed` is summed over clock slots) and `scans`
#'   (length N, expected number of attended scans per cohort).
screen_cohort_due_mat <- function(q, params, n, t, m, init = c("due", "gap")) {
  init <- match.arg(init)
  N <- n_cohorts(params, q, n)
  p <- cbind(
    annual_to_monthly(rep_len(params$p_HCC, N)),
    param_matrix(params, STAGE_P_NAMES, N)
  )
  stay <- 1 - p
  sym_prob <- cbind(0, param_matrix(params, STAGE_SYM_NAMES, N), 1)
  US_probs <- cbind(0, param_matrix(params, STAGE_SENS_NAMES, N), 1)
  qv <- rep_len(q, N) # row-wise: multiplies each column of an N x 7 matrix
  slots <- replicate(m, matrix(0, N, 7), simplify = FALSE)
  slots[[if (init == "due") m else 1]] <- cbind(rep_len(n, N), matrix(0, N, 6))
  Dx <- matrix(0, N, 7)
  scans <- numeric(N)
  for (i in 1:t) {
    for (k in seq_len(m)) {
      vec <- progress_step(slots[[k]], p, stay)
      sym <- vec * sym_prob
      Dx <- Dx + sym
      slots[[k]] <- vec - sym
    }
    due <- slots[[m]]
    att <- due * qv
    us <- att * US_probs
    Dx <- Dx + us
    scans <- scans + rowSums(att)
    new <- replicate(m, matrix(0, N, 7), simplify = FALSE)
    new[[1]] <- att - us
    for (k in seq_len(m - 1)) {
      new[[k + 1]] <- slots[[k]]
    }
    new[[m]] <- new[[m]] + (due - att)
    slots <- new
  }
  vec <- Reduce(`+`, slots)
  colnames(Dx) <- colnames(vec) <- stage_list
  list(dx = Dx, undiagnosed = vec, scans = scans)
}

#' One monthly step of the progression model on an N x 7 state matrix.
#'
#' Equivalent to `vec %*% transition_matrix(...)`, with terms added in the same
#' order (k = 1..7).
#'
#' Ref: Cumming 2024, Sec. 2.1 and Fig. 1 (arrows p_HCC, p_0, ..., p_D).
#'
#' @param vec   N x 7 undiagnosed people per state.
#' @param p     N x 6 progression probabilities out of states 1..6.
#' @param stay  `1 - p`.
#' @return N x 7 matrix after one month.
progress_step <- function(vec, p, stay) {
  new <- vec
  new[, 1] <- vec[, 1] * stay[, 1]
  for (j in 2:6) {
    new[, j] <- vec[, j - 1] * p[, j - 1] + vec[, j] * stay[, j]
  }
  new[, 7] <- vec[, 6] * p[, 6] + vec[, 7]
  new
}

#' Theoretical ultrasound detection for a cohort with a given attendance prob.
#'
#' Each month the cohort moves one step through the transition matrix, then
#' symptomatic detection removes `Sym_*` of each state, and every `m` months
#' attending people are screened with sensitivity `Sens_*`. Thin wrapper
#' around `screen_cohort_mat()` with N = 1.
#'
#' Ref: Cumming 2024, Fig. 1 (treatment layer: symptoms or incidental
#' diagnosis with probability Sym_X, else 6-monthly ultrasound that is
#' attended with probability p_Att and detects with probability Sens_X) and
#' Sec. 2.1 (p. 4, scheduled surveillance).
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
screen_cohort <- function(Pr_Att, ..., n = 100, t = 120, m = 6, params = params_def) {
  params <- modify_params(params, list(...))
  res <- screen_cohort_mat(Pr_Att, params, n, t, m)
  list(
    dx = setNames(as.numeric(res$dx), stage_list),
    undiagnosed = setNames(as.numeric(res$undiagnosed), stage_list)
  )
}


# ----- Adherence-stratified detection model -----------------------------------

#' Class shares and monthly attendance probabilities for adherence = "due".
#'
#' Requires every element of ADH_PI_NAMES and ADH_ATT_NAMES in `params` (never
#' reads p_NA, p_CA or p_Att) and checks the shares sum to 1 (tolerance 1e-8).
#'
#' Ref: plan 2026-11-24 (paper 2); not in the publication.
#'
#' @param params  Named parameter list.
#' @return list(pi = shares Low, Med, High; q = attendance probabilities).
due_class_params <- function(params) {
  miss <- setdiff(c(ADH_PI_NAMES, ADH_ATT_NAMES), names(params))
  if (length(miss) > 0) {
    stop("adherence = \"due\" requires parameter(s): ", paste(miss, collapse = ", "))
  }
  pi <- lapply(params[ADH_PI_NAMES], identity)
  tot <- Reduce(`+`, pi)
  if (any(abs(tot - 1) > 1e-8)) {
    stop("Pi_Low + Pi_Med + Pi_High must equal 1 (tolerance 1e-8)")
  }
  list(pi = unname(pi), q = unname(params[ADH_ATT_NAMES]))
}

#' Vectorised adherence-stratified diagnoses: N parameter sets at once.
#'
#' Stacks the three adherence groups as 3N rows, runs the cohort engine once,
#' then adds the groups in the order non- + partial + complete.
#'
#' `adherence = "scheduled"` (default, the publication) uses
#' `screen_cohort_mat()` with attendance 0, `p_Att`, 1 per scheduled scan;
#' `Pi_*`, `Att_*` and `init` are ignored (and not required).
#' `adherence = "due"` (paper 2, not in the publication) uses
#' `screen_cohort_due_mat()`: the groups attend each month while due with
#' probability `Att_Low`, `Att_Med`, `Att_High`, with class shares `Pi_Low`,
#' `Pi_Med`, `Pi_High` (sum 1); `p_NA`, `p_CA`, `p_Att` are not read. All six
#' are required in `params` (error otherwise).
#'
#' Ref: Cumming 2024, Sec. 2.1 (p. 4: completely adherent, partially adherent
#' and non-adherent cohorts), Table 3 (p_NA, p_CA, p_Att); Supplement
#' Sec. 3.1 (Pr(partially adherent) = 1 - Pr(non-) - Pr(completely adherent)).
#' "due": plan 2026-11-24 / not in the publication.
#'
#' @param params     Named parameter list; each element has length 1 or N.
#' @param n          Cohort size (length 1 or N).
#' @param t          Follow-up months.
#' @param m          Ultrasound interval months.
#' @param adherence  "scheduled" or "due" (see above).
#' @param init       "due" or "gap"; only used when `adherence = "due"`, see
#'                   `screen_cohort_due_mat()`.
#' @return N x 7 matrix (people diagnosed per stage), colnames `stage_list`.
screen_by_adherence_mat <- function(
  params,
  n,
  t,
  m,
  adherence = c("scheduled", "due"),
  init = c("due", "gap")
) {
  adherence <- match.arg(adherence)
  init <- match.arg(init)
  due <- if (adherence == "due") due_class_params(params)
  N <- n_cohorts(params, n)
  if (adherence == "due") {
    n_grp <- c(
      n * rep_len(due$pi[[1]], N),
      n * rep_len(due$pi[[2]], N),
      n * rep_len(due$pi[[3]], N)
    )
  } else {
    p_NA <- rep_len(params$p_NA, N)
    p_CA <- rep_len(params$p_CA, N)
    # Ref: Supplement Sec. 3.1 (p. 11, partially adherent = 1 - p_NA - p_CA)
    n_grp <- c(n * p_NA, n * (1 - p_NA - p_CA), n * p_CA)
  }
  stacked <- lapply(params, function(x) rep(rep_len(x, N), 3))
  if (adherence == "due") {
    q <- c(rep_len(due$q[[1]], N), rep_len(due$q[[2]], N), rep_len(due$q[[3]], N))
    dx <- screen_cohort_due_mat(q, stacked, n_grp, t, m, init)$dx
  } else {
    # Ref: Cumming 2024, Sec. 2.1 (attendance 0, p_Att and 1 per cohort)
    att <- c(rep(0, N), rep_len(params$p_Att, N), rep(1, N))
    dx <- screen_cohort_mat(att, stacked, n_grp, t, m)$dx
  }
  rows <- split(seq_len(3 * N), rep(1:3, each = N))
  (dx[rows[[1]], , drop = FALSE] + dx[rows[[2]], , drop = FALSE]) +
    dx[rows[[3]], , drop = FALSE]
}

#' Combine detection across three adherence groups (non-, partial, complete).
#'
#' Thin wrapper around `screen_by_adherence_mat()` with N = 1.
#'
#' Ref: Cumming 2024, Sec. 2.1 (p. 4: completely adherent, partially adherent
#' and non-adherent cohorts), Table 3 (p_NA, p_CA, p_Att); Supplement
#' Sec. 3.1 (Pr(partially adherent) = 1 - Pr(non-) - Pr(completely adherent)).
#'
#' @param ...     Parameter overrides by name (e.g. `p_NA`, `p_CA`, `p_Att`,
#'                `Sens_0`); unknown names are an error.
#' @param n       Cohort size.
#' @param t       Follow-up months.
#' @param m       Ultrasound interval months.
#' @param params  Named parameter list (default `params_def`).
#' @param adherence  "scheduled" (default) or "due"; see
#'                `screen_by_adherence_mat()`. "due" needs `Pi_*` and `Att_*`
#'                in `params` (or as overrides, which only work if `params`
#'                contains them).
#' @param init    "due" or "gap" (only for `adherence = "due"`).
#' @return Tibble with one column per `stage_list` state (people diagnosed).
screen_by_adherence <- function(
  ...,
  n = 100,
  t = 120,
  m = 6,
  params = params_def,
  adherence = c("scheduled", "due"),
  init = c("due", "gap")
) {
  adherence <- match.arg(adherence)
  init <- match.arg(init)
  params <- modify_params(params, list(...))
  dx <- screen_by_adherence_mat(params, n, t, m, adherence, init)
  as_tibble(dx, .name_repair = "minimal")
}


# ----- Life-years added by surveillance and treatment -------------------------

#' Vectorised life years added by the surveillance programme: N draws at once.
#'
#' Columns are accumulated left to right (as `sum()` does) so the result
#' matches the scalar calculation.
#'
#' Ref: Cumming 2024, Sec. 2.1 (p. 5: diagnoses per stage multiplied by the
#' additional life expectancy L_X) and Table 2 (L_0 .. L_D).
#'
#' @param params  Named parameter list; each element has length 1 or N.
#' @param n       Cohort size (length 1 or N).
#' @param t       Follow-up months.
#' @param m       Ultrasound interval months.
#' @param adherence  "scheduled" (default) or "due"; see
#'                `screen_by_adherence_mat()`.
#' @param init    "due" or "gap" (only for `adherence = "due"`).
#' @return Numeric vector of length N (life-years per cohort of n).
life_add_vec <- function(
  params,
  n,
  t,
  m,
  adherence = c("scheduled", "due"),
  init = c("due", "gap")
) {
  dx <- screen_by_adherence_mat(params, n, t, m, match.arg(adherence), match.arg(init))
  # Ref: Cumming 2024, Sec. 2.1 (p. 5), Table 2 (L_X)
  Lmat <- cbind(0, param_matrix(params, STAGE_L_NAMES, nrow(dx)), 0)
  s <- 0
  for (k in 1:7) {
    s <- s + dx[, k] * Lmat[, k]
  }
  unname(s)
}

#' Total life years added by the HCC surveillance programme.
#'
#' Sums the product of stage-specific diagnoses and stage-specific life-year
#' gains (L_*) across all adherence groups. Thin wrapper around
#' `life_add_vec()` with N = 1.
#'
#' Ref: Cumming 2024, Sec. 2.1 (p. 5: diagnoses per stage multiplied by the
#' additional life expectancy L_X) and Table 2 (L_0 .. L_D).
#'
#' @param ...     Parameter overrides by name (any of `names(params_def)`);
#'                unknown names are an error.
#' @param n       Cohort size.
#' @param t       Follow-up months.
#' @param m       Ultrasound interval months.
#' @param params  Named parameter list (default `params_def`).
#' @param adherence  "scheduled" (default) or "due"; see
#'                `screen_by_adherence_mat()`.
#' @param init    "due" or "gap" (only for `adherence = "due"`).
#' @return Scalar life-years gained per cohort of n.
life_add <- function(
  ...,
  n = 100,
  t = 120,
  m = 6,
  params = params_def,
  adherence = c("scheduled", "due"),
  init = c("due", "gap")
) {
  params <- modify_params(params, list(...))
  life_add_vec(params, n, t, m, match.arg(adherence), match.arg(init))
}


# ----- Sensitivity-analysis helpers -------------------------------------------

#' Build a one-row tibble of named parameter values for use with pmap.
#'
#' Ref: none (helper for the univariate sensitivity analysis, Sec. 2.3).
#'
#' @param Parameter  Parameter name (length 1).
#' @param value      Its value.
#' @param ...        Ignored (pmap passes the other columns of the table).
#' @return One-row tibble with one column, named `Parameter`.
dfs <- function(Parameter, value, ...) {
  tibble::as_tibble_row(setNames(value, Parameter))
}

#' Keep only the model-parameter columns of a draws tibble.
#'
#' Ref: none (code safeguard).
#'
#' Bootstrap tibbles carry scenario/result columns next to the parameter
#' columns; life_add() rejects unknown names (finding C2), so every
#' pmap(..., life_add) call selects the parameter columns explicitly.
#'
#' @param df  Tibble of parameter draws (one row per draw).
#' @return `df` restricted to columns named in `names(params_def)`.
draws <- function(df) select(df, any_of(names(params_def)))

#' Parameter draws ready for the model: the parameter columns of `df` with
#' `Sym_B = SYM_B_RATIO * Sym_A` set from each draw's Sym_A.
#'
#' Ref: Supplement Sec. 2.1.3, footnote 4 (8 x Sym_A = Sym_B).
#'
#' Sym_B is not estimated (Supplement Sec. 2.1.3, footnote 4), so every
#' bootstrap analysis, the scenarios and every sweep, must derive it from the
#' draw's Sym_A rather than keep the default Sym_B. This is the one place that
#' does it. If `df` has no Sym_A column, the default Sym_A is used.
#'
#' @param df  Tibble of parameter draws (one row per draw).
#' @return `draws(df)` with a `Sym_B` column (replacing any existing one).
prepare_draws <- function(df) {
  out <- draws(df)
  sym_a <- if ("Sym_A" %in% names(out)) out$Sym_A else params_def$Sym_A
  out$Sym_B <- SYM_B_RATIO * sym_a # Ref: Supplement Sec. 2.1.3, footnote 4
  out
}

#' Apply life_add across each row of a data-frame of parameter draws.
#'
#' Ref: none (applies `life_add()`, Sec. 2.1, to each bootstrap draw).
#'
#' @param df  Tibble of parameter draws (columns named as `params_def`).
#' @return Numeric vector of life-years, one per row.
years <- function(df) pmap_dbl(df, life_add)

#' Sum of absolute values.
#'
#' Ref: none (summary for the sensitivity plots, Fig. 2D).
#'
#' @param x  Numeric vector.
#' @return Scalar sum of `abs(x)`.
abssum <- function(x) sum(abs(x))

#' Build the parameter tibble used for univariate sensitivity tornado plots.
#'
#' Ref: Cumming 2024, Sec. 2.3 and Fig. 2D (one-at-a-time variation of the
#' modifiable parameters; +/-25% bounds when `csv` is Parameters_25.csv) and
#' Sec. 2.5 (low and high p_HCC).
#'
#' @param csv  Full path to a Parameters CSV, e.g. data/Parameters.csv.
#' @return Tibble with one row per parameter × bound combination plus computed
#'   life-year columns.
life_tibble <- function(csv) {
  lower <- load_params(csv, "LOWER")
  upper <- load_params(csv, "UPPER")
  param_tidy <- tibble(
    Parameter = names(lower), # parameters by name, in CSV order
    LOWER = unlist(lower, use.names = FALSE),
    UPPER = unlist(upper, use.names = FALSE)
  ) |>
    pivot_longer(-Parameter, names_to = "Bound")
  param_tidy |>
    mutate(
      dfs_new = pmap(param_tidy, dfs),
      Bound = as.factor(str_to_title(Bound)),
      Parameter = factor(
        Parameter,
        levels = str_sort(unique(Parameter), locale = "en")
      ),
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
      # Ref: Cumming 2024, Sec. 2.3 (change in life years vs the point
      # estimate), Sec. 2.5 (p_HCC low / high)
      years_mid_comp = years_mid - life_add(),
      years_low_comp = years_low - life_add(p_HCC = lower$p_HCC),
      years_high_comp = years_high - life_add(p_HCC = upper$p_HCC)
    )
}
