# Adherence likelihood and fit (file R/adherence.R)
# Three-class "due-until-attend" adherence model: person i has a latent class
# c in (Low, Med, High) with probabilities (Pi_Low, Pi_Med = 1 - Pi_Low - Pi_High,
# Pi_High). While
# due, they attend each month with probability q_c = (Att_Low, Att_Med,
# Att_High); after an attended scan they are not due for m months. The data
# enter only through A (scans attended while due) and F (due months without a
# scan) per person, see R/adherence_data.R. Likelihood:
#   l = sum_cells n * log( sum_c pi_c q_c^A (1 - q_c)^F ).
# All functions are pure: no file I/O, no set.seed (random draws come from the
# caller's RNG state or are passed in as arguments). Depends on stats only.
# Ref: plan 2026-11-24 (paper 2, papers/stochastic); no published reference.

# ----- Constants --------------------------------------------------------------

# Parameter names, in the order used by every function here. Classes are in the
# order Low, Med, High: pi = (Pi_Low, 1 - Pi_Low - Pi_High, Pi_High),
# q = (Att_Low, Att_Med,
# Att_High).
# Free parameters (Pi_Med = 1 - Pi_Low - Pi_High is derived). Constants
# ADH_PI_NAMES / ADH_ATT_NAMES come from R/setup.R.
ADH_FREE_NAMES <- c("Pi_Low", "Pi_High", "Att_Low", "Att_Med", "Att_High")
# All six names, in the order of the paper-2 Parameters.csv; this is the order
# of `par` returned by fit_adherence().
ADH_PAR_NAMES <- c(ADH_PI_NAMES, ADH_ATT_NAMES)

#' Full six-name parameter vector from the five free parameters.
#'
#' Pi_Med = 1 - Pi_Low - Pi_High exactly, so the shares sum to 1.
#'
#' Ref: none (helper).
#'
#' @param par  Named numeric vector holding ADH_FREE_NAMES.
#' @return Named numeric vector in ADH_PAR_NAMES order.
adh_full_par <- function(par) {
  c(
    Pi_Low = par[["Pi_Low"]],
    Pi_Med = 1 - par[["Pi_Low"]] - par[["Pi_High"]],
    Pi_High = par[["Pi_High"]],
    Att_Low = par[["Att_Low"]],
    Att_Med = par[["Att_Med"]],
    Att_High = par[["Att_High"]]
  )
}

# ----- Likelihood -------------------------------------------------------------

#' Split a parameter vector into class probabilities and attendance rates.
#'
#' Ref: none (helper).
#'
#' @param par  Named numeric vector or list with the names ADH_FREE_NAMES.
#' @return list(pi = length-3 vector (Low, Med, High), q = length-3 vector).
adh_unpack <- function(par) {
  par <- unlist(par[ADH_FREE_NAMES])
  list(
    pi = c(par[["Pi_Low"]], 1 - par[["Pi_Low"]] - par[["Pi_High"]], par[["Pi_High"]]),
    q = c(par[["Att_Low"]], par[["Att_Med"]], par[["Att_High"]])
  )
}

#' x * log(y) with the convention 0 * log(0) = 0.
#'
#' Ref: none (0^0 = 1 convention for q = 0 or 1 with A = 0 or F = 0).
#'
#' @param x,y  Numeric vectors of equal length (x >= 0).
#' @return Numeric vector.
adh_xlogy <- function(x, y) {
  out <- numeric(length(x))
  nz <- x != 0
  out[nz] <- x[nz] * log(y[nz])
  out
}

#' Cells x 3 matrix of log( pi_c q_c^A (1 - q_c)^F ).
#'
#' Ref: plan 2026-11-24, finding F2.
#'
#' @param pi,q  Length-3 vectors (see adh_unpack()).
#' @param tab   Count table with columns A, F, n.
#' @return Numeric matrix, one row per cell, one column per class.
adh_log_terms <- function(pi, q, tab) {
  lt <- vapply(
    1:3,
    function(k) {
      log(pi[k]) +
        adh_xlogy(tab$A, rep(q[k], nrow(tab))) +
        adh_xlogy(tab$F, rep(1 - q[k], nrow(tab)))
    },
    numeric(nrow(tab))
  )
  matrix(lt, nrow = nrow(tab))
}

#' Row-wise log-sum-exp of a matrix (-Inf rows stay -Inf).
#'
#' Ref: none (numerical helper).
#'
#' @param lt  Numeric matrix.
#' @return Numeric vector, one value per row.
adh_row_lse <- function(lt) {
  mx <- apply(lt, 1, max)
  safe <- ifelse(is.finite(mx), mx, 0)
  ifelse(is.finite(mx), safe + log(rowSums(exp(lt - safe))), -Inf)
}

#' Weighted log-likelihood of the adherence model.
#'
#' Ref: plan 2026-11-24, finding F2 (l = sum n log sum_c pi_c q_c^A (1-q_c)^F).
#'
#' @param par  Named numeric vector or list with the names ADH_FREE_NAMES.
#' @param tab  Count table with columns A, F, n.
#' @return Scalar log-likelihood; -Inf if pi or any q is outside [0, 1].
adherence_loglik <- function(par, tab) {
  u <- adh_unpack(par)
  if (any(!is.finite(c(u$pi, u$q))) || any(u$pi < 0) || any(u$pi > 1)) {
    return(-Inf)
  }
  if (any(u$q < 0) || any(u$q > 1)) {
    return(-Inf)
  }
  ll <- adh_row_lse(adh_log_terms(u$pi, u$q, tab))
  if (any(is.infinite(ll))) {
    return(-Inf)
  }
  sum(tab$n * ll)
}

#' Posterior class probabilities (responsibilities) for each cell.
#'
#' Ref: plan 2026-11-24, Sec. 4 S2 (E-step).
#'
#' @param par  See adherence_loglik().
#' @param tab  Count table with columns A, F, n.
#' @return Matrix, cells x 3 (columns Low, Med, High); rows sum to 1 (NaN where the
#'   cell has probability 0 under `par`).
adherence_posterior <- function(par, tab) {
  u <- adh_unpack(par)
  lt <- adh_log_terms(u$pi, u$q, tab)
  r <- exp(lt - adh_row_lse(lt))
  colnames(r) <- ADH_IDS
  r
}

# ----- Starts, ordering -------------------------------------------------------

#' Halton sequence value (van der Corput in `base`) for index i = 1, 2, ...
#'
#' Ref: Halton 1960, Numer. Math. 2:84-90 (deterministic low-discrepancy points).
#'
#' @param i     Integer index.
#' @param base  Prime base.
#' @return Number in (0, 1).
adh_halton <- function(i, base) {
  f <- 1
  r <- 0
  while (i > 0) {
    f <- f / base
    r <- r + f * (i %% base)
    i <- i %/% base
  }
  r
}

#' Deterministic multi-start grid (no RNG).
#'
#' The first start is a fixed sensible guess; the others are Halton points
#' (log-uniform sorted q's in [0.003, 0.95], pi inside the simplex).
#'
#' Ref: none (deterministic starts, as R/ files must be pure).
#'
#' @param n_starts  Number of starts (>= 1).
#' @return List of named 5-vectors in ADH_FREE_NAMES order.
adherence_starts <- function(n_starts = 20) {
  first <- c(
    Pi_Low = 1 / 3,
    Pi_High = 1 / 3,
    Att_Low = 0.02,
    Att_Med = 0.15,
    Att_High = 0.7
  )
  rest <- lapply(seq_len(max(n_starts - 1, 0)), function(i) {
    u <- vapply(c(2, 3, 5, 7, 11), function(b) adh_halton(i, b), numeric(1))
    p_na <- 0.05 + 0.6 * u[1]
    p_ca <- (0.95 - p_na) * (0.1 + 0.8 * u[2])
    q <- sort(exp(log(0.003) + u[3:5] * (log(0.95) - log(0.003))))
    c(Pi_Low = p_na, Pi_High = p_ca, Att_Low = q[1], Att_Med = q[2], Att_High = q[3])
  })
  c(list(first), rest)[seq_len(n_starts)]
}

#' Names of the five parameters that are pinned by `fix` after validation.
#'
#' Ref: none (helper).
#'
#' @param fix  Named list of fixed parameter values.
#' @return Named numeric vector (possibly empty).
adh_check_fix <- function(fix) {
  fix <- unlist(fix)
  if (length(fix) == 0) {
    return(numeric(0))
  }
  if (is.null(names(fix)) || !all(names(fix) %in% ADH_FREE_NAMES)) {
    stop("`fix` must be named with names in ADH_FREE_NAMES")
  }
  if (any(!is.finite(fix) | fix < 0 | fix > 1)) {
    stop("fixed values must lie in [0, 1]")
  }
  fix
}

#' Impose the class ordering Att_Low <= Att_Med <= Att_High.
#'
#' A class is anchored (keeps its slot) if its q, or its class probability, is
#' fixed: Pi_Low fixes slot Low, Pi_High fixes slot High (the Med share is
#' derived, so it is anchored only when both are fixed). The remaining free
#' classes are sorted by q among the free slots, permuting pi and q together.
#' Ties are left as they are.
#'
#' Ref: plan 2026-11-24, finding F3 (label switching).
#'
#' @param par       Named 5-vector.
#' @param fix_names Names of the fixed parameters.
#' @return `par` with free classes ordered.
adh_order_classes <- function(par, fix_names = character(0)) {
  u <- adh_unpack(par)
  anchored <- c(
    "Pi_Low" %in% fix_names || all(c("Pi_Low", "Pi_High") %in% fix_names),
    all(c("Pi_Low", "Pi_High") %in% fix_names),
    "Pi_High" %in% fix_names || all(c("Pi_Low", "Pi_High") %in% fix_names)
  ) |
    c("Att_Low", "Att_Med", "Att_High") %in% fix_names
  slots <- which(!anchored)
  if (length(slots) > 1) {
    ord <- slots[order(u$q[slots])]
    u$pi[slots] <- u$pi[ord]
    u$q[slots] <- u$q[ord]
  }
  c(
    Pi_Low = u$pi[1],
    Pi_High = u$pi[3],
    Att_Low = u$q[1],
    Att_Med = u$q[2],
    Att_High = u$q[3]
  )
}

# ----- EM ---------------------------------------------------------------------

#' One EM run from one start.
#'
#' E-step: responsibilities r. M-step: pi_c = sum n r_c / sum n; q_c = sum n r_c
#' A / sum n r_c (A + F) for the free q's. Stops when the log-likelihood
#' changes by less than `tol`. Fixed class probabilities are not supported
#' (use method = "optim").
#'
#' Ref: Dempster, Laird & Rubin 1977, J R Stat Soc B 39:1-38.
#'
#' @param start  Named 5-vector.
#' @param tab    Count table.
#' @param fixv   Named numeric vector of fixed parameters.
#' @param tol,maxit  Convergence tolerance and iteration cap.
#' @return list(par, loglik, iterations, converged).
adh_em_one <- function(start, tab, fixv, tol, maxit) {
  par <- start
  par[names(fixv)] <- fixv
  u <- adh_unpack(par)
  pi <- u$pi
  q <- u$q
  free_q <- !(c("Att_Low", "Att_Med", "Att_High") %in% names(fixv))
  w <- tab$n
  tot <- sum(w)
  ll_old <- -Inf
  converged <- FALSE
  it <- 0L
  while (it < maxit) {
    it <- it + 1L
    lt <- adh_log_terms(pi, q, tab)
    lse <- adh_row_lse(lt)
    ll <- sum(w * lse)
    r <- exp(lt - lse)
    if (is.finite(ll_old) && abs(ll - ll_old) < tol) {
      converged <- TRUE
      break
    }
    ll_old <- ll
    wr <- w * r
    pi <- colSums(wr) / tot
    num <- colSums(wr * tab$A)
    den <- colSums(wr * (tab$A + tab$F))
    upd <- free_q & den > 0
    q[upd] <- num[upd] / den[upd]
  }
  par <- c(
    Pi_Low = pi[1],
    Pi_High = pi[3],
    Att_Low = q[1],
    Att_Med = q[2],
    Att_High = q[3]
  )
  list(
    par = par,
    loglik = adherence_loglik(par, tab),
    iterations = it,
    converged = converged
  )
}

# ----- Direct optimisation ----------------------------------------------------

#' Transform between the unconstrained vector z and the 5 parameters.
#'
#' Class probabilities: both free -> softmax(z1, 0, z2) over (Low, Med, High); one
#' free (other fixed at f) -> (1 - f) * plogis(z); both fixed -> nothing free.
#' Free q's: plogis(z).
#'
#' Ref: plan 2026-11-24, Sec. 4 S2 (logit/softmax scale, as in R/bootstrap.R).
#'
#' @param fixv  Named numeric vector of fixed parameters.
#' @return list(n, to_par = function(z), to_z = function(par)).
adh_transform <- function(fixv) {
  fixed <- names(fixv)
  free_pi <- setdiff(c("Pi_Low", "Pi_High"), fixed)
  free_q <- setdiff(c("Att_Low", "Att_Med", "Att_High"), fixed)
  clip <- function(p) pmin(pmax(p, 1e-8), 1 - 1e-8)
  to_par <- function(z) {
    par <- numeric(5)
    names(par) <- ADH_FREE_NAMES
    par[fixed] <- fixv
    k <- 0
    if (length(free_pi) == 2) {
      e <- exp(c(z[1], 0, z[2]) - max(z[1], 0, z[2]))
      p <- e / sum(e)
      par[c("Pi_Low", "Pi_High")] <- p[c(1, 3)]
      k <- 2
    } else if (length(free_pi) == 1) {
      other <- setdiff(c("Pi_Low", "Pi_High"), free_pi)
      par[free_pi] <- (1 - fixv[[other]]) * plogis(z[1])
      k <- 1
    }
    par[free_q] <- plogis(z[k + seq_along(free_q)])
    par
  }
  to_z <- function(par) {
    par <- par[ADH_FREE_NAMES]
    z <- numeric(0)
    if (length(free_pi) == 2) {
      pa <- clip(1 - par[["Pi_Low"]] - par[["Pi_High"]])
      z <- c(log(clip(par[["Pi_Low"]]) / pa), log(clip(par[["Pi_High"]]) / pa))
    } else if (length(free_pi) == 1) {
      other <- setdiff(c("Pi_Low", "Pi_High"), free_pi)
      room <- 1 - fixv[[other]]
      z <- qlogis(clip(par[[free_pi]] / max(room, 1e-8)))
    }
    c(z, qlogis(clip(par[free_q])))
  }
  list(n = length(free_pi) + length(free_q), to_par = to_par, to_z = to_z)
}

# Penalty for a non-finite negative log-likelihood inside optim (as FIT_PENALTY
# in R/bootstrap.R).
ADH_PENALTY <- 1e10

#' One optim run from one start: Nelder-Mead, then a BFGS polish; keep the
#' better.
#'
#' Ref: R/bootstrap.R fit_one_start() (same NM + polish pattern).
#'
#' @param start  Named 5-vector.
#' @param tab    Count table.
#' @param fixv   Named numeric vector of fixed parameters.
#' @param tol,maxit  Relative tolerance and iteration cap.
#' @return list(par, loglik, iterations, converged).
adh_optim_one <- function(start, tab, fixv, tol, maxit) {
  tr <- adh_transform(fixv)
  if (tr$n == 0) {
    par <- tr$to_par(numeric(0))
    return(list(
      par = par,
      loglik = adherence_loglik(par, tab),
      iterations = 0L,
      converged = TRUE
    ))
  }
  obj <- function(z) {
    v <- -adherence_loglik(tr$to_par(z), tab)
    if (is.finite(v)) v else ADH_PENALTY
  }
  z0 <- tr$to_z(start)
  nm <- if (tr$n == 1) {
    stats::optim(z0, obj, method = "BFGS", control = list(reltol = tol, maxit = maxit))
  } else {
    stats::optim(z0, obj, control = list(reltol = tol, maxit = maxit))
  }
  pol <- tryCatch(
    stats::optim(
      nm$par,
      obj,
      method = "BFGS",
      control = list(reltol = tol, maxit = 1000)
    ),
    error = function(e) NULL
  )
  best <- if (!is.null(pol) && pol$value <= nm$value) pol else nm
  par <- tr$to_par(best$par)
  list(
    par = par,
    loglik = adherence_loglik(par, tab),
    iterations = as.integer(sum(best$counts, na.rm = TRUE)),
    converged = best$convergence == 0
  )
}

# ----- Fit --------------------------------------------------------------------

#' Maximum-likelihood fit of the adherence model, multi-start.
#'
#' Runs EM or direct optimisation (Nelder-Mead + BFGS on the logit/softmax
#' scale) from deterministic starts, keeps the best log-likelihood and orders
#' the classes Att_Low <= Att_Med <= Att_High (see adh_order_classes(): classes with a
#' fixed q or fixed pi stay in their slot; the free ones are sorted among the
#' free slots).
#'
#' Ref: plan 2026-11-24, Sec. 4 S2; pattern of fit_mle() in R/bootstrap.R.
#'
#' @param tab       Count table with columns A, F, n.
#' @param method    "em" (default) or "optim".
#' @param fix       Named list of parameters pinned at a value, e.g.
#'   list(Att_Low = 0). EM supports fixed q's only; use "optim" to fix
#'   Pi_Low or Pi_High.
#' @param starts    NULL (use adherence_starts(n_starts)) or a list of named
#'   5-vectors.
#' @param n_starts  Number of default starts.
#' @param tol       EM: log-likelihood change; optim: relative tolerance.
#' @param maxit     Iteration cap per run.
#' @return List: `par` (named, ADH_PAR_NAMES order, six values), `loglik`, `n_par_free`,
#'   `convergence` (list: converged = all runs, best_converged, iterations of
#'   the best run), `method`, `fix`, `at_bound` (named logical) and `starts`
#'   (data frame: start, loglik, converged; for spotting several optima).
fit_adherence <- function(
  tab,
  method = c("em", "optim"),
  fix = list(),
  starts = NULL,
  n_starts = 20,
  tol = 1e-10,
  maxit = 1e5
) {
  method <- match.arg(method)
  fixv <- adh_check_fix(fix)
  if (method == "em" && any(c("Pi_Low", "Pi_High") %in% names(fixv))) {
    stop("method = \"em\" cannot fix Pi_Low or Pi_High; use method = \"optim\"")
  }
  if (is.null(starts)) {
    starts <- adherence_starts(n_starts)
  }
  runner <- if (method == "em") adh_em_one else adh_optim_one
  runs <- lapply(starts, function(s) {
    s <- unlist(s)[ADH_FREE_NAMES]
    runner(s, tab, fixv, tol, maxit)
  })
  lls <- vapply(runs, function(r) r$loglik, numeric(1))
  best <- runs[[which.max(lls)]]
  par5 <- adh_order_classes(best$par, names(fixv))
  par <- adh_full_par(par5)
  list(
    par = par,
    loglik = adherence_loglik(par5, tab),
    n_par_free = length(ADH_FREE_NAMES) - length(fixv),
    convergence = list(
      converged = all(vapply(runs, function(r) r$converged, logical(1))),
      best_converged = best$converged,
      iterations = best$iterations
    ),
    method = method,
    fix = fixv,
    at_bound = par <= 1e-6 | par >= 1 - 1e-6,
    starts = data.frame(
      start = seq_along(runs),
      loglik = lls,
      converged = vapply(runs, function(r) r$converged, logical(1))
    )
  )
}

# ----- Profile-likelihood CI --------------------------------------------------

#' Profile-likelihood confidence interval for one parameter.
#'
#' The profile log-likelihood at value v refits the other free parameters with
#' `which` pinned at v (method "optim", started from the MLE and a few
#' deterministic starts). The limits solve profile(v) = loglik_max -
#' qchisq(level, 1) / 2 by uniroot() on each side of the estimate. If the
#' profile is still above the cut-off at 0 (lower) or 1 (upper), that limit is
#' the boundary. Classes pinned by the profile keep their slot (no ordering
#' constraint across it), see adh_order_classes(); a multimodal profile can
#' make the interval too narrow, which this function does not detect.
#'
#' Ref: Venzon & Moolgavkar 1988, Appl Stat 37:87-94 (profile-likelihood
#' intervals by root finding).
#'
#' @param fit    A fit_adherence() result.
#' @param tab    Count table.
#' @param which  One name from ADH_FREE_NAMES.
#' @param level  Confidence level.
#' @param fix    Parameters held fixed throughout (as in fit_adherence()).
#' @return List: which, level, estimate, lower, upper, cutoff.
adherence_profile_ci <- function(fit, tab, which, level = 0.95, fix = list()) {
  stopifnot(length(which) == 1, which %in% ADH_FREE_NAMES)
  fixv <- adh_check_fix(fix)
  if (which %in% names(fixv)) {
    stop("`which` is fixed; nothing to profile")
  }
  est <- fit$par[[which]]
  cutoff <- fit$loglik - stats::qchisq(level, 1) / 2
  starts_extra <- adherence_starts(4)
  prof <- function(v) {
    fv <- c(fixv, setNames(v, which))
    st <- c(list(replace(fit$par, which, v)), starts_extra)
    adh_fit_ll(tab, fv, st)
  }
  g <- function(v) prof(v) - cutoff
  edge <- function(bound) {
    if (est == bound) {
      return(bound)
    }
    if (g(bound) >= 0) {
      return(bound)
    }
    stats::uniroot(g, sort(c(bound, est)), tol = 1e-8)$root
  }
  lower <- edge(0)
  upper <- edge(1)
  list(
    which = which,
    level = level,
    estimate = est,
    lower = lower,
    upper = upper,
    cutoff = cutoff
  )
}

#' Best log-likelihood over starts with the optim runner (profile helper).
#'
#' Ref: none (helper).
#'
#' @param tab    Count table.
#' @param fixv   Named numeric vector of fixed parameters.
#' @param starts List of named 5-vectors.
#' @return Scalar log-likelihood.
adh_fit_ll <- function(tab, fixv, starts) {
  lls <- vapply(
    starts,
    function(s) adh_optim_one(unlist(s), tab, fixv, 1e-12, 5000)$loglik,
    numeric(1)
  )
  max(lls)
}

# ----- Bootstrap --------------------------------------------------------------

#' Draw bootstrap resamples of the people (multinomial on the cells).
#'
#' Uses the caller's RNG state (set the seed in the analysis script before
#' calling). Each row is one resample of sum(n) people; column j is the count
#' for cell j of `tab`.
#'
#' Ref: Efron & Tibshirani 1993, An Introduction to the Bootstrap, Ch. 6
#' (nonparametric bootstrap of i.i.d. units).
#'
#' @param tab  Count table with columns A, F, n.
#' @param B    Number of resamples.
#' @return Integer matrix, B x nrow(tab).
draw_bootstrap_counts <- function(tab, B) {
  t(stats::rmultinom(B, size = sum(tab$n), prob = tab$n / sum(tab$n)))
}

#' Bootstrap the fit: refit to each row of pre-drawn resampled counts.
#'
#' Pure: no RNG inside (withr is a dev-only dependency, so the draws are an
#' argument; make them with draw_bootstrap_counts()).
#'
#' Ref: Efron & Tibshirani 1993, Ch. 6.
#'
#' @param tab        Count table (cells define the columns of `counts`).
#' @param B          Number of resamples to use (<= nrow(counts)).
#' @param fit_args   Named list of further arguments for fit_adherence().
#' @param counts     Matrix of resampled counts, one row per resample.
#' @return Matrix, B x 7: the six parameters (ADH_PAR_NAMES) and `loglik`.
adherence_bootstrap <- function(tab, B, fit_args = list(), counts) {
  stopifnot(ncol(counts) == nrow(tab), B <= nrow(counts))
  out <- matrix(NA_real_, B, 7, dimnames = list(NULL, c(ADH_PAR_NAMES, "loglik")))
  for (b in seq_len(B)) {
    tb <- tab
    tb$n <- counts[b, ]
    tb <- tb[tb$n > 0, ]
    f <- do.call(fit_adherence, c(list(tab = tb), fit_args))
    out[b, ] <- c(f$par, f$loglik)
  }
  out
}

# ----- Person-level simulator -------------------------------------------------

#' Simulate (A, F) for individual people under the model.
#'
#' Each person draws a class from pi, then walks months 0..T-1 (T = their
#' exposure): due at month 0; when due, attends with probability q_c each
#' month; after attending in month s they are not due for months s+1..s+m-1
#' and due again at s + m. This generates exactly the (A, F) that
#' due_exposure(scan_months, T, m) computes: F counts due months without a
#' scan, A the scans. Uses the caller's RNG state.
#'
#' Ref: plan 2026-11-24, Sec. 4 S2/S3; consistent with due_exposure() in
#' R/adherence_data.R.
#'
#' @param par             Named parameters (ADH_FREE_NAMES).
#' @param exposure_months Numeric vector, months of follow-up per person.
#' @param m               Months not due after an attended scan.
#' @return data.frame(A, F, class) with class in 1..3 (Low, Med, High).
simulate_adherence_people <- function(par, exposure_months, m = 6) {
  u <- adh_unpack(par)
  n <- length(exposure_months)
  cls <- sample.int(3, n, replace = TRUE, prob = u$pi)
  A <- integer(n)
  Fd <- integer(n)
  for (i in seq_len(n)) {
    q <- u$q[cls[i]]
    horizon <- exposure_months[i]
    d <- 0
    repeat {
      if (d >= horizon) {
        break
      }
      wait <- if (q <= 0) Inf else stats::rgeom(1, q)
      if (d + wait >= horizon) {
        Fd[i] <- Fd[i] + (horizon - d)
        break
      }
      A[i] <- A[i] + 1L
      Fd[i] <- Fd[i] + wait
      d <- d + wait + m
    }
  }
  data.frame(A = A, F = Fd, class = cls)
}
