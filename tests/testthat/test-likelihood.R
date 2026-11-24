# Paper-structure tests for the likelihood in R/bootstrap.R
# (Cumming et al. 2024, Cancers 16:2745, and its Supplement).

# Independent construction of the 15x15 matrix from Supplementary Table S2 /
# matrix S1. States: 1 = dx BCLC 0, 2 = dx BCLC A, 3-8 = undiagnosed 0 with
# 6..1 months to US, 9-14 = undiagnosed A with 6..1 months to US,
# 15 = absorbing "progressed past A".
supp_matrix <- function(p0, pA, Sens0, SensA) {
  trans <- matrix(0, 15, 15)
  trans[1, 1] <- 1
  trans[2, 2] <- 1
  trans[15, 15] <- 1
  for (i in 3:7) {
    trans[i, i + 1] <- 1 - p0
    trans[i, i + 7] <- p0
  }
  trans[8, 1] <- (1 - p0) * Sens0
  trans[8, 2] <- p0 * SensA
  trans[8, 3] <- (1 - p0) * (1 - Sens0)
  trans[8, 9] <- p0 * (1 - SensA)
  for (i in 9:13) {
    trans[i, i + 1] <- 1 - pA
    trans[i, 15] <- pA
  }
  trans[14, 2] <- (1 - pA) * SensA
  trans[14, 9] <- (1 - pA) * (1 - SensA)
  trans[14, 15] <- pA
  trans
}

test_that("params_def$Sym_B == 8 * params_def$Sym_A (Supp. 2.1.3 footnote 4)", {
  expect_equal(
    params_def$Sym_B,
    8 * params_def$Sym_A,
    tolerance = 1e-6,
    label = paste(
      "params_def$Sym_B == 8*params_def$Sym_A",
      "(Supp. Sec 2.1.3 footnote 4; data/Parameters.csv)"
    )
  )
})

test_that("order-statistic arguments follow Supp. p.3 footnote 2", {
  cal <- calibration_data$stage
  expect_equal(
    cal$n,
    c(12, 101, 139, 138, 210),
    label = "calibration_data n for stages 0,A,B,C,D (Supp. p.3)"
  )
  expect_equal(
    cal$mlen,
    c(11, 101, 139, 137, 209),
    label = "calibration_data mlen: n if odd else n-1 (Supp. p.3 footnote 2)"
  )
  expect_equal(
    cal$j,
    c(6, 51, 70, 69, 105),
    label = "calibration_data j = (mlen+1)/2 (Supp. p.3 footnote 2)"
  )
})

test_that("likelihood_matrix has labelled states", {
  tm <- nll_matrix()
  expect_identical(
    rownames(tm),
    colnames(tm),
    label = "row and column state names agree"
  )
  expect_identical(
    rownames(tm)[c(1, 2, 3, 8, 9, 14, 15)],
    c("dx_0", "dx_A", "u0_6", "u0_1", "uA_6", "uA_1", "past_A"),
    label = "state labels at positions 1,2,3,8,9,14,15"
  )
})

test_that("independent matrix from Supp. Table S2 is row-stochastic", {
  pp <- read_params(csv_main)[["POINT"]]
  names(pp) <- rownames(read_params(csv_main))
  tm <- supp_matrix(pp[["p_0"]], pp[["p_A"]], pp[["Sens_0"]], pp[["Sens_A"]])
  expect_equal(
    rowSums(tm),
    rep(1, 15),
    tolerance = 1e-12,
    label = "row sums of independent 15x15 matrix (Supp. Table S2 / matrix S1)"
  )
})

test_that("code's 15x15 matrix equals the independent Supp. S2 matrix", {
  for (probs in list(csv_par, c(0.2, 0.1, 0.3, 0.4, 0.15, 0.3, 0.6, 0.02))) {
    expect_equal(
      unname(nll_matrix(probs)),
      supp_matrix(probs[1], probs[2], probs[6], probs[7]),
      tolerance = 0,
      ignore_attr = TRUE,
      label = paste(
        "matrix_of_trans in neg_log_lik() == Supp. Table S2 / matrix S1",
        "(exact)"
      )
    )
  }
})

test_that("t(0) %*% T^100 converges to the closed form Pr(d0) (Supp. fn 3)", {
  params <- read_params(csv_main)
  p0 <- params["p_0", "POINT"]
  pA <- params["p_A", "POINT"]
  s0 <- params["Sens_0", "POINT"]
  sA <- params["Sens_A", "POINT"]
  tm <- supp_matrix(p0, pA, s0, sA)
  t0 <- c(s0 / 6, 0, (1 - s0) / 6, rep(1 / 6, 5), rep(0, 7))
  expect_equal(sum(t0), 1, tolerance = 1e-12, label = "t(0) sums to 1")
  approx_d0 <- (t0 %*% (tm %^% 100))[1]
  closed <- s0 * (1 - (1 - p0)^6) / (6 * p0 * (1 - (1 - s0) * (1 - p0)^6))
  expect_lte(
    abs(approx_d0 - closed),
    1e-5,
    label = paste(
      "t(0) %*% T^100 [1] == closed-form Pr(d0), abs tol 1e-5",
      "(Supp. footnote 3)"
    )
  )
})

test_that("end_state_probs: exact vs closed form, T^100, T^1000 (fn 3, G8)", {
  p0 <- csv_par[["p_0"]]
  s0 <- csv_par[["Sens_0"]]
  tm <- nll_matrix()
  exact <- end_state_probs(tm, s0)
  d0_100 <- end_state_probs(tm, s0, t = 100)[["dx_0"]]
  d0_1000 <- end_state_probs(tm, s0, t = 1000)[["dx_0"]]
  closed <- s0 * (1 - (1 - p0)^6) / (6 * p0 * (1 - (1 - s0) * (1 - p0)^6))
  expect_lte(
    abs(exact[["dx_0"]] - closed),
    1e-12,
    label = "exact Pr(d0) == closed form"
  )
  expect_lte(abs(exact[["dx_0"]] - d0_1000), 1e-9, label = "exact Pr(d0) == T^1000")
  expect_equal(
    sum(exact),
    1,
    tolerance = 1e-12,
    label = "exact end-state probabilities sum to 1"
  )
  # The old approximation (documented): error 1.445e-6 at t = 100.
  expect_equal(
    abs(d0_100 - closed),
    1.445e-6,
    tolerance = 1e-3,
    label = "|T^100 - closed form| is the documented 1.445e-6"
  )
  expect_lt(abs(d0_100 - d0_1000), 1e-5, label = "T^100 vs T^1000 < 1e-5")
})

test_that("neg_log_lik is NaN outside [0,1]", {
  bad <- csv_par
  bad[1] <- 1.1
  expect_true(is.nan(neg_log_lik(bad)), label = "neg_log_lik(p_0 = 1.1) is NaN")
})

# ----- G7: the CSV point estimates are the MLE --------------------------------
# Reference only: the stale values hard-coded in analysis/01 before G7. NLL

old_hypo_probs <- c(
  0.07982013,
  0.05540597,
  0.33922926,
  0.35681215,
  0.16023503,
  0.21878495,
  0.54017653,
  0.01015650
)

test_that("fit_mle respects the box and reaches the CSV MLE (G9)", {
  fit <- fit_mle(list(NM_START))
  nm_only <- stats::optim(
    NM_START,
    neg_log_lik,
    method = "Nelder-Mead",
    control = list(maxit = 10000)
  )
  expect_identical(names(fit$par), NLL_PAR_NAMES)
  expect_true(all(fit$par >= FIT_LB & fit$par <= FIT_UB), label = "in box")
  expect_lte(fit$value, nm_only$value, label = "value <= NM-only value")
  expect_lte(
    fit$value,
    neg_log_lik(csv_par) + 1e-8,
    label = "value <= NLL(CSV) + 1e-8 (the CSV is the MLE)"
  )
  expect_lt(fit$par[["Sym_A"]], 1 / 8, label = "Sym_A < 1/8")
  expect_identical(fit$convergence, 0L)
})

test_that("fit_at_bound flags values on or outside the box", {
  par <- csv_par
  expect_false(any(fit_at_bound(par)))
  par[["p_B"]] <- 1
  par[["p_C"]] <- FIT_UB[["p_C"]] * plogis(qlogis(1 - 1e-6))
  par[["p_0"]] <- FIT_LB[["p_0"]]
  par[["Sym_A"]] <- FIT_UB[["Sym_A"]] * plogis(qlogis(1 - 1e-6))
  flagged <- names(which(fit_at_bound(par)))
  expect_setequal(flagged, c("p_B", "p_C", "p_0", "Sym_A"))
})

test_that("polish-only starts are clipped and polished", {
  start <- list(par = old_hypo_probs, nm = FALSE)
  fit <- fit_mle(list(start))
  expect_true(is.na(fit$starts$nm_value))
  expect_lt(fit$value, neg_log_lik(old_hypo_probs))
  expect_lt(fit$par[["Sym_A"]], 1 / 8)
})

test_that("refit_bootstrap does not depend on cores (2 resamples)", {
  inputs <- readRDS(paper_path("data-derived", "list_of_resampled_milton_output.rds"))[
    1:2
  ]
  hpc <- bootstrap_estimates()
  starts_fn <- function(b) list(list(par = unlist(hpc[b, ]), nm = FALSE))
  one <- refit_bootstrap(inputs, starts_fn, cores = 1)
  two <- refit_bootstrap(inputs, starts_fn, cores = 2)
  expect_identical(one, two)
  expect_identical(names(one$pars), NLL_PAR_NAMES)
  expect_identical(nrow(one$diag), 2L)
})

test_that("PUBLISHED_MLE is the published MLE in NLL_PAR_NAMES order", {
  expect_identical(names(PUBLISHED_MLE), NLL_PAR_NAMES)
  # Guards against typos: the values rounded to 2 significant digits
  expect_equal(
    unname(signif(PUBLISHED_MLE, 2)),
    c(0.08, 0.055, 0.34, 0.35, 0.16, 0.21, 0.53, 0.01)
  )
})

test_that("the refit's point fit from the fixed starts reproduces the CSV", {
  fit <- fit_mle(list(NM_START, PUBLISHED_MLE))
  expect_lt(
    max(abs(fit$par - csv_par)),
    1e-7,
    label = "max |point fit - CSV POINT| (CSV written at 8 significant digits)"
  )
})

# Expected to FAIL until the CSV is regenerated from fit_mle (plan step T3).
test_that("G9: fit_mle from c(NM_START, CSV) reproduces the CSV", {
  fit <- fit_mle(list(NM_START, csv_par))
  for (nm in NLL_PAR_NAMES) {
    expect_lt(
      abs(fit$par[[nm]] - csv_par[[nm]]),
      1e-6,
      label = paste0("|fit_mle estimate - CSV POINT| < 1e-6 for ", nm)
    )
  }
  expect_lt(
    abs(fit$value - neg_log_lik(csv_par)),
    1e-8,
    label = "|fit_mle minimum - NLL(CSV)| < 1e-8"
  )
  expect_lt(
    neg_log_lik(csv_par),
    neg_log_lik(old_hypo_probs),
    label = "NLL(CSV) < NLL(stale pre-G7 hypo_probs)"
  )
})
