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
  expect_equal(
    approx_d0,
    closed,
    tolerance = 1e-5,
    scale = 1,
    label = paste(
      "t(0) %*% T^100 [1] == closed-form Pr(d0), abs tol 1e-5",
      "(Supp. footnote 3)"
    )
  )
})

test_that("end_state_probs: T^100 vs closed form and vs T^1000 (fn 3)", {
  p0 <- csv_par[["p_0"]]
  s0 <- csv_par[["Sens_0"]]
  tm <- nll_matrix()
  d0_100 <- end_state_probs(tm, s0, t = 100)[["dx_0"]]
  d0_1000 <- end_state_probs(tm, s0, t = 1000)[["dx_0"]]
  closed <- s0 * (1 - (1 - p0)^6) / (6 * p0 * (1 - (1 - s0) * (1 - p0)^6))
  expect_lt(
    abs(d0_100 - closed),
    1e-5,
    label = "|Pr(d0) at t=100 - closed form| < 1e-5 (code matrix)"
  )
  expect_lt(
    abs(d0_100 - d0_1000),
    1e-5,
    label = "|Pr(d0) at t=100 - at t=1000| < 1e-5"
  )
  expect_equal(
    sum(end_state_probs(tm, s0)),
    1,
    tolerance = 1e-12,
    label = "end-state probabilities sum to 1"
  )
})

test_that("neg_log_lik is NaN outside [0,1]", {
  bad <- csv_par
  bad[1] <- 1.1
  expect_true(
    is.nan(neg_log_lik(bad)),
    label = "neg_log_lik(p_0 = 1.1) is NaN"
  )
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

test_that("G7: Nelder-Mead from the documented start reproduces the CSV", {
  start <- NM_START
  fit <- fit_point_estimate(start)
  for (nm in NLL_PAR_NAMES) {
    expect_lt(
      abs(fit$par[[nm]] - csv_par[[nm]]),
      1e-6,
      label = paste0("|NM estimate - CSV POINT| < 1e-6 for ", nm)
    )
  }
  expect_lt(
    abs(fit$value - neg_log_lik(csv_par)),
    1e-8,
    label = "|NM minimum - NLL(CSV)| < 1e-8"
  )
  expect_lt(
    neg_log_lik(csv_par),
    neg_log_lik(old_hypo_probs),
    label = "NLL(CSV) < NLL(stale pre-G7 hypo_probs)"
  )
})
