# Paper-structure tests for the likelihood in R/bootstrap.R
# (Cumming et al. 2024, Cancers 16:2745, and its Supplement).

# Independent construction of the 15x15 matrix from Supplementary Table S2 /
# matrix S1. States: 1 = dx BCLC 0, 2 = dx BCLC A, 3-8 = undiagnosed 0 with
# 6..1 months to US, 9-14 = undiagnosed A with 6..1 months to US,
# 15 = absorbing "progressed past A".
supp_matrix <- function(p0, pA, Sens0, SensA) {
  T <- matrix(0, 15, 15)
  T[1, 1] <- 1
  T[2, 2] <- 1
  T[15, 15] <- 1
  for (i in 3:7) {
    T[i, i + 1] <- 1 - p0
    T[i, i + 7] <- p0
  }
  T[8, 1] <- (1 - p0) * Sens0
  T[8, 2] <- p0 * SensA
  T[8, 3] <- (1 - p0) * (1 - Sens0)
  T[8, 9] <- p0 * (1 - SensA)
  for (i in 9:13) {
    T[i, i + 1] <- 1 - pA
    T[i, 15] <- pA
  }
  T[14, 2] <- (1 - pA) * SensA
  T[14, 9] <- (1 - pA) * (1 - SensA)
  T[14, 15] <- pA
  T
}

test_that("params_def$Sym_B == 8 * params_def$Sym_A (Supp. 2.1.3 footnote 4)", {
  expect_equal(
    params_def$Sym_B,
    8 * params_def$Sym_A,
    tolerance = 1e-6,
    label = "params_def$Sym_B == 8*params_def$Sym_A (Supp. Sec 2.1.3 footnote 4; data/Parameters.csv)"
  )
})

test_that("order-statistic arguments follow Supp. p.3 footnote 2", {
  # These constants are written inside neg_log_lik() (R/bootstrap.R) and are
  # not exposed, so this test checks the RULE against the constants copied
  # here. A future refactor of R/bootstrap.R should expose them (e.g. as a
  # data table) so this test can read them from the code directly.
  tab <- data.frame(
    stage = c("D", "C", "B", "A", "0"),
    n = c(210, 138, 139, 101, 12),
    code_mlen = c(209, 137, 139, 101, 11),
    code_j = c(210, 138, 140, 102, 12) / 2
  )
  # rule: mlen = n if n odd, else n - 1; j = (mlen + 1) / 2
  rule_mlen <- ifelse(tab$n %% 2 == 1, tab$n, tab$n - 1)
  rule_j <- (rule_mlen + 1) / 2
  expect_equal(
    tab$code_mlen,
    rule_mlen,
    label = "mlen = n if n odd else n-1, stages D,C,B,A,0 (Supp. p.3 footnote 2)"
  )
  expect_equal(
    tab$code_j,
    rule_j,
    label = "j = (mlen+1)/2, stages D,C,B,A,0 (Supp. p.3 footnote 2)"
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
  for (probs in list(hypo_probs, c(0.2, 0.1, 0.3, 0.4, 0.15, 0.3, 0.6, 0.02))) {
    expect_equal(
      nll_matrix(probs),
      supp_matrix(probs[1], probs[2], probs[6], probs[7]),
      tolerance = 0,
      ignore_attr = TRUE,
      label = "matrix_of_trans in neg_log_lik() == Supp. Table S2 / matrix S1 (exact)"
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
    label = "t(0) %*% T^100 [1] == closed-form Pr(d0), abs tol 1e-5 (Supp. footnote 3)"
  )
})
