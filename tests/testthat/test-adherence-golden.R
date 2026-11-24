# Adherence golden tests.
#
# RULE: fixtures/adherence_golden.rds was generated ONCE by
# fixtures/make_adherence_golden.R. It pins fit_adherence() on the committed
# table papers/stochastic/data/adherence_counts.csv (free 3-class fit and the
# Att_Low = 0 fit). Never regenerate it to make a test pass. Tolerances follow
# test-likelihood-golden.R: the log-likelihood is compared at 1e-10 and the
# parameters at 1e-6, because the optimiser amplifies round-off (and a
# parameter at a bound is only determined to the convergence tolerance).

golden_adh <- readRDS(here::here(
  "tests",
  "testthat",
  "fixtures",
  "adherence_golden.rds"
))
tab_golden <- read.csv(here::here(
  "papers",
  "stochastic",
  "data",
  "adherence_counts.csv"
))

test_that("free 3-class fit matches the golden fit", {
  f <- fit_adherence(tab_golden)
  expect_equal(f$par, golden_adh$free$par, tolerance = 1e-6, label = "free par")
  expect_equal(
    f$loglik,
    golden_adh$free$loglik,
    tolerance = 1e-10,
    label = "free loglik"
  )
})

test_that("Att_Low = 0 fit matches the golden fit", {
  f <- fit_adherence(tab_golden, fix = list(Att_Low = 0))
  expect_equal(
    f$par,
    golden_adh$att_low_0$par,
    tolerance = 1e-6,
    label = "Att_Low=0 par"
  )
  expect_equal(
    f$loglik,
    golden_adh$att_low_0$loglik,
    tolerance = 1e-10,
    label = "Att_Low=0 loglik"
  )
})
