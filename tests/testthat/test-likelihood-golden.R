# Likelihood golden tests.
#
# RULE: fixtures/likelihood_golden.rds was generated ONCE with hypogeom
# damian-t-p@2b329c4 by fixtures/make_likelihood_golden.R. It pins
# neg_log_lik(), loglik_medians() and fit_mle(NM_START) across the hypogeom
# exact-prefix change (CHANGES.md). Never regenerate it to make a test pass.
# Doubles are compared at relative tolerance 1e-12 (repo convention: libm/BLAS
# noise is ~1e-15); the fit is looser (par 1e-7, value 1e-10) because the
# optimiser amplifies round-off.

golden_ll <- readRDS(here::here(
  "tests",
  "testthat",
  "fixtures",
  "likelihood_golden.rds"
))

test_that("neg_log_lik and loglik_medians match the golden points", {
  for (pt in golden_ll$points) {
    expect_equal(
      do.call(neg_log_lik, c(list(pt$par), pt$data)),
      pt$nll,
      tolerance = 1e-12,
      label = paste0("neg_log_lik at ", pt$label)
    )
    expect_equal(
      loglik_medians(pt$par, pt$data$meds),
      pt$ll_med,
      tolerance = 1e-12,
      label = paste0("loglik_medians at ", pt$label)
    )
  }
})

test_that("fit_mle(NM_START) matches the golden fit", {
  f <- fit_mle(list(NM_START))
  expect_equal(
    f$par,
    golden_ll$fit$par,
    tolerance = 1e-7,
    label = "fit_mle(NM_START)$par"
  )
  expect_equal(
    f$value,
    golden_ll$fit$value,
    tolerance = 1e-10,
    label = "fit_mle(NM_START)$value"
  )
})
