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

# Box corner: the golden value was made with the FFT hypogeom. At the corner all
# five p_* are 1e-3, so stage 1 is a negative binomial and the exact-prefix
# hypogeom (CHANGES.md S14) is the accurate one: it matches dnbinom(33, 5, 1e-3)
# to 2e-15, while the old FFT was 5.4e-7 relative off (cdf 2.4e-6). neg_log_lik
# there moves by 5.6e-9 relative (loglik_medians by 7.5e-9), hence tolerance
# 1e-8 for this point only.
test_that("stage-1 terms at the box corner are exact", {
  # meds[1] = 38 months, minus 5 stages: x = 33
  expect_equal(calibration_data$stage$med[1] - 5, 33)
  expect_equal(
    hypogeom::dhypogeom(33, rep(1e-3, 5)),
    dnbinom(33, 5, 1e-3),
    tolerance = 1e-13
  )
  expect_equal(
    hypogeom::phypogeom(33, rep(1e-3, 5)),
    pnbinom(33, 5, 1e-3),
    tolerance = 1e-13
  )
})

test_that("neg_log_lik and loglik_medians match the golden points", {
  for (pt in golden_ll$points) {
    expect_equal(
      do.call(neg_log_lik, c(list(pt$par), pt$data)),
      pt$nll,
      tolerance = if (pt$label == "box corner") 1e-8 else 1e-12,
      label = paste0("neg_log_lik at ", pt$label)
    )
    expect_equal(
      loglik_medians(pt$par, pt$data$meds),
      pt$ll_med,
      tolerance = if (pt$label == "box corner") 1e-8 else 1e-12,
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
