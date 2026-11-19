# Invariant tests: properties that must hold for the CURRENT code.

test_that("7x7 transition matrices are row-stochastic", {
  # progress_untreated() and screen_cohort() both build their matrix with
  # transition_matrix(); the HCC probability is converted to monthly first.
  tm <- transition_matrix(
    annual_to_monthly(params_def$p_HCC),
    unlist(params_def[STAGE_P_NAMES])
  )
  expect_equal(dim(tm), c(7L, 7L), label = "transition_matrix dimensions")
  expect_equal(
    unname(rowSums(tm)),
    rep(1, 7),
    tolerance = 1e-12,
    label = "row sums of transition_matrix() (R/model_core.R)"
  )
  expect_true(
    all(tm >= 0 & tm <= 1),
    label = "entries of transition_matrix() lie in [0,1]"
  )
})

test_that("screen_cohort conserves mass", {
  res <- screen_cohort(n = 100, t = 120, m = 6, Pr_Att = 0.5)
  expect_true(
    all(res$dx >= 0) && all(res$undiagnosed >= 0),
    label = "screen_cohort: dx and undiagnosed are non-negative"
  )
  expect_equal(
    sum(res$dx) + sum(res$undiagnosed),
    100,
    tolerance = 1e-12,
    label = "screen_cohort: sum(Dx) + remaining == n at t = 120"
  )
  dx_long <- screen_cohort(n = 100, t = 2000, m = 6, Pr_Att = 0.5, p_HCC = 0.9)
  expect_equal(
    sum(dx_long$dx) + sum(dx_long$undiagnosed),
    100,
    tolerance = 1e-12,
    label = "screen_cohort: sum(Dx)+remaining == n (t=2000, p_HCC=0.9)"
  )
  expect_equal(
    sum(dx_long$undiagnosed),
    0,
    tolerance = 1e-8,
    label = "screen_cohort: remaining -> 0 in the long run (Mortality detected)"
  )
})

test_that("Pr_Att = 0 group is independent of m and Sens_*", {
  base <- screen_cohort(Pr_Att = 0)$dx
  expect_equal(
    screen_cohort(Pr_Att = 0, m = 12)$dx,
    base,
    tolerance = 1e-12,
    label = "screen_cohort(Pr_Att=0) independent of m (m=12 vs 6)"
  )
  expect_equal(
    screen_cohort(Pr_Att = 0, m = 1)$dx,
    base,
    tolerance = 1e-12,
    label = "screen_cohort(Pr_Att=0) independent of m (m=1 vs 6)"
  )
  expect_equal(
    screen_cohort(
      Pr_Att = 0,
      Sens_0 = 0.9,
      Sens_A = 0.9,
      Sens_B = 0.1,
      Sens_C = 0.1,
      Sens_D = 0.1
    )$dx,
    base,
    tolerance = 1e-12,
    label = "screen_cohort(Pr_Att=0) independent of all Sens_*"
  )
})

test_that("screen_by_adherence(p_NA=1, p_CA=0) = screen_cohort(Pr_Att=0)", {
  expect_equal(
    unlist(screen_by_adherence(p_NA = 1, p_CA = 0)),
    screen_cohort(n = 100, Pr_Att = 0)$dx,
    tolerance = 1e-12,
    label = "screen_by_adherence(NA=1,CA=0) == screen_cohort(Pr_Att=0)$dx"
  )
})

test_that("life_add is monotone increasing in Sens_0, p_CA and L_A", {
  vals <- list(
    Sens_0 = c(0.05, 0.1, 0.2, 0.4, 0.7),
    p_CA = c(0.1, 0.2, 0.3, 0.4, 0.5),
    L_A = c(1, 2, 3.6, 5, 8)
  )
  for (arg in names(vals)) {
    y <- vapply(
      vals[[arg]],
      function(v) do.call(life_add, setNames(list(v), arg)),
      numeric(1)
    )
    expect_true(
      all(diff(y) > 0),
      label = paste0(
        "life_add strictly increasing in ",
        arg,
        " over ",
        paste(vals[[arg]], collapse = ", "),
        " (got ",
        paste(signif(y, 6), collapse = ", "),
        ")"
      )
    )
  }
})

test_that("m > t gives no screening", {
  expect_equal(
    screen_cohort(Pr_Att = 1, t = 120, m = 121),
    screen_cohort(Pr_Att = 0, t = 120, m = 121),
    tolerance = 1e-12,
    label = "screen_cohort(Pr_Att=1, m=121 > t=120)$dx == no-screening"
  )
  expect_equal(
    life_add(m = 121),
    life_add(p_NA = 1, p_CA = 0),
    tolerance = 1e-12,
    label = "life_add(m > t) == life_add(p_NA=1,p_CA=0) (no screening)"
  )
})

test_that("15x15 likelihood matrix is row-stochastic", {
  tm <- nll_matrix()
  expect_equal(dim(tm), c(15L, 15L), label = "15x15 matrix dimensions")
  expect_equal(
    unname(rowSums(tm)),
    rep(1, 15),
    tolerance = 1e-12,
    label = "row sums of the 15x15 matrix inside neg_log_lik() (R/bootstrap.R)"
  )
})
