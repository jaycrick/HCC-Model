# Invariant tests: properties that must hold for the CURRENT code.

p_def <- list(
  p_HCC = p_HCC_def,
  p_0 = p_0_def,
  p_A = p_A_def,
  p_B = p_B_def,
  p_C = p_C_def,
  p_D = p_D_def
)

test_that("7x7 transition matrices are row-stochastic", {
  for (nm in c("fromf4", "us_theor_pa")) {
    tm <- if (nm == "fromf4") fromf4_tmat() else us_theor_pa_tmat()
    expect_equal(dim(tm), c(7L, 7L), label = paste0(nm, " matrix dimensions"))
    expect_equal(
      rowSums(tm),
      rep(1, 7),
      tolerance = 1e-12,
      ignore_attr = TRUE,
      label = paste0("row sums of 7x7 matrix in ", nm, "() (R/model_core.R)")
    )
    expect_true(
      all(tm >= 0 & tm <= 1),
      label = paste0("entries of 7x7 matrix in ", nm, "() lie in [0,1]")
    )
  }
})

test_that("us_theor_pa conserves mass", {
  dx <- us_theor_pa(n = 100, t = 120, m = 6, Pr_Att = 0.5)
  expect_true(
    all(dx >= 0) && sum(dx) <= 100 + 1e-9,
    label = "us_theor_pa: 0 <= sum(Dx) <= n, so remaining = n - sum(Dx) >= 0"
  )
  # In the long run with a large HCC hazard everyone reaches Mortality, where
  # detection is certain, so remaining -> 0 and sum(Dx) == n.
  dx_long <- us_theor_pa(n = 100, t = 2000, m = 6, Pr_Att = 0.5, p_HCC = 0.9)
  expect_equal(
    sum(dx_long),
    100,
    tolerance = 1e-8,
    label = "us_theor_pa: sum(Dx)+remaining == n (remaining -> 0, t=2000, p_HCC=0.9)"
  )
})

test_that("Pr_Att = 0 group is independent of m and Sens_*", {
  base <- us_theor_pa(Pr_Att = 0)
  expect_equal(
    us_theor_pa(Pr_Att = 0, m = 12),
    base,
    tolerance = 1e-12,
    label = "us_theor_pa(Pr_Att=0) independent of m (m=12 vs 6)"
  )
  expect_equal(
    us_theor_pa(Pr_Att = 0, m = 1),
    base,
    tolerance = 1e-12,
    label = "us_theor_pa(Pr_Att=0) independent of m (m=1 vs 6)"
  )
  expect_equal(
    us_theor_pa(
      Pr_Att = 0,
      Sens_0 = 0.9,
      Sens_A = 0.9,
      Sens_B = 0.1,
      Sens_C = 0.1,
      Sens_D = 0.1
    ),
    base,
    tolerance = 1e-12,
    label = "us_theor_pa(Pr_Att=0) independent of all Sens_*"
  )
})

test_that("adh_theor(p_NA=1, p_CA=0) equals us_theor_pa(Pr_Att=0)", {
  expect_equal(
    as.numeric(unlist(adh_theor(p_NA = 1, p_CA = 0))),
    as.numeric(us_theor_pa(n = 100, Pr_Att = 0)),
    tolerance = 1e-12,
    label = "adh_theor(p_NA=1,p_CA=0) == us_theor_pa(n=100,Pr_Att=0)"
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
    us_theor_pa(Pr_Att = 1, t = 120, m = 121),
    us_theor_pa(Pr_Att = 0, t = 120, m = 121),
    tolerance = 1e-12,
    label = "us_theor_pa(Pr_Att=1, m=121 > t=120) == no-screening"
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
    rowSums(tm),
    rep(1, 15),
    tolerance = 1e-12,
    label = "row sums of the 15x15 matrix inside neg_log_lik() (R/bootstrap.R)"
  )
})
