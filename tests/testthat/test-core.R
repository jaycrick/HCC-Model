# Tests for the S3 core refactor: parameter handling, the single transition
# matrix builder, annual->monthly conversion, and the progress_untreated /
# screen_cohort engine agreement (gate G1).

test_that("life_add() rejects unknown parameter names (finding C2)", {
  expect_error(
    life_add(p_att = 1),
    "unknown parameter(s): p_att",
    fixed = TRUE,
    label = "life_add(p_att = 1) errors (typo for p_Att not swallowed)"
  )
  expect_error(
    life_add(p_att = 1, foo = 2),
    "unknown parameter(s): p_att, foo",
    fixed = TRUE,
    label = "every unknown name is listed in the error message"
  )
  expect_error(
    screen_by_adherence(p_att = 1),
    "unknown parameter(s): p_att",
    fixed = TRUE,
    label = "screen_by_adherence() rejects unknown parameter names"
  )
  expect_error(
    modify_params(params_def, list(p_att = 1)),
    "unknown parameter(s): p_att",
    fixed = TRUE,
    label = "modify_params() rejects unknown parameter names"
  )
})

test_that("modify_params() overrides by name and leaves others alone", {
  p <- modify_params(params_def, list(p_Att = 0.5))
  expect_equal(p$p_Att, 0.5, label = "overridden value")
  expect_equal(
    p[setdiff(names(p), "p_Att")],
    params_def[setdiff(names(params_def), "p_Att")],
    label = "all other parameters unchanged"
  )
  expect_identical(
    modify_params(params_def, list()),
    params_def,
    label = "empty override list returns params unchanged"
  )
})

test_that("load_params() assigns by PARAMETER name and handles the BOM", {
  expect_equal(length(params_def), 24L, label = "24 parameters in the CSV")
  csv <- read.csv(csv_main, fileEncoding = "UTF-8-BOM")
  expect_identical(
    names(params_def),
    csv$PARAMETER,
    label = "names(params_def) == CSV PARAMETER column (BOM stripped)"
  )
  expect_equal(
    unlist(params_def),
    setNames(csv$POINT, csv$PARAMETER),
    label = "values are the POINT column, matched by name"
  )
  expect_equal(
    load_params(csv_main, column = "UPPER")[["p_HCC"]],
    0.04,
    label = "column = 'UPPER' reads the UPPER bound by name"
  )
})

test_that("annual_to_monthly() inverts the monthly->annual compounding", {
  # Fleurence & Hollenbeak 2007, PharmacoEconomics 25(1):3-6
  expect_equal(
    annual_to_monthly(1 - (1 - 0.1)^12),
    0.1,
    label = "annual_to_monthly(1-(1-0.1)^12) == 0.1"
  )
  expect_equal(annual_to_monthly(0), 0, label = "zero maps to zero")
  expect_equal(annual_to_monthly(1), 1, label = "one maps to one")
})

test_that("transition_matrix() is a named 7x7 row-stochastic matrix", {
  tm <- transition_matrix(
    annual_to_monthly(params_def$p_HCC),
    unlist(params_def[c("p_0", "p_A", "p_B", "p_C", "p_D")])
  )
  expect_equal(dim(tm), c(7L, 7L), label = "7x7")
  expect_identical(rownames(tm), stage_list, label = "rows named by stage_list")
  expect_identical(colnames(tm), stage_list, label = "cols named by stage_list")
  expect_equal(
    unname(rowSums(tm)),
    rep(1, 7),
    tolerance = 1e-12,
    label = "row sums are 1"
  )
  expect_true(all(tm >= 0 & tm <= 1), label = "entries in [0,1]")
  expect_equal(tm["Mortality", "Mortality"], 1, label = "Mortality absorbing")
  expect_equal(
    tm["BCLC A", "BCLC B"],
    params_def$p_A,
    label = "p_A is the BCLC A -> BCLC B probability"
  )
})

test_that("screen_cohort() exposes the undiagnosed state (mass conservation)", {
  for (pa in c(0, 0.449, 1)) {
    res <- screen_cohort(Pr_Att = pa, n = 100, t = 120, m = 6)
    expect_identical(
      names(res$dx),
      stage_list,
      label = "dx named by stage_list"
    )
    expect_identical(
      names(res$undiagnosed),
      stage_list,
      label = "undiagnosed named by stage_list"
    )
    expect_equal(
      sum(res$dx) + sum(res$undiagnosed),
      100,
      tolerance = 1e-12,
      label = paste0(
        "sum(Dx) + remaining == n at t = 120 (Pr_Att = ",
        pa,
        ")"
      )
    )
  }
})

test_that("progress_untreated() matches screen_cohort() w/o detection (G1)", {
  # With Pr_Att = 0 and every Sym_* = 0 nobody is ever diagnosed except by
  # death, so the undiagnosed states plus the deaths counted in dx must equal
  # the plain Markov progression. This FAILS while progress_untreated applies
  # a second annual->monthly conversion to the monthly p_0..p_D (C1).
  no_sym <- as.list(setNames(
    rep(0, 5),
    paste0("Sym_", c("0", "A", "B", "C", "D"))
  ))
  for (tt in c(12, 120)) {
    sc <- do.call(
      screen_cohort,
      c(list(Pr_Att = 0, n = 100, t = tt, m = 6), no_sym)
    )
    expect_equal(
      unlist(progress_untreated(n = 100, t = tt)),
      sc$undiagnosed + sc$dx,
      tolerance = 1e-10,
      label = paste0(
        "progress_untreated(100, ",
        tt,
        ") == undiagnosed + dx of screen_cohort (Pr_Att=0, Sym_*=0)"
      )
    )
  }
})
