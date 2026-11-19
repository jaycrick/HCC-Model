# Golden (characterisation) tests.
#
# RULE: fixtures/golden.rds was generated ONCE from the unmodified code by
# fixtures/make_golden.R. These tests pin the code AS IT IS, bugs included.
# The fixture is never regenerated except under an approved gate (G1..G7),
# and the gate is recorded in the commit message. Doubles are compared at
# relative tolerance 1e-12. This is portable across platforms (libm/BLAS noise
# is ~1e-15), and real model changes are >= 1e-3.

golden <- readRDS(golden_path())
rows20 <- resampled_pars_fin[1:20, ]

test_that("life_add() anchors are bit-for-bit stable", {
  expect_equal(
    life_add(),
    golden$life_add_baseline,
    tolerance = 1e-12,
    label = "life_add() baseline anchor (paper Supp Fig S3; README)"
  )
  expect_equal(
    life_add(p_NA = 1, p_CA = 0),
    golden$life_add_no_screening,
    tolerance = 1e-12,
    label = "life_add(p_NA=1, p_CA=0) no-screening anchor (README)"
  )
  # plan anchors, independent of the fixture
  expect_equal(
    life_add(),
    68.023264456,
    tolerance = 1e-9,
    label = "life_add() vs plan anchor 68.023264456"
  )
  expect_equal(
    life_add(p_NA = 1, p_CA = 0),
    14.8567745055,
    tolerance = 1e-9,
    label = "life_add(p_NA=1,p_CA=0) vs plan anchor 14.8567745055"
  )
})

test_that("screen_by_adherence() 7-vector is stable", {
  expect_equal(
    unlist(screen_by_adherence()),
    golden$screen_by_adherence_default,
    tolerance = 1e-12,
    label = "screen_by_adherence() default stage vector (R/model_core.R)"
  )
})

test_that("progress_untreated(n=100,t=120) is stable (post gate G1)", {
  # Gate G1 (finding C1) removed a second annual->monthly conversion of the
  # already-monthly p_0..p_D: Mortality 0.09 -> 13.32. Fixture entry
  # regenerated under G1; see CHANGES.md.
  expect_equal(
    unlist(progress_untreated(n = 100, t = 120)),
    golden$progress_untreated_100_120,
    tolerance = 1e-12,
    label = "progress_untreated(n=100,t=120) (fixture regenerated, G1)"
  )
})

test_that("neg_log_lik(hypo_probs) value is stable", {
  expect_equal(
    neg_log_lik(hypo_probs, meds = meds),
    golden$neg_log_lik_hypo_probs,
    tolerance = 1e-12,
    label = "neg_log_lik(hypo_probs) NLL (R/bootstrap.R)"
  )
})

# life-years of one scenario on rows20 at base risk, via the scenarios API
scen_ly <- function(id) {
  # nolint next: object_usage_linter.
  long <- run_scenarios(rows20, risk = risk_multipliers["base"])
  long$ly[long$scenario == id]
}

test_that("scenarios reproduce the golden fixture (tolerance 1e-12)", {
  ids <- c(
    base = "baseline",
    `1a` = "real_adh",
    `1b` = "opt_adh",
    `2` = "diag_sens",
    `3` = "treat"
  )
  for (nm in names(ids)) {
    expect_equal(
      scen_ly(ids[[nm]]),
      golden[[paste0("scen_", nm)]],
      tolerance = 1e-12,
      label = paste0(
        "scen_",
        nm,
        " (",
        ids[[nm]],
        ") on first 20 rows of resampled_pars_fin (golden fixture)"
      )
    )
  }
})

test_that("scenarios match shipped resampled_pars_fin columns", {
  # NOTE: the shipped columns were computed on another machine/run and differ
  # from a fresh computation by ~2e-16 relative (floating-point noise), so a
  # tolerance of exactly 0 is false here; 1e-12 is the documented
  # float-reordering tolerance. The golden test above compares against a locally generated fixture.
  stored <- scenario_table$column
  names(stored) <- scenario_table$id
  for (id in names(stored)) {
    expect_equal(
      scen_ly(id),
      rows20[[stored[[id]]]],
      tolerance = 1e-12,
      label = paste0(
        id,
        " vs stored column '",
        stored[[id]],
        "' in data-derived/resampled_pars_fin.rds (rows 1-20)"
      )
    )
  }
})

test_that("life_tibble() is stable for both parameter CSVs", {
  expect_equal(
    life_tibble(csv_main),
    golden$life_tibble_main,
    tolerance = 1e-12,
    label = "life_tibble(data/Parameters.csv) (R/model_core.R)"
  )
  expect_equal(
    life_tibble(csv_25),
    golden$life_tibble_25,
    tolerance = 1e-12,
    label = "life_tibble(data/Parameters_25.csv) (R/model_core.R)"
  )
})

test_that("param_names.rds names are stable", {
  # NB: analysis scripts mutate this tracked file; the fixture pins the
  # names as they were when S1 was run (working-tree copy).
  expect_identical(
    names(readRDS(here::here("data-derived", "param_names.rds"))),
    golden$param_names,
    label = "names(data-derived/param_names.rds)"
  )
})
