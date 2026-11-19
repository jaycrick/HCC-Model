# Scenarios as data (R/scenarios.R): reproduction of the shipped bootstrap
# columns, point estimates against plain life_add() calls, and the
# "one place" property of the scenario definitions.

rows20 <- resampled_pars_fin[1:20, ]

test_that("run_scenarios reproduces all stored scenario columns", {
  wide <- scenario_columns(run_scenarios(rows20))
  legacy <- names(resampled_pars_fin)
  legacy <- legacy[!legacy %in% names(draws(resampled_pars_fin))]
  expect_identical(names(wide), legacy)
  for (nm in legacy) {
    # Cross differences (scen_1a_3 ...) subtract two ~70 values, so the stored
    # file's ~1e-14 absolute noise is a larger *relative* error on near-zero
    # elements; testthat compares the vector at mean-relative tolerance.
    expect_equal(
      wide[[nm]],
      rows20[[nm]],
      tolerance = 1e-12,
      label = paste0("run_scenarios column ", nm, " (rows 1-20)")
    )
  }
})

test_that("point estimates via run_scenarios equal plain life_add() (04)", {
  pt <- run_scenarios(tibble::tibble(.rows = 1), derive_sym_b = FALSE)
  d <- params_def
  plain <- list(
    baseline = list(),
    real_adh = list(
      p_NA = d$p_NA - 0.05,
      p_CA = d$p_CA + 0.1,
      p_Att = d$p_Att + 0.1
    ),
    opt_adh = list(
      p_NA = d$p_NA - 0.1,
      p_CA = d$p_CA + 0.15,
      p_Att = d$p_Att + 0.2
    ),
    diag_sens = list(
      Sens_0 = d$Sens_0 + 0.05,
      Sens_A = d$Sens_A + 0.05,
      Sens_B = d$Sens_B + 0.05
    ),
    treat = list(L_A = d$L_A + 1, L_B = d$L_B + 0.5, L_C = d$L_C + 0.5)
  )
  for (r in names(risk_multipliers)) {
    hcc <- d$p_HCC * c(base = 1, low_risk = 1 / 2, high_risk = 2)[[r]]
    for (s in names(plain)) {
      expect_equal(
        pt$ly[pt$scenario == s & pt$risk == r],
        do.call(life_add, c(list(p_HCC = hcc), plain[[s]])),
        tolerance = 1e-12,
        label = paste("point estimate", s, r)
      )
    }
  }
})

test_that("each scenario changes exactly its own parameters", {
  changed <- function(id) {
    new <- scenario_defs[[id]](params_def)
    names(params_def)[!mapply(identical, params_def, new)]
  }
  expect_identical(changed("baseline"), character(0))
  expect_setequal(changed("real_adh"), ADHERENCE_PARAMS)
  expect_setequal(changed("opt_adh"), ADHERENCE_PARAMS)
  expect_setequal(changed("diag_sens"), c("Sens_0", "Sens_A", "Sens_B"))
  expect_setequal(changed("treat"), c("L_A", "L_B", "L_C"))
})

test_that("scenario_table, defs and risk multipliers are consistent", {
  expect_identical(scenario_table$id, names(scenario_defs))
  expect_identical(names(risk_multipliers), c("base", "low_risk", "high_risk"))
  expect_equal(unname(risk_multipliers), c(1, 0.5, 2))
})
