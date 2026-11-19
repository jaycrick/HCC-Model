# Published-result anchors: Cumming et al. 2024, Cancers 16:2745.
# Values are checked against the shipped data-derived/resampled_pars_fin.rds
# and the current code, after rounding to the paper's printed precision.

rp <- resampled_pars_fin
base_la <- life_add()
q <- function(x, p) unname(quantile(x, p))

test_that("Abstract / Sec 3.3 scenario point estimates", {
  pts <- c(
    realistic = life_add(
      p_NA = params_def$p_NA - 0.05,
      p_CA = params_def$p_CA + 0.10,
      p_Att = params_def$p_Att + 0.10
    ),
    optimistic = life_add(
      p_NA = params_def$p_NA - 0.10,
      p_CA = params_def$p_CA + 0.15,
      p_Att = params_def$p_Att + 0.20
    ),
    diagnosis = life_add(
      Sens_0 = params_def$Sens_0 + 0.05,
      Sens_A = params_def$Sens_A + 0.05,
      Sens_B = params_def$Sens_B + 0.05
    ),
    treatment = life_add(
      L_A = params_def$L_A + 1,
      L_B = params_def$L_B + 0.5,
      L_C = params_def$L_C + 0.5
    )
  ) -
    base_la
  paper <- c(
    realistic = 9.8,
    optimistic = 17.2,
    diagnosis = 7.0,
    treatment = 9.0
  )
  for (nm in names(paper)) {
    expect_equal(
      round(pts[[nm]], 1),
      paper[[nm]],
      label = paste0(
        "Abstract/Sec 3.3 point estimate, scenario '",
        nm,
        "' (years per 100)"
      )
    )
  }
})

test_that("Abstract / Sec 3.3 scenario 95% CIs", {
  cols <- c(
    realistic = "scen_1a_diff",
    optimistic = "scen_1b_diff",
    diagnosis = "scen_2_diff",
    treatment = "scen_3_diff"
  )
  lo <- c(realistic = 7.9, optimistic = 13.9, diagnosis = 3.6, treatment = 7.5)
  hi <- c(
    realistic = 11.6,
    optimistic = 20.3,
    diagnosis = 13.8,
    treatment = 10.3
  )
  for (nm in names(cols)) {
    expect_equal(
      round(q(rp[[cols[[nm]]]], 0.025), 1),
      lo[[nm]],
      label = paste0("Abstract/Sec 3.3 95% CI lower, scenario '", nm, "'")
    )
    expect_equal(
      round(q(rp[[cols[[nm]]]], 0.975), 1),
      hi[[nm]],
      label = paste0("Abstract/Sec 3.3 95% CI upper, scenario '", nm, "'")
    )
  }
})

test_that("Table 4: percent of bootstrap rows where X_diff > Y_diff", {
  pct <- function(x, y) round(100 * mean(rp[[x]] > rp[[y]]), 1)
  expect_equal(
    pct("scen_1b_diff", "scen_2_diff"),
    98.4,
    label = "Table 4: optimistic > diagnosis"
  )
  expect_equal(
    pct("scen_1a_diff", "scen_2_diff"),
    80.1,
    label = "Table 4: realistic > diagnosis"
  )
  expect_equal(
    pct("scen_3_diff", "scen_2_diff"),
    74.7,
    label = "Table 4: treatment > diagnosis"
  )
  expect_equal(
    pct("scen_1b_diff", "scen_3_diff"),
    100.0,
    label = "Table 4: optimistic > treatment"
  )
  expect_equal(
    pct("scen_1a_diff", "scen_3_diff"),
    73.6,
    label = "Table 4: realistic > treatment"
  )
  expect_equal(
    pct("scen_1b_diff", "scen_1a_diff"),
    100.0,
    label = "Table 4: optimistic > realistic"
  )
})

test_that("Sec 3.1: mean months in each stage (Fig 2A)", {
  params <- read_params(csv_main)
  st <- STAGE_P_NAMES
  expect_equal(
    round(1 / params[st, "POINT"], 1),
    c(12.5, 18.0, 2.9, 2.8, 6.2),
    ignore_attr = TRUE,
    label = "Sec 3.1: mean months in stage = 1/p (POINT), stages 0,A,B,C,D"
  )
  ci_lo <- vapply(st, function(s) 1 / q(rp[[s]], 0.975), numeric(1))
  ci_hi <- vapply(st, function(s) 1 / q(rp[[s]], 0.025), numeric(1))
  expect_equal(
    unname(round(ci_lo, 1)),
    c(4.2, 12.2, 1.4, 1.7, 6.1),
    ignore_attr = TRUE,
    label = "Sec 3.1: lower CI of mean months = 1/q97.5(p), stages 0,A,B,C,D"
  )
  expect_equal(
    unname(round(ci_hi, 1)),
    c(32.6, 23.5, 4.9, 4.5, 7.8),
    ignore_attr = TRUE,
    label = "Sec 3.1: upper CI of mean months = 1/q2.5(p), stages 0,A,B,C,D"
  )
  expect_equal(
    round(100 * mean(rp$p_A < rp$p_0), 0),
    73,
    label = "Sec 3.1: stage A longer than stage 0 in 73% of bootstrap estimates"
  )
})

test_that("Sec 3.4: SE of Baseline and partial correlations", {
  expect_equal(
    round(sd(rp$Baseline), 1),
    5.7,
    label = "Sec 3.4: SE of Baseline = sd(Baseline)"
  )
  lg <- as.data.frame(lapply(
    rp[NLL_PAR_NAMES],
    qlogis
  ))
  lg$Baseline <- rp$Baseline
  pc <- ppcor::pcor(lg)$estimate
  expect_equal(
    round(pc["Baseline", "p_0"], 2),
    -0.87,
    label = "Sec 3.4: partial correlation Baseline ~ logit(p_0)"
  )
  expect_equal(
    round(pc["Baseline", "Sens_0"], 2),
    0.92,
    label = "Sec 3.4: partial correlation Baseline ~ logit(Sens_0)"
  )
})

test_that("Sec 3.2: +/-25% univariate sensitivity", {
  lt <- life_tibble(csv_25)
  get <- function(par, bound) {
    lt$years_mid_comp[lt$Parameter == par & lt$Bound == bound]
  }
  chk <- function(par, bound, val, dp) {
    expect_equal(
      round(get(par, bound), dp),
      val,
      label = paste0(
        "Sec 3.2 +/-25%: ",
        bound,
        " ",
        par,
        " (years_mid_comp, ",
        dp,
        " dp)"
      )
    )
  }
  chk("p_Att", "Upper", 3.5, 1)
  chk("p_NA", "Lower", 4.4, 1)
  chk("p_CA", "Upper", 3.6, 1)
  chk("Sens_0", "Upper", 6.6, 1)
  chk("Sens_A", "Upper", 1.7, 1)
  chk("Sens_B", "Upper", 0.2, 1)
  chk("Sens_C", "Upper", 0.005, 3)
  chk("L_0", "Upper", 10.2, 1)
  chk("L_A", "Upper", 5.7, 1)
  chk("L_B", "Upper", 0.95, 2)
  chk("L_C", "Upper", 0.2, 1)
})
