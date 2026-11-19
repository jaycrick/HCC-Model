# Scenarios as data (file R/scenarios.R)
# Intervention scenarios as data: the single source of truth for what each
# scenario changes, how it is labelled and what it is called in the legacy
# result columns. Requires R/setup.R and R/model_core.R.
#
# Definitions: Cumming et al. 2024, Cancers 16:2745, Sec. 2.4. Risk levels
# (low / high HCC incidence): Supplementary Figs. S6-S7.
# To change a scenario, edit its entry in `scenario_defs` (and nothing else).

# ----- Scenario definitions ---------------------------------------------------

#' Scenario modifiers: pure functions `params -> params`.
#'
#' Each modifier receives the full parameter list (defaults, or defaults with a
#' bootstrap draw applied) and returns it with the scenario's changes made.
#' Increments are relative to the incoming value, so a bootstrapped Sens_0 or
#' Sens_A is shifted, not replaced.
#'
#' - baseline:  current practice (no change)
#' - real_adh:  realistic adherence improvement (paper Sec. 2.4)
#' - opt_adh:   optimistic adherence improvement
#' - diag_sens: ultrasound sensitivity +0.05 for BCLC 0, A and B
#' - treat:     better treatment, i.e. longer life gained per stage at diagnosis
#'
#' NOTE (finding P6): the paper text lists only L_A and L_B for the treatment
#' scenario, but the code, and the published 9.0 life-years, also include
#' L_C + 0.5. L_C is kept here so published results reproduce.
scenario_defs <- list(
  baseline = function(p) p,
  real_adh = function(p) {
    modify_params(
      p,
      list(p_NA = p$p_NA - 0.05, p_CA = p$p_CA + 0.10, p_Att = p$p_Att + 0.10)
    )
  },
  opt_adh = function(p) {
    modify_params(
      p,
      list(p_NA = p$p_NA - 0.10, p_CA = p$p_CA + 0.15, p_Att = p$p_Att + 0.20)
    )
  },
  diag_sens = function(p) {
    modify_params(
      p,
      list(
        Sens_0 = p$Sens_0 + 0.05,
        Sens_A = p$Sens_A + 0.05,
        Sens_B = p$Sens_B + 0.05
      )
    )
  },
  treat = function(p) {
    modify_params(
      p,
      list(L_A = p$L_A + 1, L_B = p$L_B + 0.5, L_C = p$L_C + 0.5) # L_C: see P6
    )
  }
)

#' Multipliers on the annual HCC incidence p_HCC for the three risk levels
#' (Supplementary Figs. S6-S7): half (low risk) and double (high risk).
risk_multipliers <- c(base = 1, low_risk = 0.5, high_risk = 2)

#' The one mapping from internal scenario ids to everything legacy.
#'
#' - id:         name in `scenario_defs`
#' - column:     legacy column in data-derived/resampled_pars_fin.rds
#' - label:      axis label in the scenario figures
#' - ledger_key: prefix of the entries in data-derived/param_names.rds
#'               (NA for the baseline, which has no improvement to report)
scenario_table <- tibble::tribble(
  ~id,
  ~column,
  ~label,
  ~ledger_key,
  "baseline",
  "Baseline",
  "Current Scenario",
  NA,
  "real_adh",
  "scen_1a",
  "Realistic Adherence\nImprovements",
  "real_adh",
  "opt_adh",
  "scen_1b",
  "Optimistic Adherence\nImprovements",
  "opt_adh",
  "diag_sens",
  "scen_2",
  "Diagnosis Sensitivity\nImprovements",
  "diag_scen",
  "treat",
  "scen_3",
  "Treatment Improvements",
  "treat_scens"
)

#' Legacy cross differences between two improvement scenarios, `a - b`.
#' The legacy suffixes are confusing: "_3" means "vs scen_2" (diagnosis
#' sensitivity) and "_4" means "vs scen_3" (treatment).
cross_diff_table <- tibble::tribble(
  ~column,
  ~a,
  ~b,
  "scen_1a_3",
  "real_adh",
  "diag_sens",
  "scen_1a_4",
  "real_adh",
  "treat",
  "scen_1b_3",
  "opt_adh",
  "diag_sens",
  "scen_1b_4",
  "opt_adh",
  "treat"
)

# The table and the definitions must list the same scenarios in the same order.
stopifnot(identical(scenario_table$id, names(scenario_defs)))

# ----- Running scenarios ------------------------------------------------------

#' Life-years for every scenario x risk level x parameter draw.
#'
#' @param param_rows    Tibble of parameter draws, one row per draw; columns
#'                      that are not model parameters are ignored, and
#'                      parameters not given take their `params_def` value.
#' @param defs          Named list of scenario modifiers (`scenario_defs`);
#'                      must contain `baseline`.
#' @param risk          Named numeric multipliers on the draw's p_HCC
#'                      (`risk_multipliers`).
#' @param derive_sym_b  TRUE (default) for bootstrap draws: `prepare_draws()`
#'                      sets `Sym_B = SYM_B_RATIO * Sym_A` from the draw's
#'                      Sym_A, because Sym_B is not estimated. FALSE for point
#'                      estimates at default parameters, where Sym_B is taken
#'                      from `params_def` as `life_add()` does.
#' @return Long tibble with columns draw (row index), scenario, risk, ly (life
#'   years per cohort of 100) and diff (ly minus the baseline's ly for the same
#'   draw and risk level).
run_scenarios <- function(
  param_rows,
  defs = scenario_defs,
  risk = risk_multipliers,
  derive_sym_b = TRUE
) {
  stopifnot("defs needs a 'baseline' scenario" = "baseline" %in% names(defs))
  rows <- if (derive_sym_b) prepare_draws(param_rows) else draws(param_rows)
  one_draw <- function(i, mult, def) {
    p <- modify_params(params_def, as.list(rows[i, ]))
    p$p_HCC <- p$p_HCC * mult
    life_add(params = def(p))
  }
  out <- tidyr::expand_grid(
    draw = seq_len(nrow(param_rows)),
    scenario = names(defs),
    risk = names(risk)
  )
  out$ly <- purrr::pmap_dbl(
    list(out$draw, out$scenario, out$risk),
    function(d, s, r) one_draw(d, risk[[r]], defs[[s]])
  )
  base_ly <- out$ly[out$scenario == "baseline"]
  base_key <- paste(
    out$draw[out$scenario == "baseline"],
    out$risk[out$scenario == "baseline"]
  )
  out$diff <- out$ly - base_ly[match(paste(out$draw, out$risk), base_key)]
  out
}

# ----- Legacy wide format -----------------------------------------------------

#' Pivot `run_scenarios()` output to the legacy wide columns of
#' data-derived/resampled_pars_fin.rds.
#'
#' Per risk level (suffix "" for base, "_low_risk", "_high_risk"): one column
#' per scenario (`scenario_table$column`) and its `_diff` against baseline; the
#' base risk level additionally gets the cross differences of
#' `cross_diff_table` (scen_1a_3 etc.). Column order matches the legacy file.
#'
#' @param long  Output of `run_scenarios()`.
#' @return Tibble with one row per draw.
scenario_columns <- function(long) {
  col_of <- setNames(scenario_table$column, scenario_table$id)
  improvements <- setdiff(scenario_table$id, "baseline")
  risks <- unique(long$risk)
  suffix <- function(r) ifelse(r == "base", "", paste0("_", r))
  # one named list of vectors (per draw) for a risk level and a column of `long`
  by_scenario <- function(r, value, ids) {
    d <- long[long$risk == r, ]
    out <- lapply(ids, function(id) d[[value]][d$scenario == id])
    setNames(out, paste0(col_of[ids], if (value == "diff") "_diff", suffix(r)))
  }
  base_ly <- by_scenario("base", "ly", scenario_table$id)
  cross <- lapply(seq_len(nrow(cross_diff_table)), function(k) {
    ly <- setNames(base_ly, scenario_table$id)
    ly[[cross_diff_table$a[k]]] - ly[[cross_diff_table$b[k]]]
  })
  setNames(cross, cross_diff_table$column)
  cols <- c(
    unlist(lapply(risks, by_scenario, "ly", scenario_table$id), FALSE),
    by_scenario("base", "diff", improvements),
    setNames(cross, cross_diff_table$column),
    unlist(
      lapply(setdiff(risks, "base"), by_scenario, "diff", improvements),
      FALSE
    )
  )
  tibble::as_tibble(cols)
}
