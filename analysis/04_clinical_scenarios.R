# analysis/04_clinical_scenarios.R
#
# Clinical scenario comparisons under baseline, low-risk, and high-risk
# HCC incidence assumptions.
# Requires data-derived/resampled_pars_fin.rds (from 01_bootstrap_params.R)
# (read-only).
#
# Figures written to figures/:
#   scenarios_all.jpg
#   scenarios_new.{jpg,png,pdf}
#   increasing_Adh_plot_10_year.jpg
#   increasing_US_Sens_plot_10_year.jpg
#   scenarios_new_low_risk.{jpg,png,pdf}
#   scenarios_new_high_risk.{jpg,png,pdf}

# ----- Initialise -------------------------------------------------------------
source(here::here("R", "setup.R"))
source(here::here("R", "model_core.R"))
source(here::here("R", "scenarios.R"))
source(here::here("R", "plot_theme.R"))
source(here::here("R", "plot_sweeps.R"))
source(here::here("R", "plot_scenarios.R"))
source(here::here("R", "plot_calibration.R"))

# ----- Regression anchors -----------------------------------------------------
noUS_scen <- life_add(p_NA = 1, p_CA = 0) # ≈ 14.856
Baseline <- life_add() # ≈ 68.023

# Regression anchors, derived from the current Parameters.csv. (The comments
# "14.19459" and "66.4987" in the original Clinical Scenarios.R reflected a
# prior version of the file; 68.02326 matches the current one.)
stopifnot(isTRUE(all.equal(noUS_scen, 14.85677, tolerance = 0.01)))
stopifnot(isTRUE(all.equal(Baseline, 68.02326, tolerance = 0.01)))
message("Regression checks passed:")
message("  noUS_scen = ", round(noUS_scen, 5))
message("  Baseline  = ", round(Baseline, 5))

# ----- Point estimates for every scenario and risk level ----------------------
# Scenario definitions: R/scenarios.R. At default parameters Sym_B comes from
# params_def (as life_add() does), not derived from Sym_A.
point <- run_scenarios(tibble(.rows = 1), derive_sym_b = FALSE)
resampled_pars <- readRDS(here("data-derived", "resampled_pars_fin.rds"))

# ----- All scenarios (absolute life years) ------------------------------------
plot_scenarios_all(point, resampled_pars) |>
  save_figure("scenarios_all", 6, 4, "jpg", dpi = 300)

# ----- Diminishing returns: adherence -----------------------------------------
list(
  list(
    name = "Partial Adherence Attendence",
    param = "p_Att",
    max = 50,
    sign = 1
  ),
  list(name = "Non-Adherence", param = "p_NA", max = 30, sign = -1),
  list(name = "Complete Adherence", param = "p_CA", max = 40, sign = 1)
) |>
  response_curves() |>
  plot_response_curves(
    "Effect of Adherence Increase\non Intervention Outcomes",
    "Percent change in Category",
    "Adherence Category"
  ) |>
  save_figure("increasing_Adh_plot_10_year", 6, 3, "jpg", dpi = 300)

# ----- Diminishing returns: ultrasound sensitivity ----------------------------
list(
  list(name = "BCLC 0", param = "Sens_0", max = 50, sign = 1),
  list(name = "BCLC A", param = "Sens_A", max = 30, sign = 1),
  list(name = "BCLC B", param = "Sens_B", max = 15, sign = 1)
) |>
  response_curves() |>
  plot_response_curves(
    "Effect of US Sensitivities on Intervention Outcomes",
    "Increase in Sensitivity",
    "US Stage"
  ) |>
  save_figure("increasing_US_Sens_plot_10_year", 6, 4, "jpg", dpi = 300)

# ----- Scenario improvements (difference from baseline), by risk level --------
# Writes the point estimates and 95% limits to the ledger and saves the figure.
risk_figures <- tibble(
  risk = names(risk_multipliers),
  xlab = c(
    "Scenario",
    "Scenario Under Low HCC Risk Assumption",
    "Scenario Under High HCC Risk Assumption"
  ),
  ymax = c(30, 15, 60)
)

for (k in seq_len(nrow(risk_figures))) {
  risk <- risk_figures$risk[k]
  sfx <- ifelse(risk == "base", "", paste0("_", risk))
  d <- scenario_diff_data(point, resampled_pars, risk)

  key <- d$points$ledger_key
  record_values(setNames(
    c(d$points[[2]], d$points$lowCI, d$points$highCI),
    c(paste0(key, sfx), paste0(key, "_lower", sfx), paste0(key, "_upper", sfx))
  ))

  plot_scenario_violins(d, risk_figures$xlab[k], risk_figures$ymax[k]) |>
    save_figure(paste0("scenarios_new", sfx), 6, 4, dpi = FIGURE_DPI_ALL)
}
