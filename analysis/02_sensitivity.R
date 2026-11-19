# analysis/02_sensitivity.R
#
# Univariate sensitivity analyses for the HCC surveillance model.
# Reads Parameters.csv and Parameters_25.csv from data/.
# Writes figures to figures/ and records values in output/manuscript_values.csv.
#
# Figures written:
#   sens_anal_univariate.jpg
#   sens_anal_f4_0_upper.jpg
#   sens_anal_f4_0_lower.jpg
#   sens_anal_f4_0_both.jpg
#   sens_anal_mod.{jpg,png,pdf}
#   sens_anal_25_univariate.jpg
#   sens_anal_25_f4_0_upper.jpg
#   sens_anal_25_f4_0_lower.jpg

# ----- Initialise -------------------------------------------------------------
source(here::here("R", "setup.R"))
source(here::here("R", "model_core.R"))
source(here::here("R", "plotting.R"))

# Parameter groups and labels (the facet and axis labels) are in R/setup.R.
x_lab_mod <- "Additional Life Years Gained"
all_three <- c(jpg = 300, png = 600, pdf = 600)

# ----- Sensitivity data -------------------------------------------------------
param_tbl <- life_tibble(here("data", "Parameters.csv"))

# ----- Univariate sensitivity tornado -----------------------------------------
param_tbl |>
  order_parameters("years_mid_comp") |>
  add_group() |>
  dplyr::filter(years_mid_comp != 0) |>
  plot_tornado(
    "years_mid_comp",
    variable_labels,
    title = "Univariate Sensitivity Analysis",
    facet = "group"
  ) |>
  save_figure("sens_anal_univariate", 5.5, 8, "jpg", dpi = 300)

# ----- Upper and lower HCC-incidence graphs -----------------------------------
param_tbl |>
  plot_bound_tornado(
    "high",
    title = TeX(r'(Upper $\textit{p}_{HCC}$ Sensitivity Analysis)')
  ) |>
  save_figure("sens_anal_f4_0_upper", 5.5, 9, "jpg", dpi = 300)

param_tbl |>
  plot_bound_tornado(
    "low",
    title = TeX(r"(Lower $\textit{p}_{HCC}$ Sensitivity Analysis)")
  ) |>
  save_figure("sens_anal_f4_0_lower", 5.5, 9, "jpg", dpi = 300)

# ----- Combined upper and lower graph -----------------------------------------
param_tbl |>
  order_parameters("years_low_comp") |>
  add_group() |>
  dplyr::filter(years_low_comp != 0, Parameter != "p_HCC") |>
  dplyr::select(Parameter, Bound, years_low_comp, years_high_comp, Type) |>
  tidyr::pivot_longer(
    c("years_low_comp", "years_high_comp"),
    names_to = "years_comp"
  ) |>
  plot_tornado(
    "value",
    variable_labels,
    title = "Sensitivity Analysis",
    facet = "grid"
  ) |>
  save_figure("sens_anal_f4_0_both", 5.5, 8, "jpg", dpi = 300)

# ----- 25%-variation univariate sensitivity -----------------------------------
param_tbl_25 <- life_tibble(here("data", "Parameters_25.csv"))
param_chang_tbl <- modifiable_sensitivity(param_tbl_25)

# Ledger: change in life-years from a 25% change in each modifiable parameter
merged_names <- paste(
  param_chang_tbl$Parameter,
  param_chang_tbl$Bound,
  sep = "_"
)
record_values(setNames(param_chang_tbl$years_mid_comp, merged_names))

plot_tornado(
  param_chang_tbl,
  "years_mid_comp",
  x_lab = x_lab_mod,
  y_lab = NULL,
  facet = "group"
) |>
  save_figure("sens_anal_mod", 5, 4, dpi = all_three)

# ----- Full 25% univariate (all parameters) -----------------------------------
param_tbl_25 |>
  order_parameters("years_mid_comp") |>
  add_group() |>
  plot_tornado(
    "years_mid_comp",
    variable_labels,
    x_lab = x_lab_mod,
    y_lab = NULL,
    facet = "group"
  ) |>
  save_figure("sens_anal_25_univariate", 5, 8, "jpg", dpi = 300)

# ----- 25% upper and lower HCC-incidence graphs -------------------------------
param_tbl_25 |>
  plot_bound_tornado(
    "high",
    exclude = NON_MODIFIABLE_PARAMS,
    title = TeX(
      r"(Upper $\textit{p}_{HCC}$ Estimated change in mortality from 25% change in each modifiable parameter)"
    )
  ) |>
  save_figure("sens_anal_25_f4_0_upper", 5.5, 9, "jpg", dpi = 300)

param_tbl_25 |>
  plot_bound_tornado(
    "low",
    exclude = NON_MODIFIABLE_PARAMS,
    title = TeX(
      r"(Lower $\textit{p}_{HCC}$ Estimated change in mortality from 25% change\nin each modifiable parameter)"
    )
  ) |>
  save_figure("sens_anal_25_f4_0_lower", 5.5, 9, "jpg", dpi = 300)

# ----- Regression anchors (verified against current Parameters.csv) -----------
# life_add() = 68.02326 with current default parameters
stopifnot(isTRUE(all.equal(life_add(), 68.02326, tolerance = 0.01)))
message("Regression check passed: life_add() = ", round(life_add(), 4))
