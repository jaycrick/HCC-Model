# analysis/01_bootstrap_params.R
#
# Bootstrap parameter estimation for the HCC surveillance model.
#
# IMPORTANT: This script is computationally expensive (bootstrap sampling
# across 10,000 draws with parallel optimisation).  The pre-computed
# outputs are shipped in data-derived/ and do NOT need to be regenerated
# to reproduce the published figures — run_all.R skips this script by default.
#
# To regenerate from scratch (requires HPC-scale resources):
#   Rscript analysis/01_bootstrap_params.R
#
# Outputs written to data-derived/:
#   resampled_pars.rds     — 8 bootstrapped transition parameters (10 000 rows)
#   resampled_pars_fin.rds — as above plus scenario life-year columns
#
# The manuscript values are recorded to output/manuscript_values.csv via
# record_values() (R/setup.R); data-derived/param_names.rds is a frozen,
# read-only reference and is no longer written.
#
# The saveRDS() calls below that write resampled_pars*.rds are the ONLY writes
# to data-derived/ in the pipeline. They deliberately regenerate inputs and run
# only when this script is run directly (run_all.R --bootstrap).
#
# Sections: (1) estimate, (2) summarise, (3) figures. Reusable pieces live in
# R/bootstrap.R (likelihood, loaders) and R/plot_*.R (figure builders).
#
# Figures written to figures/:
#   Bootstrap_sample_ranges.{jpg,png,pdf}
#   HCC_Rates_Range_Plot.jpg
#   HCC_Rates_Hist.{jpg,png,pdf}
#   HCC_Rates_Violin.{jpg,png,pdf}
#   HCC_Time_Range_Plot.{jpg,png,pdf}
#   US_sensitivity_Cal.{jpg,png,pdf}
#   bootstrap_var_base.{jpg,png,pdf}
#   plot_par_cor_logit.{jpg,png,pdf}
#   plot_cor_base_logit.{jpg,png,pdf}
#   multipanel_cals.{png,pdf}

# ----- Initialise -------------------------------------------------------------
source(here::here("R", "setup.R"))
source(here::here("R", "model_core.R"))
source(here::here("R", "bootstrap.R"))
source(here::here("R", "scenarios.R"))
source(here::here("R", "plot_theme.R"))
source(here::here("R", "plot_sweeps.R"))
source(here::here("R", "plot_scenarios.R"))
source(here::here("R", "plot_calibration.R"))

library(hypogeom) # Written for this project by D. Pavlyshyn; see README

dir.create(here("figures"), showWarnings = FALSE)

# ===== 1. ESTIMATE =====
# Point estimates are the published values (data/Parameters.csv = paper
# Tables 1-2). The optimiser below re-derives them as a check; it never
# overwrites them (finding G7).
point_par <- unlist(params_def[NLL_PAR_NAMES])
fit <- fit_point_estimate(NM_START)
message(
  "Nelder-Mead vs CSV: max |diff| = ",
  signif(max(abs(fit$par - point_par)), 3),
  "; NLL ",
  signif(fit$value, 10),
  " vs ",
  signif(neg_log_lik(point_par), 10)
)

# The sampling and parallel optimisation ran on an HPC cluster; load its
# output. The two saveRDS() calls in this script regenerate inputs and run
# only under --bootstrap.
set.seed(314)
resampled_pars <- bootstrap_estimates()
saveRDS(resampled_pars, file = here("data-derived", "resampled_pars.rds"))
resampled_pars <- readRDS(here("data-derived", "resampled_pars.rds"))

set.seed(315)
resampled_pars <- bind_cols(
  resampled_pars,
  scenario_columns(run_scenarios(resampled_pars))
)
saveRDS(resampled_pars, file = here("data-derived", "resampled_pars_fin.rds"))
resampled_pars <- readRDS(here("data-derived", "resampled_pars_fin.rds"))

# ===== 2. SUMMARISE =====
rates_tbl <- rates_table(resampled_pars, point_par)
mean_ids <- paste0("mean_", STAGE_IDS)
scen_diff <- \(a, b) mean(resampled_pars[[a]] > resampled_pars[[b]])

record_values(c(
  muA_gt_mu0 = mean(resampled_pars$p_A < resampled_pars$p_0),
  setNames(1 / rates_tbl$Point / 12, mean_ids),
  setNames(1 / rates_tbl$Upper / 12, paste0(mean_ids, "_Lower")),
  setNames(1 / rates_tbl$Lower / 12, paste0(mean_ids, "_Upper")),
  opt_vs_diag = scen_diff("scen_1b_diff", "scen_2_diff"),
  opt_vs_treat = scen_diff("scen_1b_diff", "scen_3_diff"),
  opt_vs_real = scen_diff("scen_1b_diff", "scen_1a_diff"),
  real_vs_diag = scen_diff("scen_1a_diff", "scen_2_diff"),
  real_vs_treat = scen_diff("scen_1a_diff", "scen_3_diff"),
  treat_vs_diag = scen_diff("scen_3_diff", "scen_2_diff"),
  sd_noint = sd(resampled_pars$Baseline),
  mean_noint = mean(resampled_pars$Baseline)
))

# Partial correlations of the logit parameters with life years gained.
pars_wide <- resampled_pars |> dplyr::select(1:8, Baseline)
pars_logit <- pars_wide |> mutate(across(-Baseline, stats::qlogis))
pcor_logit <- ppcor::pcor(pars_logit)$estimate
record_values(c(
  pcor_logit_0 = unname(pcor_logit["Baseline", "p_0"]),
  pcor_logit_US_0 = unname(pcor_logit["Baseline", "Sens_0"]),
  pcor_logit_US_A = unname(pcor_logit["Baseline", "Sens_A"])
))

# Results table: bootstrap median and 95% limits; "point" is replaced by the
# point-estimate scenarios (base risk, Sym_B derived from Sym_A) where they
# exist.
results_tbl <- resampled_pars |>
  pivot_longer(everything(), names_to = "result") |>
  group_by(result) |>
  summarise(
    point = median(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975)
  )
point_wide <- scenario_columns(run_scenarios(
  as_tibble_row(as.list(point_par)),
  risk = risk_multipliers["base"]
))
scen_ids <- c("1a", "1b", "2", "3")
point_results <- c(
  paste0("scen_", scen_ids, "_diff"),
  paste0("scen_", scen_ids),
  "scen_1a_3",
  "scen_1a_4",
  "scen_1b_3",
  "scen_1b_4"
)
results_tbl <- results_tbl |>
  rows_update(
    tibble(
      result = point_results,
      point = unlist(point_wide[1, point_results], use.names = FALSE)
    ),
    by = "result"
  )

if (interactive()) {
  print(table(resampled_pars$p_0 < resampled_pars$p_A))
  print(summary(resampled_pars$Sens_0))
  print(summary(resampled_pars$Sens_A))
  print(summary(resampled_pars$scen_1a_diff))
  print(quantile(resampled_pars$scen_3_diff, c(0.025, 0.975)))
  print(summary.lm(lm(Baseline ~ ., data = pars_logit))$r.squared)
  pairs(head(pars_logit, 1000))
  corrplot::corrplot(cor(pars_logit), method = "color", type = "lower")
  print(results_tbl, n = Inf)
}

# ===== 3. FIGURES =====
save_all <- function(plot, name, w, h) {
  save_figure(plot, name, w, h, formats = FIGURE_FORMATS, dpi = FIGURE_DPI_ALL)
}

save_all(
  plot_bootstrap_counts(bootstrap_inputs_long()),
  "Bootstrap_sample_ranges",
  6,
  4
)
save_figure(
  plot_rates_distribution(resampled_pars, rates_tbl, "box"),
  "HCC_Rates_Range_Plot",
  6,
  3,
  "jpg",
  dpi = 300
)
save_all(plot_rates_hist(resampled_pars, rates_tbl), "HCC_Rates_Hist", 6, 7)
save_all(
  plot_rates_distribution(resampled_pars, rates_tbl, "violin"),
  "HCC_Rates_Violin",
  6,
  4
)
HCC_Time_Range_Plot <- plot_time_range(resampled_pars, rates_tbl)
save_all(HCC_Time_Range_Plot, "HCC_Time_Range_Plot", 6, 4)

set.seed(315)
US_sensitivity_Cal <- plot_sensitivity_hex(resampled_pars, point_par)
save_all(US_sensitivity_Cal, "US_sensitivity_Cal", 6, 3.5)
save_all(plot_baseline_hist(resampled_pars), "bootstrap_var_base", 6, 3)

save_figure(
  plot_pcor_heatmap(pcor_logit),
  "plot_par_cor_logit",
  6,
  5,
  "jpg",
  dpi = 300
)
plot_cor_base_logit <- plot_pcor_baseline(pcor_logit)
save_all(plot_cor_base_logit, "plot_cor_base_logit", 6, 3)

# Multipanel: calibration plots plus the 25%-variation univariate sensitivity.
plot_univar_mod <- plot_tornado(
  modifiable_sensitivity(life_tibble(here("data", "Parameters_25.csv"))),
  "years_mid_comp",
  x_lab = "Additional Life Years Gained",
  y_lab = NULL,
  facet = "group"
)
multipanel_cals <- cowplot::plot_grid(
  HCC_Time_Range_Plot +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)),
  US_sensitivity_Cal + theme(legend.position = "bottom"),
  plot_cor_base_logit +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      legend.position = "bottom"
    ),
  plot_univar_mod,
  labels = "AUTO",
  rel_widths = c(1, 2)
)
save_figure(
  multipanel_cals,
  "multipanel_cals",
  6,
  8,
  formats = c("png", "pdf"),
  dpi = c(png = 600, pdf = 600)
)
