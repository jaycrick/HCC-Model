# analysis/04_clinical_scenarios.R
#
# Clinical scenario comparisons under baseline, low-risk, and high-risk
# HCC incidence assumptions.
# Requires data-derived/resampled_pars_fin.rds (from 01_bootstrap_params.R)
# and data-derived/param_names.rds (updated by 02_sensitivity.R).
#
# Regression anchors (checked against original output):
#   noUS_scen ≈ 14.19459   (no ultrasound surveillance)
#   Baseline  ≈ 66.4987    (current practice, default parameters)
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

dir.create(here("figures"), showWarnings = FALSE)

# ----- Regression anchors -----------------------------------------------------
noUS_scen <- life_add(p_NA = 1, p_CA = 0) # ≈ 14.856
Baseline <- life_add() # ≈ 68.023

# Regression anchors — derived from the current Parameters.csv.
# (Note: comments "14.19459" and "66.4987" in the original Clinical Scenarios.R
# reflected a prior version of Parameters.csv; the values from the original
# Sensitivity_analyses.R:304 comment "68.02326" match the current file.)
stopifnot(isTRUE(all.equal(noUS_scen, 14.85677, tolerance = 0.01)))
stopifnot(isTRUE(all.equal(Baseline, 68.02326, tolerance = 0.01)))
message("Regression checks passed:")
message("  noUS_scen = ", round(noUS_scen, 5))
message("  Baseline  = ", round(Baseline, 5))

# ----- Point estimates for every scenario and risk level -----------------------
# Scenario definitions: R/scenarios.R. At default parameters Sym_B comes from
# params_def (as life_add() does), not derived from Sym_A.
point <- run_scenarios(tibble(.rows = 1), derive_sym_b = FALSE)

# ----- Load bootstrap outputs -------------------------------------------------
resampled_pars <- readRDS(here("data-derived", "resampled_pars_fin.rds"))

# ----- All scenarios (absolute life years) ------------------------------------
# Order: treatment, diagnosis, optimistic, realistic adherence, baseline.
all_tbl <- scenario_table[c(5, 4, 3, 2, 1), ]
violin_scenarios <- tibble(
  Scenario = rep(all_tbl$label, each = nrow(resampled_pars)),
  `Life Years per Cohort of 100` = unlist(
    resampled_pars[all_tbl$column],
    use.names = FALSE
  )
)

scenarios <- tibble(
  Scenario = all_tbl$label,
  `Life Years per Cohort of 100` = point$ly[point$risk == "base"][match(
    all_tbl$id,
    point$scenario[point$risk == "base"]
  )],
  type = ifelse(all_tbl$id == "baseline", "base", "scen"),
  lowCI = vapply(
    all_tbl$column,
    \(cl) as.double(quantile(resampled_pars[[cl]], 0.025)),
    numeric(1)
  ),
  highCI = vapply(
    all_tbl$column,
    \(cl) as.double(quantile(resampled_pars[[cl]], 0.975)),
    numeric(1)
  )
)

all_scens_plot <- ggplot(
  data = scenarios,
  aes(fct_inorder(Scenario), `Life Years per Cohort of 100`)
) +
  geom_col() +
  geom_violin(
    data = violin_scenarios,
    mapping = aes(Scenario, `Life Years per Cohort of 100`)
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_fill_viridis(discrete = TRUE) +
  labs(x = "Scenario") +
  theme(legend.position = "none")

all_scens_plot

ggsave(
  plot = all_scens_plot,
  here("figures", "scenarios_all.jpg"),
  height = 4,
  width = 6
)

# ----- Diminishing returns: adherence -----------------------------------------
`Partial Adherence Attendence` <- map_dbl(
  0:50 / 100,
  ~ life_add(p_Att = params_def$p_Att + .x)
)
`Non-Adherence` <- map_dbl(0:30 / 100, ~ life_add(p_NA = params_def$p_NA - .x))
`Complete Adherence` <- map_dbl(
  0:40 / 100,
  ~ life_add(p_CA = params_def$p_CA + .x)
)

adh_inc_tbl <- tibble(`Percent change in Category` = 0:50 / 100) %>%
  merge(
    tibble(
      `Percent change in Category` = 0:50 / 100,
      `Partial Adherence Attendence`
    ),
    all = TRUE
  ) %>%
  merge(
    tibble(`Percent change in Category` = 0:30 / 100, `Non-Adherence`),
    all = TRUE
  ) %>%
  merge(
    tibble(`Percent change in Category` = 0:40 / 100, `Complete Adherence`),
    all = TRUE
  ) %>%
  pivot_longer(
    !`Percent change in Category`,
    names_to = "Adherence Category",
    values_to = "Life Extended"
  )

plot_adh_inc <- ggplot(
  adh_inc_tbl,
  aes(
    `Percent change in Category`,
    `Life Extended`,
    color = `Adherence Category`
  )
) +
  geom_line(linewidth = 1.5) +
  theme_bw(base_family = "serif", base_size = 12) +
  ggtitle("Effect of Adherence Increase\non Intervention Outcomes") +
  scale_x_continuous(labels = label_percent(), expand = expansion(0)) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  )

plot_adh_inc

ggsave(
  here("figures", "increasing_Adh_plot_10_year.jpg"),
  plot_adh_inc,
  width = 6,
  height = 3
)

# ----- Diminishing returns: ultrasound sensitivity ---------------------------
`BCLC 0` <- map_dbl(0:50 / 100, ~ life_add(Sens_0 = params_def$Sens_0 + .x))
`BCLC A` <- map_dbl(0:30 / 100, ~ life_add(Sens_A = params_def$Sens_A + .x))
`BCLC B` <- map_dbl(0:15 / 100, ~ life_add(Sens_B = params_def$Sens_B + .x))

US_inc_tbl <- tibble(`Increase in Sensitivity` = 0:50 / 100) %>%
  merge(
    tibble(`Increase in Sensitivity` = 0:50 / 100, `BCLC 0`),
    all = TRUE
  ) %>%
  merge(
    tibble(`Increase in Sensitivity` = 0:30 / 100, `BCLC A`),
    all = TRUE
  ) %>%
  merge(
    tibble(`Increase in Sensitivity` = 0:15 / 100, `BCLC B`),
    all = TRUE
  ) %>%
  pivot_longer(
    !`Increase in Sensitivity`,
    names_to = "US Stage",
    values_to = "Life Extended"
  )

plot_US_inc <- ggplot(
  US_inc_tbl,
  aes(`Increase in Sensitivity`, `Life Extended`, color = `US Stage`)
) +
  geom_line(linewidth = 1.5) +
  theme_bw(base_family = "serif", base_size = 12) +
  ggtitle("Effect of US Sensitivities on Intervention Outcomes") +
  scale_x_continuous(labels = label_percent(), expand = expansion(0)) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  )

plot_US_inc

ggsave(
  here("figures", "increasing_US_Sens_plot_10_year.jpg"),
  plot_US_inc,
  width = 6,
  height = 4
)

# ----- Scenario improvements (difference from baseline), by risk level ---------
# Writes the point estimates and 95% limits to the ledger and saves the figure.
risk_figures <- tibble(
  risk = c("base", "low_risk", "high_risk"),
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

  param_names <- readRDS(here("data-derived", "param_names.rds"))
  key <- d$points$ledger_key
  param_names[paste0(key, sfx)] <- d$points[[2]]
  param_names[paste0(key, "_lower", sfx)] <- d$points$lowCI
  param_names[paste0(key, "_upper", sfx)] <- d$points$highCI
  saveRDS(param_names, file = here("data-derived", "param_names.rds"))

  fig <- plot_scenario_violins(d, risk_figures$xlab[k], risk_figures$ymax[k])
  print(fig)
  for (ext in c("jpg", "png", "pdf")) {
    ggsave(
      plot = fig,
      here("figures", paste0("scenarios_new", sfx, ".", ext)),
      width = 6,
      height = 4,
      dpi = if (ext == "jpg") 300 else 600
    )
  }
}
