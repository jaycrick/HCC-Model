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

# ----- Scenario 1: Realistic adherence improvement ---------------------------
scen_1a <- life_add(
  p_NA = p_NA_def - 0.05,
  p_CA = p_CA_def + 0.10,
  p_Att = p_Att_def + 0.10,
  Sens_0 = Sens_0_def,
  Sens_A = Sens_A_def,
  Sens_B = Sens_B_def,
  Sens_C = Sens_C_def,
  Sens_D = Sens_D_def,
  L_0 = L_0_def,
  L_A = L_A_def,
  L_B = L_B_def,
  L_C = L_C_def,
  L_D = L_D_def
)

# ----- Scenario 2: Optimistic adherence improvement --------------------------
scen_1b <- life_add(
  p_NA = p_NA_def - 0.10,
  p_CA = p_CA_def + 0.15,
  p_Att = p_Att_def + 0.20,
  Sens_0 = Sens_0_def,
  Sens_A = Sens_A_def,
  Sens_B = Sens_B_def,
  Sens_C = Sens_C_def,
  Sens_D = Sens_D_def,
  L_0 = L_0_def,
  L_A = L_A_def,
  L_B = L_B_def,
  L_C = L_C_def,
  L_D = L_D_def
)

# ----- Scenario 3: Improved ultrasound sensitivity ---------------------------
scen_2 <- life_add(
  p_NA = p_NA_def,
  p_CA = p_CA_def,
  p_Att = p_Att_def,
  Sens_0 = Sens_0_def + 0.05,
  Sens_A = Sens_A_def + 0.05,
  Sens_B = Sens_B_def + 0.05,
  Sens_C = Sens_C_def,
  Sens_D = Sens_D_def,
  L_0 = L_0_def,
  L_A = L_A_def,
  L_B = L_B_def,
  L_C = L_C_def,
  L_D = L_D_def
)

# ----- Scenario 4: Better treatments -----------------------------------------
scen_3 <- life_add(
  p_NA = p_NA_def,
  p_CA = p_CA_def,
  p_Att = p_Att_def,
  Sens_0 = Sens_0_def,
  Sens_A = Sens_A_def,
  Sens_B = Sens_B_def,
  Sens_C = Sens_C_def,
  Sens_D = Sens_D_def,
  L_0 = L_0_def,
  L_A = L_A_def + 1,
  L_B = L_B_def + 0.5,
  L_C = L_C_def + 0.5,
  L_D = L_D_def
)

# ----- Load bootstrap outputs -------------------------------------------------
resampled_pars <- readRDS(here("data-derived", "resampled_pars_fin.rds"))

# ----- All scenarios (absolute life years) ------------------------------------
violin_scenarios <- tribble(
  ~Scenario,
  ~`Life Years per Cohort of 100`,
  "Treatment Improvements",
  resampled_pars$scen_3,
  "Diagnosis Sensitivity\nImprovements",
  resampled_pars$scen_2,
  "Optimistic Adherence\nImprovements",
  resampled_pars$scen_1b,
  "Realistic Adherence\nImprovements",
  resampled_pars$scen_1a,
  "Current Scenario",
  resampled_pars$Baseline
) %>%
  unnest(`Life Years per Cohort of 100`)

scenarios <- tribble(
  ~Scenario,
  ~`Life Years per Cohort of 100`,
  ~type,
  ~lowCI,
  ~highCI,
  "Treatment Improvements",
  scen_3,
  "scen",
  as.double(quantile(resampled_pars$scen_3, 0.025)),
  as.double(quantile(resampled_pars$scen_3, 0.975)),
  "Diagnosis Sensitivity\nImprovements",
  scen_2,
  "scen",
  as.double(quantile(resampled_pars$scen_2, 0.025)),
  as.double(quantile(resampled_pars$scen_2, 0.975)),
  "Optimistic Adherence\nImprovements",
  scen_1b,
  "scen",
  as.double(quantile(resampled_pars$scen_1b, 0.025)),
  as.double(quantile(resampled_pars$scen_1b, 0.975)),
  "Realistic Adherence\nImprovements",
  scen_1a,
  "scen",
  as.double(quantile(resampled_pars$scen_1a, 0.025)),
  as.double(quantile(resampled_pars$scen_1a, 0.975)),
  "Current Scenario",
  Baseline,
  "base",
  as.double(quantile(resampled_pars$Baseline, 0.025)),
  as.double(quantile(resampled_pars$Baseline, 0.975))
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

# ----- Scenario improvements (difference from baseline) ----------------------
new_scenarios <- tribble(
  ~Scenario,
  ~`Life Years per Cohort of 100`,
  ~type,
  ~lowCI,
  ~highCI,
  "Treatment Improvements",
  scen_3 - Baseline,
  "Point\nEstimate",
  as.double(quantile(resampled_pars$scen_3_diff, 0.025)),
  as.double(quantile(resampled_pars$scen_3_diff, 0.975)),
  "Diagnosis Sensitivity\nImprovements",
  scen_2 - Baseline,
  "Point\nEstimate",
  as.double(quantile(resampled_pars$scen_2_diff, 0.025)),
  as.double(quantile(resampled_pars$scen_2_diff, 0.975)),
  "Optimistic Adherence\nImprovements",
  scen_1b - Baseline,
  "Point\nEstimate",
  as.double(quantile(resampled_pars$scen_1b_diff, 0.025)),
  as.double(quantile(resampled_pars$scen_1b_diff, 0.975)),
  "Realistic Adherence\nImprovements",
  scen_1a - Baseline,
  "Point\nEstimate",
  as.double(quantile(resampled_pars$scen_1a_diff, 0.025)),
  as.double(quantile(resampled_pars$scen_1a_diff, 0.975))
)

param_names <- readRDS(here("data-derived", "param_names.rds"))
param_names[c("treat_scens", "diag_scen", "opt_adh", "real_adh")] <-
  new_scenarios[[2]]
param_names[c(
  "treat_scens_lower",
  "diag_scen_lower",
  "opt_adh_lower",
  "real_adh_lower"
)] <- new_scenarios[[4]]
param_names[c(
  "treat_scens_upper",
  "diag_scen_upper",
  "opt_adh_upper",
  "real_adh_upper"
)] <- new_scenarios[[5]]
saveRDS(param_names, file = here("data-derived", "param_names.rds"))

new_violin_scenarios <- tribble(
  ~Scenario,
  ~`Life Years per Cohort of 100`,
  "Treatment Improvements",
  resampled_pars$scen_3_diff,
  "Diagnosis Sensitivity\nImprovements",
  resampled_pars$scen_2_diff,
  "Optimistic Adherence\nImprovements",
  resampled_pars$scen_1b_diff,
  "Realistic Adherence\nImprovements",
  resampled_pars$scen_1a_diff
) %>%
  unnest(`Life Years per Cohort of 100`)

new_scens_plot <- ggplot(
  data = new_violin_scenarios,
  mapping = aes(Scenario, `Life Years per Cohort of 100`)
) +
  geom_violin(trim = TRUE) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_x_discrete(
    limits = c(
      "Realistic Adherence\nImprovements",
      "Optimistic Adherence\nImprovements",
      "Diagnosis Sensitivity\nImprovements",
      "Treatment Improvements"
    )
  ) +
  labs(x = "Scenario") +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)
  ) +
  geom_point(
    data = new_scenarios,
    aes(Scenario, `Life Years per Cohort of 100`, shape = type),
    size = 3,
    show.legend = TRUE
  ) +
  scale_y_continuous(limits = c(0, 30), expand = expansion(0))

new_scens_plot

ggsave(
  plot = new_scens_plot,
  here("figures", "scenarios_new.jpg"),
  width = 6,
  height = 4
)
ggsave(
  plot = new_scens_plot,
  here("figures", "scenarios_new.png"),
  width = 6,
  height = 4,
  dpi = 600
)
ggsave(
  plot = new_scens_plot,
  here("figures", "scenarios_new.pdf"),
  width = 6,
  height = 4,
  dpi = 600
)

# ----- Diminishing returns: adherence -----------------------------------------
`Partial Adherence Attendence` <- map_dbl(
  0:50 / 100,
  ~ life_add(p_Att = p_Att_def + .x)
)
`Non-Adherence` <- map_dbl(0:30 / 100, ~ life_add(p_NA = p_NA_def - .x))
`Complete Adherence` <- map_dbl(0:40 / 100, ~ life_add(p_CA = p_CA_def + .x))

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
`BCLC 0` <- map_dbl(0:50 / 100, ~ life_add(Sens_0 = Sens_0_def + .x))
`BCLC A` <- map_dbl(0:30 / 100, ~ life_add(Sens_A = Sens_A_def + .x))
`BCLC B` <- map_dbl(0:15 / 100, ~ life_add(Sens_B = Sens_B_def + .x))

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

# ----- Scenarios under alternative HCC risk assumptions ----------------------
Baseline_low_risk <- life_add(p_HCC = p_HCC_def / 2)
Baseline_high_risk <- life_add(p_HCC = p_HCC_def * 2)

# Low-risk scenario calculations
scen_1a_low_risk <- life_add(
  p_HCC = p_HCC_def / 2,
  p_NA = p_NA_def - 0.05,
  p_CA = p_CA_def + 0.1,
  p_Att = p_Att_def + 0.1,
  Sens_0 = Sens_0_def,
  Sens_A = Sens_A_def,
  Sens_B = Sens_B_def,
  Sens_C = Sens_C_def,
  Sens_D = Sens_D_def,
  L_0 = L_0_def,
  L_A = L_A_def,
  L_B = L_B_def,
  L_C = L_C_def,
  L_D = L_D_def
)

scen_1b_low_risk <- life_add(
  p_HCC = p_HCC_def / 2,
  p_NA = p_NA_def - 0.1,
  p_CA = p_CA_def + 0.15,
  p_Att = p_Att_def + 0.2,
  Sens_0 = Sens_0_def,
  Sens_A = Sens_A_def,
  Sens_B = Sens_B_def,
  Sens_C = Sens_C_def,
  Sens_D = Sens_D_def,
  L_0 = L_0_def,
  L_A = L_A_def,
  L_B = L_B_def,
  L_C = L_C_def,
  L_D = L_D_def
)

scen_2_low_risk <- life_add(
  p_HCC = p_HCC_def / 2,
  p_NA = p_NA_def,
  p_CA = p_CA_def,
  p_Att = p_Att_def,
  Sens_0 = Sens_0_def + 0.05,
  Sens_A = Sens_A_def + 0.05,
  Sens_B = Sens_B_def + 0.05,
  Sens_C = Sens_C_def,
  Sens_D = Sens_D_def,
  L_0 = L_0_def,
  L_A = L_A_def,
  L_B = L_B_def,
  L_C = L_C_def,
  L_D = L_D_def
)

scen_3_low_risk <- life_add(
  p_HCC = p_HCC_def / 2,
  p_NA = p_NA_def,
  p_CA = p_CA_def,
  p_Att = p_Att_def,
  Sens_0 = Sens_0_def,
  Sens_A = Sens_A_def,
  Sens_B = Sens_B_def,
  Sens_C = Sens_C_def,
  Sens_D = Sens_D_def,
  L_0 = L_0_def,
  L_A = L_A_def + 1,
  L_B = L_B_def + 0.5,
  L_C = L_C_def + 0.5,
  L_D = L_D_def
)

# High-risk scenario calculations
scen_1a_high_risk <- life_add(
  p_HCC = p_HCC_def * 2,
  p_NA = p_NA_def - 0.05,
  p_CA = p_CA_def + 0.1,
  p_Att = p_Att_def + 0.1,
  Sens_0 = Sens_0_def,
  Sens_A = Sens_A_def,
  Sens_B = Sens_B_def,
  Sens_C = Sens_C_def,
  Sens_D = Sens_D_def,
  L_0 = L_0_def,
  L_A = L_A_def,
  L_B = L_B_def,
  L_C = L_C_def,
  L_D = L_D_def
)

scen_1b_high_risk <- life_add(
  p_HCC = p_HCC_def * 2,
  p_NA = p_NA_def - 0.1,
  p_CA = p_CA_def + 0.15,
  p_Att = p_Att_def + 0.2,
  Sens_0 = Sens_0_def,
  Sens_A = Sens_A_def,
  Sens_B = Sens_B_def,
  Sens_C = Sens_C_def,
  Sens_D = Sens_D_def,
  L_0 = L_0_def,
  L_A = L_A_def,
  L_B = L_B_def,
  L_C = L_C_def,
  L_D = L_D_def
)

scen_2_high_risk <- life_add(
  p_HCC = p_HCC_def * 2,
  p_NA = p_NA_def,
  p_CA = p_CA_def,
  p_Att = p_Att_def,
  Sens_0 = Sens_0_def + 0.05,
  Sens_A = Sens_A_def + 0.05,
  Sens_B = Sens_B_def + 0.05,
  Sens_C = Sens_C_def,
  Sens_D = Sens_D_def,
  L_0 = L_0_def,
  L_A = L_A_def,
  L_B = L_B_def,
  L_C = L_C_def,
  L_D = L_D_def
)

scen_3_high_risk <- life_add(
  p_HCC = p_HCC_def * 2,
  p_NA = p_NA_def,
  p_CA = p_CA_def,
  p_Att = p_Att_def,
  Sens_0 = Sens_0_def,
  Sens_A = Sens_A_def,
  Sens_B = Sens_B_def,
  Sens_C = Sens_C_def,
  Sens_D = Sens_D_def,
  L_0 = L_0_def,
  L_A = L_A_def + 1,
  L_B = L_B_def + 0.5,
  L_C = L_C_def + 0.5,
  L_D = L_D_def
)

# ----- Low-risk scenario plot ------------------------------------------------
new_scenarios_low_risk <- tribble(
  ~Scenario,
  ~`Life Years per Cohort of 100`,
  ~type,
  ~lowCI,
  ~highCI,
  "Treatment Improvements",
  scen_3_low_risk - Baseline_low_risk,
  "Point\nEstimate",
  as.double(quantile(resampled_pars$scen_3_diff_low_risk, 0.025)),
  as.double(quantile(resampled_pars$scen_3_diff_low_risk, 0.975)),
  "Diagnosis Sensitivity\nImprovements",
  scen_2_low_risk - Baseline_low_risk,
  "Point\nEstimate",
  as.double(quantile(resampled_pars$scen_2_diff_low_risk, 0.025)),
  as.double(quantile(resampled_pars$scen_2_diff_low_risk, 0.975)),
  "Optimistic Adherence\nImprovements",
  scen_1b_low_risk - Baseline_low_risk,
  "Point\nEstimate",
  as.double(quantile(resampled_pars$scen_1b_diff_low_risk, 0.025)),
  as.double(quantile(resampled_pars$scen_1b_diff_low_risk, 0.975)),
  "Realistic Adherence\nImprovements",
  scen_1a_low_risk - Baseline_low_risk,
  "Point\nEstimate",
  as.double(quantile(resampled_pars$scen_1a_diff_low_risk, 0.025)),
  as.double(quantile(resampled_pars$scen_1a_diff_low_risk, 0.975))
)

param_names <- readRDS(here("data-derived", "param_names.rds"))
param_names[c(
  "treat_scens_low_risk",
  "diag_scen_low_risk",
  "opt_adh_low_risk",
  "real_adh_low_risk"
)] <- new_scenarios_low_risk[[2]]
param_names[c(
  "treat_scens_lower_low_risk",
  "diag_scen_lower_low_risk",
  "opt_adh_lower_low_risk",
  "real_adh_lower_low_risk"
)] <- new_scenarios_low_risk[[4]]
param_names[c(
  "treat_scens_upper_low_risk",
  "diag_scen_upper_low_risk",
  "opt_adh_upper_low_risk",
  "real_adh_upper_low_risk"
)] <- new_scenarios_low_risk[[5]]
saveRDS(param_names, file = here("data-derived", "param_names.rds"))

new_violin_scenarios_low_risk <- tribble(
  ~Scenario,
  ~`Life Years per Cohort of 100`,
  "Treatment Improvements",
  resampled_pars$scen_3_diff_low_risk,
  "Diagnosis Sensitivity\nImprovements",
  resampled_pars$scen_2_diff_low_risk,
  "Optimistic Adherence\nImprovements",
  resampled_pars$scen_1b_diff_low_risk,
  "Realistic Adherence\nImprovements",
  resampled_pars$scen_1a_diff_low_risk
) %>%
  unnest(`Life Years per Cohort of 100`)

new_scens_plot_low_risk <- ggplot(
  data = new_violin_scenarios_low_risk,
  mapping = aes(Scenario, `Life Years per Cohort of 100`)
) +
  geom_violin(trim = TRUE) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_x_discrete(
    limits = c(
      "Realistic Adherence\nImprovements",
      "Optimistic Adherence\nImprovements",
      "Diagnosis Sensitivity\nImprovements",
      "Treatment Improvements"
    )
  ) +
  labs(x = "Scenario Under Low HCC Risk Assumption") +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)
  ) +
  geom_point(
    data = new_scenarios_low_risk,
    aes(Scenario, `Life Years per Cohort of 100`, shape = type),
    size = 3,
    show.legend = TRUE
  ) +
  scale_y_continuous(limits = c(0, 15), expand = expansion(0))

new_scens_plot_low_risk

ggsave(
  plot = new_scens_plot_low_risk,
  here("figures", "scenarios_new_low_risk.jpg"),
  width = 6,
  height = 4
)
ggsave(
  plot = new_scens_plot_low_risk,
  here("figures", "scenarios_new_low_risk.png"),
  width = 6,
  height = 4,
  dpi = 600
)
ggsave(
  plot = new_scens_plot_low_risk,
  here("figures", "scenarios_new_low_risk.pdf"),
  width = 6,
  height = 4,
  dpi = 600
)

# ----- High-risk scenario plot -----------------------------------------------
new_scenarios_high_risk <- tribble(
  ~Scenario,
  ~`Life Years per Cohort of 100`,
  ~type,
  ~lowCI,
  ~highCI,
  "Treatment Improvements",
  scen_3_high_risk - Baseline_high_risk,
  "Point\nEstimate",
  as.double(quantile(resampled_pars$scen_3_diff_high_risk, 0.025)),
  as.double(quantile(resampled_pars$scen_3_diff_high_risk, 0.975)),
  "Diagnosis Sensitivity\nImprovements",
  scen_2_high_risk - Baseline_high_risk,
  "Point\nEstimate",
  as.double(quantile(resampled_pars$scen_2_diff_high_risk, 0.025)),
  as.double(quantile(resampled_pars$scen_2_diff_high_risk, 0.975)),
  "Optimistic Adherence\nImprovements",
  scen_1b_high_risk - Baseline_high_risk,
  "Point\nEstimate",
  as.double(quantile(resampled_pars$scen_1b_diff_high_risk, 0.025)),
  as.double(quantile(resampled_pars$scen_1b_diff_high_risk, 0.975)),
  "Realistic Adherence\nImprovements",
  scen_1a_high_risk - Baseline_high_risk,
  "Point\nEstimate",
  as.double(quantile(resampled_pars$scen_1a_diff_high_risk, 0.025)),
  as.double(quantile(resampled_pars$scen_1a_diff_high_risk, 0.975))
)

param_names <- readRDS(here("data-derived", "param_names.rds"))
param_names[c(
  "treat_scens_high_risk",
  "diag_scen_high_risk",
  "opt_adh_high_risk",
  "real_adh_high_risk"
)] <- new_scenarios_high_risk[[2]]
param_names[c(
  "treat_scens_lower_high_risk",
  "diag_scen_lower_high_risk",
  "opt_adh_lower_high_risk",
  "real_adh_lower_high_risk"
)] <- new_scenarios_high_risk[[4]]
param_names[c(
  "treat_scens_upper_high_risk",
  "diag_scen_upper_high_risk",
  "opt_adh_upper_high_risk",
  "real_adh_upper_high_risk"
)] <- new_scenarios_high_risk[[5]]
saveRDS(param_names, file = here("data-derived", "param_names.rds"))

new_violin_scenarios_high_risk <- tribble(
  ~Scenario,
  ~`Life Years per Cohort of 100`,
  "Treatment Improvements",
  resampled_pars$scen_3_diff_high_risk,
  "Diagnosis Sensitivity\nImprovements",
  resampled_pars$scen_2_diff_high_risk,
  "Optimistic Adherence\nImprovements",
  resampled_pars$scen_1b_diff_high_risk,
  "Realistic Adherence\nImprovements",
  resampled_pars$scen_1a_diff_high_risk
) %>%
  unnest(`Life Years per Cohort of 100`)

new_scens_plot_high_risk <- ggplot(
  data = new_violin_scenarios_high_risk,
  mapping = aes(Scenario, `Life Years per Cohort of 100`)
) +
  geom_violin(trim = TRUE) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_x_discrete(
    limits = c(
      "Realistic Adherence\nImprovements",
      "Optimistic Adherence\nImprovements",
      "Diagnosis Sensitivity\nImprovements",
      "Treatment Improvements"
    )
  ) +
  labs(x = "Scenario Under High HCC Risk Assumption") +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)
  ) +
  geom_point(
    data = new_scenarios_high_risk,
    aes(Scenario, `Life Years per Cohort of 100`, shape = type),
    size = 3,
    show.legend = TRUE
  ) +
  scale_y_continuous(limits = c(0, 60), expand = expansion(0))

new_scens_plot_high_risk

ggsave(
  plot = new_scens_plot_high_risk,
  here("figures", "scenarios_new_high_risk.jpg"),
  width = 6,
  height = 4
)
ggsave(
  plot = new_scens_plot_high_risk,
  here("figures", "scenarios_new_high_risk.png"),
  width = 6,
  height = 4,
  dpi = 600
)
ggsave(
  plot = new_scens_plot_high_risk,
  here("figures", "scenarios_new_high_risk.pdf"),
  width = 6,
  height = 4,
  dpi = 600
)
