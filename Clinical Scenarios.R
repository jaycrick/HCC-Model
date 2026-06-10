# Init. ----------------------------
source(here("Sensitivity_analyses.R"))

library(conflicted)
library(tidyverse)
library(viridis)
library(here)
library(scales)

# Baseline no USs at all

noUS_scen <- life_add(p_NA = 1, p_CA = 0) # 14.19459

Baseline <- life_add() # 66.4987

# Scenario 1 ------ Evidence based Adherence Increase

scen_1a <- life_add(
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

# Scenario 2 ------ Gold Star Adherence Increase

scen_1b <- life_add(
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

# Scenario 3 ------ Improved US

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

# Scenario 4 ------ Better Treatments

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

# Graph all scens -------

resampled_pars <- readRDS(file = "resampled_pars_fin.rds")

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
  as.double(quantile(resampled_pars$Baseline, 0.975)) # ,
  # "No USs", noUS_scen, "none"
)

all_scens_plot <- ggplot(
  data = scenarios,
  aes(
    fct_inorder(Scenario),
    `Life Years per Cohort of 100`,
    # fill = fct_inorder(Scenario)
  )
) +
  geom_col() +
  # geom_errorbar(aes(ymin = lowCI, ymax = highCI), width = 0.5, size = 1) +
  geom_violin(
    data = violin_scenarios,
    mapping = aes(Scenario, `Life Years per Cohort of 100`)
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_fill_viridis(discrete = T) +
  labs(x = "Scenario") +
  theme(legend.position = "none")

all_scens_plot

ggsave(
  plot = all_scens_plot,
  here("figures", "scenarios_all.jpg"),
  height = 4,
  width = 6
)

# Graph of just new scenarios ----

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


param_names <- readRDS(file = "param_names.rds")
param_names[c(
  "treat_scens",
  "diag_scen",
  "opt_adh",
  "real_adh"
)] <- new_scenarios[[2]]
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
saveRDS(param_names, file = "param_names.rds")


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

ggplot(violin_scenarios, aes(Scenario, `Life Years per Cohort of 100`)) +
  geom_violin()

new_scens_plot <- ggplot(
  data = new_violin_scenarios,
  mapping = aes(Scenario, `Life Years per Cohort of 100`)
) +
  # geom_col() +
  geom_violin(trim = T) +
  # geom_errorbar(aes(ymin = lowCI, ymax = highCI), width = 0.5, size = 1) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_x_discrete(
    limits = c(
      "Realistic Adherence\nImprovements",
      "Optimistic Adherence\nImprovements",
      "Diagnosis Sensitivity\nImprovements",
      "Treatment Improvements"
    )
  ) +
  # scale_fill_viridis(discrete = T) +
  labs(x = "Scenario") +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)
  ) +
  # stat_summary(fun.data = "mean_sdl", mult = 1, geom = "pointrange", color = "red") +
  geom_point(
    data = new_scenarios,
    aes(
      Scenario,
      `Life Years per Cohort of 100`,
      shape = type
    ),
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


# Testing diminishing returns ----

## Adherence Tests
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
    all = T
  ) %>%
  merge(
    tibble(`Percent change in Category` = 0:30 / 100, `Non-Adherence`),
    all = T
  ) %>%
  merge(
    tibble(`Percent change in Category` = 0:40 / 100, `Complete Adherence`),
    all = T
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
  geom_line(size = 1.5) +
  theme_bw(base_family = "serif", base_size = 12) +
  ggtitle("Effect of Adherence Increase\non Intervention Outcomes") +
  scale_x_continuous(
    labels = label_percent(),
    expand = expansion(0)
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  )

plot_adh_inc

ggsave(
  here("Figures", "increasing_Adh_plot_10_year.jpg"),
  plot_adh_inc,
  width = 6,
  height = 3
)


## US tests ----

`BCLC 0` <- map_dbl(0:50 / 100, ~ life_add(Sens_0 = Sens_0_def + .x))
`BCLC A` <- map_dbl(0:30 / 100, ~ life_add(Sens_A = Sens_A_def + .x))
`BCLC B` <- map_dbl(0:15 / 100, ~ life_add(Sens_B = Sens_B_def + .x))

US_inc_tbl <- tibble(`Increase in Sensitivity` = 0:50 / 100) %>%
  merge(tibble(`Increase in Sensitivity` = 0:50 / 100, `BCLC 0`), all = T) %>%
  merge(tibble(`Increase in Sensitivity` = 0:30 / 100, `BCLC A`), all = T) %>%
  merge(tibble(`Increase in Sensitivity` = 0:15 / 100, `BCLC B`), all = T) %>%
  pivot_longer(
    !`Increase in Sensitivity`,
    names_to = "US Stage",
    values_to = "Life Extended"
  )

plot_US_inc <- ggplot(
  US_inc_tbl,
  aes(`Increase in Sensitivity`, `Life Extended`, color = `US Stage`)
) +
  geom_line(size = 1.5) +
  theme_bw(base_family = "serif", base_size = 12) +
  ggtitle("Effect of US Sensitivities on Intervention Outcomes") +
  scale_x_continuous(
    labels = label_percent(),
    expand = expansion(0)
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  )

plot_US_inc

ggsave(
  here("Figures", "increasing_US_Sens_plot_10_year.jpg"),
  plot_US_inc,
  width = 6,
  height = 4
)


# Scenarios with changed HCC risk -----------------------------------------

Baseline_low_risk <- life_add(p_HCC = p_HCC_def / 2)
Baseline_high_risk <- life_add(p_HCC = p_HCC_def * 2)

## Scenario 1 ------ Evidence based Adherence Increase

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

# Scenario 2 ------ Gold Star Adherence Increase

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

# Scenario 3 ------ Improved US

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

# Scenario 4 ------ Better Treatments

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

## Graph of low risk scenarios ----

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


param_names <- readRDS(file = "param_names.rds")
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
saveRDS(param_names, file = "param_names.rds")


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
  geom_violin(trim = T) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_x_discrete(
    limits = c(
      "Realistic Adherence\nImprovements",
      "Optimistic Adherence\nImprovements",
      "Diagnosis Sensitivity\nImprovements",
      "Treatment Improvements"
    )
  ) +
  # scale_fill_viridis(discrete = T) +
  labs(x = "Scenario Under Low HCC Risk Assumption") +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)
  ) +
  geom_point(
    data = new_scenarios_low_risk,
    aes(
      Scenario,
      `Life Years per Cohort of 100`,
      shape = type
    ),
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


## Graph of high risk scenarios ----

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


param_names <- readRDS(file = "param_names.rds")
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
saveRDS(param_names, file = "param_names.rds")


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
  geom_violin(trim = T) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_x_discrete(
    limits = c(
      "Realistic Adherence\nImprovements",
      "Optimistic Adherence\nImprovements",
      "Diagnosis Sensitivity\nImprovements",
      "Treatment Improvements"
    )
  ) +
  # scale_fill_viridis(discrete = T) +
  labs(x = "Scenario Under High HCC Risk Assumption") +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)
  ) +
  geom_point(
    data = new_scenarios_high_risk,
    aes(
      Scenario,
      `Life Years per Cohort of 100`,
      shape = type
    ),
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
