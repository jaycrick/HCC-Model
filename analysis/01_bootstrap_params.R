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
#   param_names.rds        — named list of key reported values for the manuscript
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

library(hypogeom) # Written for this project by D. Pavlyshyn; see README for install
library(parallel)
library(DescTools)
library(corrplot)
library(patchwork)
library(cowplot)

dir.create(here("figures"), showWarnings = FALSE)

set.seed(314) # bootstrap-specific seed (different from model seed)

# ----- Observed data ----------------------------------------------------------
med_0 <- 38 # n = 12
med_a <- 25 # n = 101
med_b <- 10 # n = 139
med_c <- 7 # n = 138
med_d <- 4 # n = 210

meds <- c(med_0, med_a, med_b, med_c, med_d)
init_hypo_probs <- c(0.08, 0.06, 0.34, 0.35, 0.16, 0.22, 0.54, 0.01)

hypo_probs <- c(
  0.07982013,
  0.05540597,
  0.33922926,
  0.35681215,
  0.16023503,
  0.21878495,
  0.54017653,
  0.01015650
)

# ----- Point estimate optimisation --------------------------------------------
optim(
  init_hypo_probs,
  function(x) {
    neg_log_lik(x, meds = meds[1:5], Dx_0 = 36, Dx_A = 51, Dx_sym_early = 52)
  },
  method = "Nelder-Mead",
  control = list(maxit = 10000)
)
optim(
  init_hypo_probs,
  function(x) {
    neg_log_lik(x, meds = meds[1:5], Dx_0 = 36, Dx_A = 51)
  },
  method = "BFGS",
  control = list(maxit = 10000)
)

set.seed(314)

# ----- Load pre-computed bootstrap outputs ------------------------------------
# The bootstrap sampling + parallel optimisation was run on an HPC cluster
# (see the commented block below for the sampling/optimisation code).
# Outputs are loaded from data-derived/.

milton_list_of_resampled <- readRDS(
  here("data-derived", "list_of_resampled_milton_output.rds")
)

bootstrapped_data <- tibble(all = milton_list_of_resampled) %>%
  unnest_wider(
    all,
    names_sep = "Yes",
    names_repair = ~ c(
      "m_0",
      "m_A",
      "m_B",
      "m_C",
      "m_D",
      "D_0",
      "D_A",
      "D_alpha"
    )
  ) %>%
  suppressMessages()

bootstrapped_data_long <- bootstrapped_data %>%
  pivot_longer(c(
    "m_0",
    "m_A",
    "m_B",
    "m_C",
    "m_D",
    "D_0",
    "D_A",
    "D_alpha"
  )) %>%
  mutate(name = map_chr(name, ~ paste(Boot_labels[.x])))

ggplot(bootstrapped_data_long, aes(x = value)) +
  geom_bar() +
  facet_wrap(
    vars(name),
    scales = "free",
    labeller = label_parsed,
    nrow = 3
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_y_continuous(expand = expansion(c(0, 0.05)))

ggsave(here("figures", "Bootstrap_sample_ranges.jpg"), width = 6, height = 4)
ggsave(
  here("figures", "Bootstrap_sample_ranges.png"),
  width = 6,
  height = 4,
  dpi = 600
)
ggsave(
  here("figures", "Bootstrap_sample_ranges.pdf"),
  width = 6,
  height = 4,
  dpi = 600
)

# ----- Process MLE outputs from cluster job -----------------------------------
list_of_pars <- readRDS(here("data-derived", "miltonoutput.rds"))

resampled_pars <- tibble(list_of_pars) %>%
  unnest_wider(list_of_pars) %>%
  select(par) %>%
  unnest_wider(
    par,
    names_sep = "k",
    names_repair = ~ c(
      "p_0",
      "p_A",
      "p_B",
      "p_C",
      "p_D",
      "Sens_0",
      "Sens_A",
      "Sym_A"
    )
  ) %>%
  suppressMessages()

saveRDS(resampled_pars, file = here("data-derived", "resampled_pars.rds"))

resampled_pars <- readRDS(here("data-derived", "resampled_pars.rds"))

table(resampled_pars$p_0 < resampled_pars$p_A)
table(resampled_pars$p_0 < resampled_pars$p_B)
table(resampled_pars$p_0 < resampled_pars$p_C)
table(resampled_pars$p_0 < resampled_pars$p_D)

table(resampled_pars$p_A < resampled_pars$p_0)
table(resampled_pars$p_A < resampled_pars$p_B)
table(resampled_pars$p_A < resampled_pars$p_C)
table(resampled_pars$p_A < resampled_pars$p_D)

# Initialise param_names ledger (or update if it already exists)
if (file.exists(here("data-derived", "param_names.rds"))) {
  param_names <- readRDS(here("data-derived", "param_names.rds"))
} else {
  param_names <- list()
}
param_names[c("muA_gt_mu0")] <- table(
  resampled_pars$p_A < resampled_pars$p_0
)[["TRUE"]] /
  10000
saveRDS(param_names, file = here("data-derived", "param_names.rds"))

# ----- Rates table for plotting -----------------------------------------------
variable_labels <- c(
  "p_0" = "BCLC Stage\n0 to A HCC",
  "p_A" = "BCLC Stage\nA to B HCC",
  "p_B" = "BCLC Stage\nB to C HCC",
  "p_C" = "BCLC Stage\nC to D HCC",
  "p_D" = "BCLC Stage\nD HCC to Death"
)

rates_tbl <- tibble(
  rate = factor(c("p_0", "p_A", "p_B", "p_C", "p_D")),
  Point = hypo_probs[1:5],
  Lower = c(
    quantile(resampled_pars$`p_0`, 0.025),
    quantile(resampled_pars$p_A, 0.025),
    quantile(resampled_pars$p_B, 0.025),
    quantile(resampled_pars$p_C, 0.025),
    quantile(resampled_pars$p_D, 0.025)
  ),
  Upper = c(
    quantile(resampled_pars$`p_0`, 0.975),
    quantile(resampled_pars$p_A, 0.975),
    quantile(resampled_pars$p_B, 0.975),
    quantile(resampled_pars$p_C, 0.975),
    quantile(resampled_pars$p_D, 0.975)
  )
)

ggplot(
  resampled_pars[1:5] %>% pivot_longer(everything()),
  aes(fct_inorder(name), value)
) +
  geom_boxplot(outlier.shape = 1, outlier.alpha = 0.2) +
  geom_point(
    mapping = aes(x = rate, y = Point),
    data = rates_tbl,
    inherit.aes = FALSE,
    colour = "red"
  ) +
  scale_x_discrete(labels = variable_labels) +
  labs(x = NULL, y = "Probability of\nProgressing per Month") +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0))

ggsave(here("figures", "HCC_Rates_Range_Plot.jpg"), width = 6, height = 3)

resampled_rates_long <- resampled_pars[1:5] %>%
  pivot_longer(everything(), names_to = "rate") %>%
  mutate(rate = map_chr(rate, ~ paste(Boot_labels[.x])))

ggplot(resampled_rates_long, aes(value)) +
  geom_histogram(bins = 20) +
  geom_vline(
    mapping = aes(xintercept = Point),
    data = rates_tbl %>%
      mutate(rate = map_chr(rate, ~ paste(Boot_labels[as.character(.x)]))),
    colour = "red"
  ) +
  facet_wrap(vars(rate), nrow = 3, scales = "free", labeller = label_parsed) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_y_continuous(expand = expansion(c(0, 0.05)))

ggsave(here("figures", "HCC_Rates_Hist.jpg"), width = 6, height = 7)
ggsave(here("figures", "HCC_Rates_Hist.png"), dpi = 600, width = 6, height = 7)
ggsave(here("figures", "HCC_Rates_Hist.pdf"), dpi = 600, width = 6, height = 7)

ggplot(
  resampled_pars[1:5] %>% pivot_longer(everything()),
  aes(fct_inorder(name), value)
) +
  geom_violin(scale = "width") +
  geom_point(
    mapping = aes(x = rate, y = Point),
    data = rates_tbl,
    inherit.aes = FALSE,
    colour = "red"
  ) +
  scale_x_discrete(labels = Boot_labels) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0))

ggsave(here("figures", "HCC_Rates_Violin.jpg"), width = 6, height = 4)
ggsave(
  here("figures", "HCC_Rates_Violin.png"),
  dpi = 600,
  width = 6,
  height = 4
)
ggsave(
  here("figures", "HCC_Rates_Violin.pdf"),
  dpi = 600,
  width = 6,
  height = 4
)

mean_years_tbl <- resampled_pars[1:5] %>%
  pivot_longer(everything()) %>%
  mutate(mean_years = 1 / value / 12, name = fct_inorder(name))

HCC_Time_Range_Plot <- ggplot(
  mean_years_tbl %>%
    group_by(name) %>%
    summarise(
      mean025 = quantile(mean_years, 0.025),
      mean975 = quantile(mean_years, 0.975)
    ),
  aes(name, ymin = mean025, ymax = mean975)
) +
  geom_errorbar() +
  geom_point(
    mapping = aes(x = rate, y = 1 / Point / 12),
    data = rates_tbl,
    inherit.aes = FALSE
  ) +
  scale_x_discrete(
    labels = c(
      p_0 = "BCLC Stage 0",
      p_A = "BCLC Stage A",
      p_B = "BCLC Stage B",
      p_C = "BCLC Stage C",
      p_D = "BCLC Stage D"
    )
  ) +
  labs(x = NULL, y = "Mean Number of\nYears In BCLC Stage") +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_y_continuous(limits = c(0, 3), expand = expansion(0))

HCC_Time_Range_Plot

param_names <- readRDS(here("data-derived", "param_names.rds"))
param_names[c("mean_0", "mean_A", "mean_B", "mean_C", "mean_D")] <-
  1 / rates_tbl$Point / 12
param_names[c(
  "mean_0_Lower",
  "mean_A_Lower",
  "mean_B_Lower",
  "mean_C_Lower",
  "mean_D_Lower"
)] <- 1 / rates_tbl$Upper / 12
param_names[c(
  "mean_0_Upper",
  "mean_A_Upper",
  "mean_B_Upper",
  "mean_C_Upper",
  "mean_D_Upper"
)] <- 1 / rates_tbl$Lower / 12
saveRDS(param_names, file = here("data-derived", "param_names.rds"))

ggsave(
  here("figures", "HCC_Time_Range_Plot.jpg"),
  HCC_Time_Range_Plot,
  width = 6,
  height = 4
)
ggsave(
  here("figures", "HCC_Time_Range_Plot.png"),
  HCC_Time_Range_Plot,
  width = 6,
  height = 4,
  dpi = 600
)
ggsave(
  here("figures", "HCC_Time_Range_Plot.pdf"),
  HCC_Time_Range_Plot,
  width = 6,
  height = 4,
  dpi = 600
)

set.seed(315)

US_sensitivity_Cal <- ggplot(resampled_pars, aes(x = Sens_0, y = Sens_A)) +
  geom_area(
    aes(x = seq(0, 1, length.out = 10000), y = seq(0, 1, length.out = 10000)),
    alpha = 0.25
  ) +
  geom_hex(alpha = 0.75, bins = c(100, 80)) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_fill_viridis(
    option = "C",
    trans = "log",
    breaks = c(1, 2, 5, 10, 20, 35)
  ) +
  scale_x_continuous(
    labels = label_percent(),
    expand = expansion(0),
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    labels = label_percent(),
    expand = expansion(0),
    limits = c(0, 1)
  ) +
  geom_point(
    aes(x = params_def$Sens_0, y = params_def$Sens_A),
    shape = 4,
    stroke = 1,
    size = 1,
    colour = "black"
  ) +
  labs(
    x = "BCLC Stage 0 HCC\nUltrasound Sensitivity",
    y = "BCLC Stage A HCC\nUltrasound Sensitivity",
    fill = "Number of\nSamples"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  )

US_sensitivity_Cal

ggsave(
  here("figures", "US_sensitivity_Cal.jpg"),
  US_sensitivity_Cal,
  height = 3.5,
  width = 6
)
ggsave(
  here("figures", "US_sensitivity_Cal.png"),
  US_sensitivity_Cal,
  height = 3.5,
  width = 6,
  dpi = 600
)
ggsave(
  here("figures", "US_sensitivity_Cal.pdf"),
  US_sensitivity_Cal,
  height = 3.5,
  width = 6,
  dpi = 600
)

summary(resampled_pars$Sens_0)
quantile(resampled_pars$Sens_0, c(0.025, 0.975))

summary(resampled_pars$Sens_A)
quantile(resampled_pars$Sens_A, c(0.025, 0.975))

# ----- Scenario life-year calculations ----------------------------------------
# Scenario definitions, risk levels and legacy column names: R/scenarios.R.
resampled_pars <- bind_cols(
  resampled_pars,
  scenario_columns(run_scenarios(resampled_pars))
)

head(resampled_pars)

summary(resampled_pars$scen_1a_diff)
quantile(resampled_pars$scen_3_diff, c(0.025, 0.975))

pairs(
  head(resampled_pars, 1000) %>%
    select(c("scen_1a_diff", "scen_2_diff", "scen_3_diff"))
)
table(resampled_pars$scen_1a_diff > resampled_pars$scen_3_diff)

pairs(
  head(resampled_pars, 1000) %>% select(c("p_0", "p_A", "p_B", "p_C", "p_D"))
)

param_names <- readRDS(here("data-derived", "param_names.rds"))
param_names[c(
  "opt_vs_diag",
  "opt_vs_treat",
  "opt_vs_real",
  "real_vs_diag",
  "real_vs_treat",
  "treat_vs_diag"
)] <- c(
  table(resampled_pars$scen_1b_diff > resampled_pars$scen_2_diff)[["TRUE"]] /
    10000,
  table(resampled_pars$scen_1b_diff > resampled_pars$scen_3_diff)[["TRUE"]] /
    10000,
  table(resampled_pars$scen_1b_diff > resampled_pars$scen_1a_diff)[["TRUE"]] /
    10000,
  table(resampled_pars$scen_1a_diff > resampled_pars$scen_2_diff)[["TRUE"]] /
    10000,
  table(resampled_pars$scen_1a_diff > resampled_pars$scen_3_diff)[["TRUE"]] /
    10000,
  table(resampled_pars$scen_3_diff > resampled_pars$scen_2_diff)[["TRUE"]] /
    10000
)
saveRDS(param_names, file = here("data-derived", "param_names.rds"))

saveRDS(resampled_pars, file = here("data-derived", "resampled_pars_fin.rds"))

resampled_pars <- readRDS(here("data-derived", "resampled_pars_fin.rds"))

var_names <- colnames(resampled_pars)[c(1:9)]

# ----- Uncertainty analysis ---------------------------------------------------
ggplot(resampled_pars, aes(x = Baseline, y = "")) +
  geom_violin() +
  theme_bw(base_family = "serif", base_size = 12)

ggplot(resampled_pars, aes(x = Baseline)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.5) +
  geom_vline(xintercept = life_add(), colour = "red") +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_y_continuous(expand = expansion(0), limits = c(0, 0.08)) +
  geom_function(
    fun = dnorm,
    args = list(
      mean = mean(resampled_pars$Baseline),
      sd = sd(resampled_pars$Baseline)
    ),
    colour = "blue",
    linewidth = 1
  ) +
  labs(x = "Life Years Gained", y = "Density")

ggsave(here("figures", "bootstrap_var_base.jpg"), width = 6, height = 3)
ggsave(
  here("figures", "bootstrap_var_base.png"),
  width = 6,
  height = 3,
  dpi = 600
)
ggsave(
  here("figures", "bootstrap_var_base.pdf"),
  width = 6,
  height = 3,
  dpi = 600
)

param_names <- readRDS(here("data-derived", "param_names.rds"))
param_names[c("sd_noint", "mean_noint")] <- c(
  sd(resampled_pars$Baseline),
  mean(resampled_pars$Baseline)
)
saveRDS(param_names, file = here("data-derived", "param_names.rds"))

ggplot(resampled_pars, aes(x = Baseline)) +
  geom_boxplot() +
  theme_bw(base_family = "serif", base_size = 12)

# ----- Sensitivity / uncertainty analysis with partial correlation ------------
library(ppcor)

pars_wide <- resampled_pars %>%
  dplyr::select(1:8, Baseline)

pars_logit <- resampled_pars %>%
  dplyr::select(1:8, Baseline) %>%
  mutate(across(-Baseline, Logit))

pairs(head(pars_logit, 1000))
pairs(head(pars_wide, 1000))

lm_fit <- lm(
  pars_wide$Baseline ~ .,
  data = pars_wide %>% dplyr::select(-Baseline)
)
lm_logit <- lm(
  pars_logit$Baseline ~ .,
  data = pars_logit %>% dplyr::select(-Baseline)
)

summary.lm(lm_fit)$r.squared
summary.lm(lm_logit)$r.squared

cor_matrix <- cor(pars_wide)
cor_matrix_logit <- cor(pars_logit)
pcor_matrix_logit <- pcor(pars_logit)$estimate

corrplot(
  cor_matrix_logit,
  method = "color",
  type = "lower",
  diag = TRUE,
  addCoef.col = "black",
  title = "Logit Pearson correlation coefficient"
)

corrplot(
  pcor_matrix_logit,
  method = "color",
  type = "lower",
  diag = TRUE,
  addCoef.col = "black",
  title = "Logit Partial correlation coefficient"
)

long_pcor_logit <- pcor_matrix_logit %>%
  as_tibble() %>%
  mutate("row" = fct_inorder(var_names)) %>%
  pivot_longer(-row, names_to = "col") %>%
  mutate(col = fct_inorder(col))

plot_par_cor_logit <- ggplot(
  long_pcor_logit,
  aes(x = col, y = row, fill = value)
) +
  geom_tile() +
  scale_fill_gradient2(limit = c(-1, 1), low = "red", high = "#075AFF") +
  scale_y_discrete(expand = expansion(0), limits = rev) +
  scale_x_discrete(expand = expansion(0), position = "top") +
  theme_bw(base_family = "serif", base_size = 12) +
  geom_text(aes(label = round(value, 2)), color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0)) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(title = "Partial Rank\nCoefficient"))

plot_par_cor_logit

ggsave(
  here("figures", "plot_par_cor_logit.jpg"),
  plot = plot_par_cor_logit,
  width = 6,
  height = 5
)

pcor_base_logit <- long_pcor_logit %>%
  dplyr::filter(row == "Baseline" & col != "Baseline") %>%
  rename(Variable = col)

plot_cor_base_logit <- ggplot(
  pcor_base_logit,
  aes(x = Variable, y = value, fill = value)
) +
  geom_col(colour = "black") +
  scale_fill_gradient2(limit = c(-1, 1), low = "red", high = "#075AFF") +
  theme_bw(base_family = "serif", base_size = 12) +
  labs(y = "Partial Correlation Coefficient", x = NULL) +
  guides(fill = guide_colourbar(title = NULL)) +
  scale_y_continuous(limits = c(-1, 1), expand = expansion(0)) +
  scale_x_discrete(labels = TeX_labels)

plot_cor_base_logit

ggsave(
  here("figures", "plot_cor_base_logit.jpg"),
  plot = plot_cor_base_logit,
  width = 6,
  height = 3
)
ggsave(
  here("figures", "plot_cor_base_logit.png"),
  plot = plot_cor_base_logit,
  width = 6,
  height = 3,
  dpi = 600
)
ggsave(
  here("figures", "plot_cor_base_logit.pdf"),
  plot = plot_cor_base_logit,
  width = 6,
  height = 3,
  dpi = 600
)

detach("package:ppcor", unload = TRUE)
detach("package:MASS", unload = TRUE)

param_names <- readRDS(here("data-derived", "param_names.rds"))
param_names[c("pcor_logit_0", "pcor_logit_US_0", "pcor_logit_US_A")] <- c(
  pcor_matrix_logit["Baseline", "p_0"],
  pcor_matrix_logit["Baseline", "Sens_0"],
  pcor_matrix_logit["Baseline", "Sens_A"]
)
saveRDS(param_names, file = here("data-derived", "param_names.rds"))

# ----- Results summary table --------------------------------------------------
default_param_tbl <- tibble(params = list(hypo_probs)) %>%
  unnest_wider(
    params,
    names_sep = "t",
    names_repair = ~ c(
      "p_0",
      "p_A",
      "p_B",
      "p_C",
      "p_D",
      "Sens_0",
      "Sens_A",
      "Sym_A"
    )
  )

results_tbl <- resampled_pars %>%
  pivot_longer(everything(), names_to = "result") %>%
  group_by(result) %>%
  summarise(
    point = median(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975)
  )

# Point estimates at the hypo_probs parameters (Sym_B derived from Sym_A, as
# for the bootstrap draws), base risk level only.
point_wide <- scenario_columns(
  run_scenarios(default_param_tbl, risk = risk_multipliers["base"])
)
point_results <- c(
  "scen_1a_diff",
  "scen_1b_diff",
  "scen_2_diff",
  "scen_3_diff",
  "scen_1a",
  "scen_1b",
  "scen_2",
  "scen_3",
  "scen_1a_3",
  "scen_1a_4",
  "scen_1b_3",
  "scen_1b_4"
)
results_tbl <- results_tbl %>%
  rows_update(
    tibble(
      result = point_results,
      point = unlist(point_wide[1, point_results], use.names = FALSE)
    ),
    by = "result"
  )

# ----- Multipanel calibration figure ------------------------------------------
# Assembles bootstrap calibration plots with the univariate sensitivity plot
# (recreated here from model_core functions).

# Recreate plot_univar_mod (25%-variation univariate sensitivity)
param_tbl_25 <- life_tibble(here("data", "Parameters_25.csv"))
param_chang_tbl <- param_tbl_25 %>%
  dplyr::filter(
    years_mid_comp != 0,
    !Parameter %in%
      c(
        "p_HCC",
        "p_0",
        "p_A",
        "p_B",
        "p_C",
        "p_D",
        "Sym_0",
        "Sym_A",
        "Sym_B",
        "Sym_C",
        "Sym_D"
      )
  ) %>%
  mutate(
    Parameter = fct_reorder(.f = Parameter, .x = years_mid_comp, abssum),
    Type = factor(
      word(Parameter, sep = fixed("_")),
      levels = c("Pr", "Sens", "L", "p", "Sym")
    )
  ) %>%
  dplyr::select(Parameter, Bound, years_mid_comp, Type)

facet_labels <- c(
  L = "Treatment Quality",
  Pr = "Adherence\nto Screening",
  p = "Monthly Risk",
  Sym = "Monthly probability\nof symptoms",
  Sens = "Ultrasound\nSensitivity"
)

plot_univar_mod <- ggplot(
  param_chang_tbl,
  aes(y = Parameter, x = years_mid_comp)
) +
  geom_col(position = "identity", colour = "black") +
  labs(x = "Additional Life Years Gained", y = NULL, title = NULL) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_y_discrete(labels = TeX_labels) +
  facet_wrap(
    ~Type,
    ncol = 1,
    scales = "free_y",
    labeller = labeller(Type = facet_labels)
  )

# Variant plots with rotated labels for multipanel layout
HCC_Time_Range_Plot_side_labs <- ggplot(
  mean_years_tbl %>%
    group_by(name) %>%
    summarise(
      mean025 = quantile(mean_years, 0.025),
      mean975 = quantile(mean_years, 0.975)
    ),
  aes(name, ymin = mean025, ymax = mean975)
) +
  geom_errorbar() +
  geom_point(
    mapping = aes(x = rate, y = 1 / Point / 12),
    data = rates_tbl,
    inherit.aes = FALSE
  ) +
  scale_x_discrete(
    labels = c(
      p_0 = "BCLC Stage 0",
      p_A = "BCLC Stage A",
      p_B = "BCLC Stage B",
      p_C = "BCLC Stage C",
      p_D = "BCLC Stage D"
    )
  ) +
  labs(x = NULL, y = "Mean Number of\nYears In BCLC Stage") +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_y_continuous(limits = c(0, 3), expand = expansion(0)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

plot_cor_base_logit_side_labs <- ggplot(
  pcor_base_logit,
  aes(x = Variable, y = value, fill = value)
) +
  geom_col(colour = "black") +
  scale_fill_gradient2(limit = c(-1, 1), low = "red", high = "#075AFF") +
  theme_bw(base_family = "serif", base_size = 12) +
  labs(y = "Partial Correlation Coefficient", x = NULL) +
  guides(fill = guide_colourbar(title = NULL)) +
  scale_y_continuous(limits = c(-1, 1), expand = expansion(0)) +
  scale_x_discrete(labels = TeX_labels) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "bottom"
  )

US_sensitivity_Cal_side_labs <- ggplot(
  resampled_pars,
  aes(x = Sens_0, y = Sens_A)
) +
  geom_area(
    aes(x = seq(0, 1, length.out = 10000), y = seq(0, 1, length.out = 10000)),
    alpha = 0.25
  ) +
  geom_hex(alpha = 0.75, bins = c(100, 80)) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_fill_viridis(
    option = "C",
    trans = "log",
    breaks = c(1, 2, 5, 10, 20, 35)
  ) +
  scale_x_continuous(
    labels = label_percent(),
    expand = expansion(0),
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    labels = label_percent(),
    expand = expansion(0),
    limits = c(0, 1)
  ) +
  geom_point(
    aes(x = params_def$Sens_0, y = params_def$Sens_A),
    shape = 4,
    stroke = 1,
    size = 1,
    colour = "black"
  ) +
  labs(
    x = "BCLC Stage 0 HCC\nUltrasound Sensitivity",
    y = "BCLC Stage A HCC\nUltrasound Sensitivity",
    fill = "Number of\nSamples"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

multipanel_cals <- plot_grid(
  HCC_Time_Range_Plot_side_labs,
  US_sensitivity_Cal_side_labs,
  plot_cor_base_logit_side_labs,
  plot_univar_mod,
  labels = "AUTO",
  rel_widths = c(1, 2)
)

multipanel_cals

ggsave(here("figures", "multipanel_cals.png"), width = 6, height = 8, dpi = 600)
ggsave(here("figures", "multipanel_cals.pdf"), width = 6, height = 8, dpi = 600)
