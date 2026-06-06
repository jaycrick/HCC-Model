# install.packages("pacman")

library(cli)
# devtools::install_github("damian-t-p/hypogeom")

library(hypogeom)
# library(ppcor)
library(expm)
library(here)
library(viridis)
library(scales)
library(parallel)
library(DescTools)
library(corrplot)
library(latex2exp)
library(tidyverse)

load("Sensitivity_analyses.RData")

set.seed(314)

med_0 <- 38 # n = 12
med_a <- 25 # n = 101
med_b <- 10 # n = 139
med_c <- 7 # n = 138
med_d <- 4 # n = 210


meds <- c(med_0, med_a, med_b, med_c, med_d)
init_hypo_probs <- c(0.08, 0.06, 0.34, 0.35, 0.16, 0.22, 0.54, 0.01)

hypo_probs <- c(0.07982013, 0.05540597, 0.33922926, 0.35681215, 0.16023503, 0.21878495, 0.54017653, 0.01015650)


neg_log_lik <- function(hypo_probs, meds = c(med_0, med_a, med_b, med_c, med_d), Dx_0 = 36, Dx_A = 51, Dx_sym_early = 52) {
  while (prod((0 <= hypo_probs) * (hypo_probs <= 1)) == 0) {
    return(NaN)
  }

  p_0a <- hypo_probs[1]
  p_Ab <- hypo_probs[2]
  p_Bc <- hypo_probs[3]
  US_00 <- hypo_probs[6]
  US_Aa <- hypo_probs[7]
  Sym_Aa <- hypo_probs[8]
  Sym_Bb <- 8 * Sym_Aa
  matrix_of_trans <- matrix(
    c(
      1, rep(0, 14),
      0, 1, rep(0, 13),
      rep(0, 3), 1 - `p_0a`, rep(0, 5), `p_0a`, rep(0, 5),
      rep(0, 4), 1 - `p_0a`, rep(0, 5), `p_0a`, rep(0, 4),
      rep(0, 5), 1 - `p_0a`, rep(0, 5), `p_0a`, rep(0, 3),
      rep(0, 6), 1 - `p_0a`, rep(0, 5), `p_0a`, rep(0, 2),
      rep(0, 7), 1 - `p_0a`, rep(0, 5), `p_0a`, rep(0, 1),
      (1 - `p_0a`) * US_00, `p_0a` * US_Aa, (1 - `p_0a`) * (1 - US_00), rep(0, 5), `p_0a` * (1 - US_Aa), rep(0, 6),
      rep(0, 9), 1 - p_Ab, rep(0, 4), p_Ab,
      rep(0, 10), 1 - p_Ab, rep(0, 3), p_Ab,
      rep(0, 11), 1 - p_Ab, rep(0, 2), p_Ab,
      rep(0, 12), 1 - p_Ab, 0, p_Ab,
      rep(0, 13), 1 - p_Ab, p_Ab,
      0, (1 - p_Ab) * US_Aa, rep(0, 6), (1 - p_Ab) * (1 - US_Aa), rep(0, 5), p_Ab,
      rep(0, 14), 1
    ),
    byrow = T,
    nrow = 15
  )

  end_vector <- c(US_00 / 6, 0, (1 - US_00) / 6, rep(1 / 6, 5), rep(0, 7)) %*% (matrix_of_trans %^% 100)
  end_vector[1]

  NLL <- -(dOrder(meds[5] - 1, distn = "hypogeom", probs = hypo_probs[5], mlen = 209, j = 210 / 2, log = T) +
    dOrder(meds[4] - 2, distn = "hypogeom", probs = hypo_probs[4:5], mlen = 137, j = 138 / 2, log = T) +
    dOrder(meds[3] - 3, distn = "hypogeom", probs = hypo_probs[3:5], mlen = 139, j = 140 / 2, log = T) +
    dOrder(meds[2] - 4, distn = "hypogeom", probs = hypo_probs[2:5], mlen = 101, j = 102 / 2, log = T) +
    dOrder(meds[1] - 5, distn = "hypogeom", probs = hypo_probs[1:5], mlen = 11, j = 12 / 2, log = T) +
    dbinom(Dx_0, 109,
      end_vector[1],
      log = T
    ) +
    dbinom(Dx_A, 109 - Dx_0,
      end_vector[2] / (1 - end_vector[1]),
      log = T
    ) +
    dbinom(Dx_sym_early, 158,
      Sym_Aa / (1 - (1 - p_Ab) * (1 - Sym_Aa)) + (1 - Sym_Aa / (1 - (1 - p_Ab) * (1 - Sym_Aa))) * Sym_Bb / (1 - (1 - p_Bc) * (1 - Sym_Bb)),
      log = T
    )
  )
  return(NLL)
}

optim(init_hypo_probs, function(x) {
  neg_log_lik(x, meds = meds[1:5], Dx_0 = 36, Dx_A = 51, Dx_sym_early = 52)
}, method = "Nelder-Mead", control = list(maxit = 10000))
optim(init_hypo_probs, function(x) {
  neg_log_lik(x, meds = meds[1:5], Dx_0 = 36, Dx_A = 51)
}, method = "BFGS", control = list(maxit = 10000))
# optim(init_hypo_probs, function(x) {
#   neg_log_lik(x, meds = meds[1:5], Dx_0 = 36, Dx_A = 51)
# }, method = "CG", control = list(maxit = 10000))

set.seed(314)

# list_of_resampled <- list()
#
# for (i in 1:10000) {
#   sam_d <- 1 + rhypogeom(209, hypo_probs[5]) # should be 210
#   sam_c <- 2 + rhypogeom(137, hypo_probs[4:5]) # should be 138
#   sam_b <- 3 + rhypogeom(139, hypo_probs[3:5])
#   sam_a <- 4 + rhypogeom(101, hypo_probs[2:5])
#   sam_0 <- 5 + rhypogeom(11, hypo_probs[1:5]) # should be 12
#
#   sam_meds <- c(median(sam_0), median(sam_a), median(sam_b), median(sam_c), median(sam_d))
#
#   Dx_0 <- rbinom(1, 109, 36 / 109)
#   Dx_A <- rbinom(1, 109 - Dx_0, 51 / 73)
#   Dx_sym_early <- rbinom(1, 158, 52 / 158)
#
#   list_of_resampled[[i]] <- c(sam_meds, Dx_0, Dx_A, Dx_sym_early)
# }


# n.cores <- detectCores()
# clust <- makeCluster(n.cores)
# clusterExport(clust, c("list_of_resampled", "init_hypo_probs", "neg_log_lik", "%^%", "dOrder", "dhypogeom", "phypogeom"))
#
# list_of_pars <- parLapply(clust, list_of_resampled, function(lor) {
#   optim(init_hypo_probs, function(x) {
#     neg_log_lik(x, meds = lor[1:5], Dx_0 = lor[6], Dx_A = lor[7], Dx_sym_early = lor[8])
#   }, method = "Nelder-Mead", control = list(maxit = 10000))
# })
# stopCluster(clust)

milton_list_of_resampled <- readRDS("list_of_resampled_milton_output.rds")

bootstrapped_data <- tibble(all = milton_list_of_resampled) %>%
  unnest_wider(all, names_sep = "Yes", names_repair = ~ c("m_0", "m_A", "m_B", "m_C", "m_D", "D_0", "D_A", "D_alpha")) %>%
  suppressMessages()

bootstrapped_data_long <- bootstrapped_data %>%
  pivot_longer(c("m_0", "m_A", "m_B", "m_C", "m_D", "D_0", "D_A", "D_alpha")) %>%
  mutate(name = map_chr(name, ~ paste(Boot_labels[.x])))

Tex_labeller <- function(value) {
  return(TeX_labels[value])
}

ggplot(bootstrapped_data_long, aes(x = value)) +
  geom_bar() +
  facet_wrap(
    vars(name),
    scales = "free", labeller = label_parsed, nrow = 3
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_y_continuous(expand = expansion(c(0, 0.05)))

ggsave(here("Figures", "Bootstrap_sample_ranges.jpg"), width = 6, height = 4)
ggsave(here("Figures", "Bootstrap_sample_ranges.png"), width = 6, height = 4, dpi = 600)
ggsave(here("Figures", "Bootstrap_sample_ranges.pdf"), width = 6, height = 4, dpi = 600)





list_of_pars <- readRDS("miltonoutput.rds")

resampled_pars <- tibble(list_of_pars) %>%
  unnest_wider(list_of_pars) %>%
  select(par) %>%
  unnest_wider(par, names_sep = "k", names_repair = ~ c("p_0", "p_A", "p_B", "p_C", "p_D", "Sens_0", "Sens_A", "Sym_A")) %>%
  suppressMessages()

saveRDS(resampled_pars, file = "resampled_pars.rds")

resampled_pars <- readRDS(file = "resampled_pars.rds")

table(resampled_pars$p_0 < resampled_pars$p_A)
table(resampled_pars$p_0 < resampled_pars$p_B)
table(resampled_pars$p_0 < resampled_pars$p_C)
table(resampled_pars$p_0 < resampled_pars$p_D)

table(resampled_pars$p_A < resampled_pars$p_0)
table(resampled_pars$p_A < resampled_pars$p_B)
table(resampled_pars$p_A < resampled_pars$p_C)
table(resampled_pars$p_A < resampled_pars$p_D)

param_names <- readRDS(file = "param_names.rds")
param_names[c("muA_gt_mu0")] <- table(resampled_pars$p_A < resampled_pars$p_0)["TRUE"] / 10000
saveRDS(param_names, file = "param_names.rds")

variable_labels <- c(
  "L_0" = "BCLC Stage 0\nHCC diagnosis",
  "L_A" = "BCLC Stage A\nHCC diagnosis",
  "L_B" = "BCLC Stage B\nHCC diagnosis",
  "L_C" = "BCLC Stage C\nHCC diagnosis",
  "L_D" = "BCLC Stage D\nHCC diagnosis",
  "p_CA" = "Completely adherent",
  "p_NA" = "Non-adherent",
  "p_Att" = "Attendence probability\nwithin partial adherents",
  "p_HCC" = "Developing HCC",
  "p_0" = "BCLC Stage\n0 to A HCC",
  "p_A" = "BCLC Stage\nA to B HCC",
  "p_B" = "BCLC Stage\nB to C HCC",
  "p_C" = "BCLC Stage\nC to D HCC",
  "p_D" = "BCLC Stage\nD HCC to Death",
  "Sym_0" = "Stage 0",
  "Sym_A" = "Stage A",
  "Sym_B" = "Stage B",
  "Sym_C" = "Stage C",
  "Sym_D" = "Stage D",
  "Sens_0" = "Stage 0",
  "Sens_A" = "Stage A",
  "Sens_B" = "Stage B",
  "Sens_C" = "Stage C"
)

rates_tbl <- tibble(
  rate = factor(c("p_0", "p_A", "p_B", "p_C", "p_D")), Point = hypo_probs[1:5],
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

ggplot(resampled_pars[1:5] %>% pivot_longer(everything()), aes(fct_inorder(name), value)) +
  geom_boxplot(outlier.shape = 1, outlier.alpha = 0.2) +
  geom_point(mapping = aes(x = rate, y = Point), data = rates_tbl, inherit.aes = F, colour = "red") +
  scale_x_discrete(labels = variable_labels) +
  labs(
    x = NULL,
    y = "Probability of\nProgressing per Month"
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0))

ggsave(here("Figures", "HCC_Rates_Range_Plot.jpg"), width = 6, height = 3)

resampled_rates_long <- resampled_pars[1:5] %>%
  pivot_longer(everything(), names_to = "rate") %>%
  mutate(rate = map_chr(rate, ~ paste(Boot_labels[.x])))

ggplot(resampled_rates_long, aes(value)) +
  geom_histogram(bins = 20) +
  geom_vline(mapping = aes(xintercept = Point), data = rates_tbl %>%
    mutate(rate = map_chr(rate, ~ paste(Boot_labels[as.character(.x)]))), colour = "red") +
  facet_wrap(vars(rate), nrow = 3, scales = "free", labeller = label_parsed) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_y_continuous(expand = expansion(c(0, 0.05)))

ggsave(here("Figures", "HCC_Rates_Hist.jpg"), width = 6, height = 7)
ggsave(here("Figures", "HCC_Rates_Hist.png"), dpi = 600, width = 6, height = 7)
ggsave(here("Figures", "HCC_Rates_Hist.pdf"), dpi = 600, width = 6, height = 7)

ggplot(resampled_pars[1:5] %>% pivot_longer(everything()), aes(fct_inorder(name), value)) +
  geom_violin(scale = "width") +
  geom_point(mapping = aes(x = rate, y = Point), data = rates_tbl, inherit.aes = F, colour = "red") +
  scale_x_discrete(labels = Boot_labels) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0))

ggsave(here("Figures", "HCC_Rates_Violin.jpg"), width = 6, height = 4)
ggsave(here("Figures", "HCC_Rates_Violin.png"), dpi = 600, width = 6, height = 4)
ggsave(here("Figures", "HCC_Rates_Violin.pdf"), dpi = 600, width = 6, height = 4)

mean_years_tbl <- resampled_pars[1:5] %>%
  pivot_longer(everything()) %>%
  mutate(
    mean_years = 1 / value / 12,
    name = fct_inorder(name)
  )

HCC_Time_Range_Plot <- ggplot(
  mean_years_tbl %>%
    group_by(name) %>%
    summarise(mean025 = quantile(mean_years, 0.025), mean975 = quantile(mean_years, 0.975)),
  aes(name, ymin = mean025, ymax = mean975)
) +
  geom_errorbar() +
  geom_point(mapping = aes(x = rate, y = 1 / Point / 12), data = rates_tbl, inherit.aes = F) +
  scale_x_discrete(labels = c(
    p_0 = "BCLC Stage 0",
    p_A = "BCLC Stage A",
    p_B = "BCLC Stage B",
    p_C = "BCLC Stage C",
    p_D = "BCLC Stage D"
  )) +
  labs(
    x = NULL,
    y = "Mean Number of\nYears In BCLC Stage"
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_y_continuous(limits = c(0, 3), expand = expansion(0))

HCC_Time_Range_Plot

param_names <- readRDS(file = "param_names.rds")
param_names[c("mean_0", "mean_A", "mean_B", "mean_C", "mean_D")] <- 1 / rates_tbl$Point / 12
param_names[c("mean_0_Lower", "mean_A_Lower", "mean_B_Lower", "mean_C_Lower", "mean_D_Lower")] <- 1 / rates_tbl$Upper / 12
param_names[c("mean_0_Upper", "mean_A_Upper", "mean_B_Upper", "mean_C_Upper", "mean_D_Upper")] <- 1 / rates_tbl$Lower / 12
saveRDS(param_names, file = "param_names.rds")

ggsave(here("Figures", "HCC_Time_Range_Plot.jpg"), HCC_Time_Range_Plot, width = 6, height = 4)
ggsave(here("Figures", "HCC_Time_Range_Plot.png"), HCC_Time_Range_Plot, width = 6, height = 4, dpi = 600)
ggsave(here("Figures", "HCC_Time_Range_Plot.pdf"), HCC_Time_Range_Plot, width = 6, height = 4, dpi = 600)

set.seed(315)

US_sensitivity_Cal <- ggplot(resampled_pars, aes(x = Sens_0, y = Sens_A)) +
  # geom_point() +
  # geom_density2d_filled(alpha = 0.75) +
  # geom_raster(stat = "density") +
  geom_area(aes(x = seq(0, 1, length.out = 10000), y = seq(0, 1, length.out = 10000)), alpha = 0.25) +
  geom_hex(alpha = 0.75, bins = c(100, 80)) +
  # geom_density2d() +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_fill_viridis(
    # discrete = T,
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
  geom_point(aes(x = Sens_0_def, y = Sens_A_def), shape = 4, stroke = 1, size = 1, colour = "black") +
  labs(
    # title = "Ultrasound Sensitivities Calibration",
    x = "BCLC Stage 0 HCC\nUltrasound Sensitivity",
    y = "BCLC Stage A HCC\nUltrasound Sensitivity",
    fill = "Number of\nSamples"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  )

US_sensitivity_Cal

ggsave(filename = here("Figures", "US_sensitivity_Cal.jpg"), US_sensitivity_Cal, height = 3.5, width = 6)
ggsave(filename = here("Figures", "US_sensitivity_Cal.png"), US_sensitivity_Cal, height = 3.5, width = 6, dpi = 600)
ggsave(filename = here("Figures", "US_sensitivity_Cal.pdf"), US_sensitivity_Cal, height = 3.5, width = 6, dpi = 600)

summary(resampled_pars$Sens_0)
quantile(resampled_pars$Sens_0, c(0.025, 0.975))

summary(resampled_pars$Sens_A)
quantile(resampled_pars$Sens_A, c(0.025, 0.975))


base_scen_fun <- function(`p_0`, p_A, p_B, p_C, p_D, Sens_0, Sens_A, Sym_A = Sym_A_def, Sym_B = Sym_B_def, p_HCC = p_HCC_def, ...) {
  life_add(
    p_HCC = p_HCC,
    `p_0` = `p_0`,
    p_A = p_A,
    p_B = p_B,
    p_C = p_C,
    p_D = p_D,
    Sens_0 = Sens_0,
    Sens_A = Sens_A,
    Sym_A = Sym_A,
    Sym_B = Sym_A * 8
  )
}

scen_1a_fun <- function(`p_0`, p_A, p_B, p_C, p_D, Sens_0, Sens_A, Sym_A = Sym_A_def, Sym_B = Sym_B_def, p_HCC = p_HCC_def, ...) {
  life_add(
    p_HCC = p_HCC,
    p_NA = p_NA_def - 0.05,
    p_CA = p_CA_def + 0.1,
    p_Att = p_Att_def + 0.1,
    `p_0` = `p_0`,
    p_A = p_A,
    p_B = p_B,
    p_C = p_C,
    p_D = p_D,
    Sens_0 = Sens_0,
    Sens_A = Sens_A,
    Sym_A = Sym_A,
    Sym_B = Sym_A * 8
  )
}

scen_1b_fun <- function(`p_0`, p_A, p_B, p_C, p_D, Sens_0, Sens_A, Sym_A = Sym_A_def, Sym_B = Sym_B_def, p_HCC = p_HCC_def, ...) {
  life_add(
    p_HCC = p_HCC,
    p_NA = p_NA_def - 0.1,
    p_CA = p_CA_def + 0.15,
    p_Att = p_Att_def + 0.2,
    `p_0` = `p_0`,
    p_A = p_A,
    p_B = p_B,
    p_C = p_C,
    p_D = p_D,
    Sens_0 = Sens_0,
    Sens_A = Sens_A,
    Sym_A = Sym_A,
    Sym_B = Sym_A * 8
  )
}

scen_2_fun <- function(`p_0`, p_A, p_B, p_C, p_D, Sens_0, Sens_A, Sym_A = Sym_A_def, Sym_B = Sym_B_def, p_HCC = p_HCC_def, ...) {
  life_add(
    p_HCC = p_HCC,
    `p_0` = `p_0`,
    p_A = p_A,
    p_B = p_B,
    p_C = p_C,
    p_D = p_D,
    Sens_0 = Sens_0 + 0.05,
    Sens_A = Sens_A + 0.05,
    Sens_B = Sens_B_def + 0.05,
    Sym_A = Sym_A,
    Sym_B = Sym_A * 8
  )
}

scen_3_fun <- function(`p_0`, p_A, p_B, p_C, p_D, Sens_0, Sens_A, Sym_A = Sym_A_def, Sym_B = Sym_B_def, p_HCC = p_HCC_def, ...) {
  life_add(
    p_HCC = p_HCC,
    L_A = L_A_def + 1,
    L_B = L_B_def + 0.5,
    L_C = L_C_def + 0.5,
    `p_0` = `p_0`,
    p_A = p_A,
    p_B = p_B,
    p_C = p_C,
    p_D = p_D,
    Sens_0 = Sens_0,
    Sens_A = Sens_A,
    Sym_A = Sym_A,
    Sym_B = Sym_A * 8
  )
}

resampled_pars <- resampled_pars %>%
  mutate(
    Baseline = pmap_dbl(., base_scen_fun),
    scen_1a = pmap_dbl(., scen_1a_fun),
    scen_1b = pmap_dbl(., scen_1b_fun),
    scen_2 = pmap_dbl(., scen_2_fun),
    scen_3 = pmap_dbl(., scen_3_fun),
    Baseline_low_risk = pmap_dbl(mutate(., p_HCC = p_HCC_def / 2), base_scen_fun),
    scen_1a_low_risk = pmap_dbl(mutate(., p_HCC = p_HCC_def / 2), scen_1a_fun),
    scen_1b_low_risk = pmap_dbl(mutate(., p_HCC = p_HCC_def / 2), scen_1b_fun),
    scen_2_low_risk = pmap_dbl(mutate(., p_HCC = p_HCC_def / 2), scen_2_fun),
    scen_3_low_risk = pmap_dbl(mutate(., p_HCC = p_HCC_def / 2), scen_3_fun),
    Baseline_high_risk = pmap_dbl( mutate(., p_HCC = p_HCC_def * 2), base_scen_fun),
    scen_1a_high_risk = pmap_dbl(mutate(., p_HCC = p_HCC_def * 2), scen_1a_fun),
    scen_1b_high_risk = pmap_dbl(mutate(., p_HCC = p_HCC_def * 2), scen_1b_fun),
    scen_2_high_risk = pmap_dbl(mutate(., p_HCC = p_HCC_def * 2), scen_2_fun),
    scen_3_high_risk = pmap_dbl(mutate(., p_HCC = p_HCC_def * 2), scen_3_fun)
  )

head(resampled_pars)

resampled_pars <- resampled_pars %>%
  mutate(
    scen_1a_diff = scen_1a - Baseline,
    scen_1b_diff = scen_1b - Baseline,
    scen_2_diff = scen_2 - Baseline,
    scen_3_diff = scen_3 - Baseline,
    scen_1a_3 = scen_1a - scen_2,
    scen_1a_4 = scen_1a - scen_3,
    scen_1b_3 = scen_1b - scen_2,
    scen_1b_4 = scen_1b - scen_3,
    scen_1a_diff_low_risk = scen_1a_low_risk - Baseline_low_risk,
    scen_1b_diff_low_risk = scen_1b_low_risk - Baseline_low_risk,
    scen_2_diff_low_risk = scen_2_low_risk - Baseline_low_risk,
    scen_3_diff_low_risk = scen_3_low_risk - Baseline_low_risk,
    scen_1a_diff_high_risk = scen_1a_high_risk - Baseline_high_risk,
    scen_1b_diff_high_risk = scen_1b_high_risk - Baseline_high_risk,
    scen_2_diff_high_risk = scen_2_high_risk - Baseline_high_risk,
    scen_3_diff_high_risk = scen_3_high_risk - Baseline_high_risk
  )

summary(resampled_pars$scen_1a_diff)
quantile(resampled_pars$scen_3_diff, c(0.025, 0.975))

pairs(head(resampled_pars, 1000) %>% select(c("scen_1a_diff", "scen_2_diff", "scen_3_diff")))
table(resampled_pars$scen_1a_diff > resampled_pars$scen_3_diff)

pairs(head(resampled_pars, 1000) %>% select(c("p_0", "p_A", "p_B", "p_C", "p_D")))


param_names <- readRDS(file = "param_names.rds")
param_names[c("opt_vs_diag", "opt_vs_treat", "opt_vs_real", "real_vs_diag", "real_vs_treat", "treat_vs_diag")] <- c(table(resampled_pars$scen_1b_diff > resampled_pars$scen_2_diff)[["TRUE"]] / 10000, table(resampled_pars$scen_1b_diff > resampled_pars$scen_3_diff)[["TRUE"]] / 10000, table(resampled_pars$scen_1b_diff > resampled_pars$scen_1a_diff)[["TRUE"]] / 10000, table(resampled_pars$scen_1a_diff > resampled_pars$scen_2_diff)[["TRUE"]] / 10000, table(resampled_pars$scen_1a_diff > resampled_pars$scen_3_diff)[["TRUE"]] / 10000, table(resampled_pars$scen_3_diff > resampled_pars$scen_2_diff)[["TRUE"]] / 10000)
saveRDS(param_names, file = "param_names.rds")


saveRDS(resampled_pars, file = "resampled_pars_fin.rds")

resampled_pars <- readRDS(file = "resampled_pars_fin.rds")

var_names <- colnames(resampled_pars)[c(1:9)]

# Uncertainty Analysis -------

ggplot(resampled_pars, aes(x = Baseline, y = "")) +
  geom_violin() +
  # geom_jitter(position = position_jitter(seed = 1, width = 0.2)) +
  theme_bw(base_family = "serif", base_size = 12)

ggplot(resampled_pars, aes(x = Baseline)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.5) +
  geom_vline(xintercept = life_add(), colour = "red") +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_y_continuous(expand = expansion(0), limits = c(0, 0.08)) +
  geom_function(fun = dnorm, args = list(mean = mean(resampled_pars$Baseline), sd = sd(resampled_pars$Baseline)), colour = "blue", linewidth = 1) +
  labs(x = "Life Years Gained", y = "Density")


ggsave(here("figures", "bootstrap_var_base.jpg"), width = 6, height = 3)
ggsave(here("figures", "bootstrap_var_base.png"), width = 6, height = 3, dpi = 600)
ggsave(here("figures", "bootstrap_var_base.pdf"), width = 6, height = 3, dpi = 600)

param_names <- readRDS(file = "param_names.rds")
param_names[c("sd_noint", "mean_noint")] <- c(sd(resampled_pars$Baseline), mean(resampled_pars$Baseline))
saveRDS(param_names, file = "param_names.rds")


ggplot(resampled_pars, aes(x = Baseline)) +
  geom_boxplot() +
  theme_bw(base_family = "serif", base_size = 12)

# Sensitivity/Uncertainty Analysis ----------

pars_def <- resampled_pars %>%
  dplyr::select(1:8, Baseline)

pars_logit <- resampled_pars %>%
  dplyr::select(1:8, Baseline) %>%
  mutate(across(-Baseline, Logit))

pairs(head(pars_logit, 1000))
pairs(head(pars_def, 1000))

pars_exp <- resampled_pars %>%
  dplyr::select(1:8, Baseline) %>%
  mutate(across(-Baseline, exp))

lm_def <- lm(pars_def$Baseline ~ ., data = pars_def %>% dplyr::select(-Baseline))
lm_logit <- lm(pars_logit$Baseline ~ ., data = pars_logit %>% dplyr::select(-Baseline))

summary.lm(lm_def)$r.squared
summary.lm(lm_logit)$r.squared

library(ppcor)
cor_matrix <- cor(pars_def)
# pcor_matrix <- pcor(pars_def)$estimate #error with inverse
row.names(pcor_matrix) <- var_names
# colnames(pcor_matrix) <- var_names

cor_matrix_logit <- cor(pars_logit)
pcor_matrix_logit <- pcor(pars_logit)$estimate



corrplot(
  cor_matrix,
  method = "color",
  type = "lower",
  diag = T,
  addCoef.col = "black",
  title = "Pearson correlation coefficient"
)

corrplot(
  cor_matrix_logit,
  method = "color",
  type = "lower",
  diag = T,
  addCoef.col = "black",
  title = "Logit Pearson correlation coefficient"
)

# corrplot(
#   pcor_matrix,
#   method = "color",
#   type = "lower",
#   diag = T,
#   addCoef.col = "black",
#   title = "Partial correlation coefficient"
# )

corrplot(
  pcor_matrix_logit,
  method = "color",
  type = "lower",
  diag = T,
  addCoef.col = "black",
  title = "Logit Partial correlation coefficient"
)

long_pcor <- pcor_matrix %>%
  as_tibble() %>%
  mutate("row" = fct_inorder(var_names)) %>%
  pivot_longer(-row, names_to = "col") %>%
  mutate(col = fct_inorder(col))

ggplot(long_pcor, aes(x = col, y = row, fill = value)) +
  geom_tile() +
  scale_fill_gradient2() +
  scale_y_discrete(
    expand = expansion(0),
    limits = rev,
    labels = TeX_labels
  ) +
  scale_x_discrete(
    expand = expansion(0),
    position = "top",
    labels = TeX_labels
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  geom_text(aes(label = round(value, 2)), color = "black") +
  theme(
    axis.text.x = element_text(
      angle = 60,
      size = 12,
      hjust = 0
    ),
    axis.text.y = element_text(
      size = 12
    )
  )


long_pcor_logit <- pcor_matrix_logit %>%
  as_tibble() %>%
  mutate("row" = fct_inorder(var_names)) %>%
  pivot_longer(-row, names_to = "col") %>%
  mutate(col = fct_inorder(col))


plot_par_cor_logit <- ggplot(long_pcor_logit, aes(x = col, y = row, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    limit = c(-1, 1),
    low = "red",
    high = "#075AFF"
  ) +
  scale_y_discrete(
    expand = expansion(0),
    limits = rev
  ) +
  scale_x_discrete(
    expand = expansion(0),
    position = "top"
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  geom_text(aes(label = round(value, 2)), color = "black") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      # size = 12,
      hjust = 0
    ),
    # axis.text.y = element_text(
    #   size = 12
    # )
  ) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(title = "Partial Rank\nCoefficient"))

plot_par_cor_logit

ggsave(plot = plot_par_cor_logit, here("figures/plot_par_cor_logit.jpg"), width = 6, height = 5)

pcor_base_logit <- long_pcor_logit %>%
  filter(row == "Baseline" & col != "Baseline") %>%
  rename(Variable = col)




plot_cor_base_logit <- ggplot(pcor_base_logit, aes(x = Variable, y = value, fill = value)) +
  geom_col(colour = "black") +
  scale_fill_gradient2(
    limit = c(-1, 1),
    low = "red",
    high = "#075AFF"
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  # theme(
  #   axis.text.x = element_text(
  #     angle = 30,
  #     hjust = 1
  #   )
  # ) +
  labs(
    y = "Partial Correlation Coefficient",
    x = NULL
  ) +
  guides(fill = guide_colourbar(title = NULL)) +
  scale_y_continuous(limits = c(-1, 1), expand = expansion(0)) +
  scale_x_discrete(labels = TeX_labels)

plot_cor_base_logit

ggsave(plot = plot_cor_base_logit, here("figures", "plot_cor_base_logit.jpg"), width = 6, height = 3)
ggsave(plot = plot_cor_base_logit, here("figures", "plot_cor_base_logit.png"), width = 6, height = 3, dpi = 600)
ggsave(plot = plot_cor_base_logit, here("figures", "plot_cor_base_logit.pdf"), width = 6, height = 3, dpi = 600)

detach("package:ppcor", unload = TRUE)
detach("package:MASS", unload = TRUE)

param_names <- readRDS(file = "param_names.rds")
param_names[c("pcor_logit_0", "pcor_logit_US_0", "pcor_logit_US_A")] <- c(pcor_matrix_logit["Baseline", "p_0"], pcor_matrix_logit["Baseline", "Sens_0"], pcor_matrix_logit["Baseline", "Sens_A"])
saveRDS(param_names, file = "param_names.rds")


default_param_tbl <- tibble(params = list(hypo_probs)) %>% unnest_wider(params, names_sep = "t", names_repair = ~ c("p_0", "p_A", "p_B", "p_C", "p_D", "Sens_0", "Sens_A", "Sym_A"))

results_tbl <- resampled_pars %>%
  # select(-c(gradient, value, convergence, `function`)) %>%
  pivot_longer(everything(), names_to = "result") %>%
  group_by(result) %>%
  summarise(point = median(value), lower = quantile(value, 0.025), upper = quantile(value, 0.975))

results_tbl <- results_tbl %>%
  rows_update(
    tibble(
      result = c(
        "scen_1a_diff", "scen_1b_diff", "scen_2_diff", "scen_3_diff",
        "scen_1a", "scen_1b", "scen_2", "scen_3",
        "scen_1a_3", "scen_1a_4", "scen_1b_3", "scen_1b_4"
      ),
      point = c(
        pmap_dbl(default_param_tbl, scen_1a_fun) - pmap_dbl(default_param_tbl, base_scen_fun),
        pmap_dbl(default_param_tbl, scen_1b_fun) - pmap_dbl(default_param_tbl, base_scen_fun),
        pmap_dbl(default_param_tbl, scen_2_fun) - pmap_dbl(default_param_tbl, base_scen_fun),
        pmap_dbl(default_param_tbl, scen_3_fun) - pmap_dbl(default_param_tbl, base_scen_fun),
        pmap_dbl(default_param_tbl, scen_1a_fun), pmap_dbl(default_param_tbl, scen_1b_fun), pmap_dbl(default_param_tbl, scen_2_fun), pmap_dbl(default_param_tbl, scen_3_fun),
        pmap_dbl(default_param_tbl, scen_1a_fun) - pmap_dbl(default_param_tbl, scen_2_fun),
        pmap_dbl(default_param_tbl, scen_1a_fun) - pmap_dbl(default_param_tbl, scen_3_fun),
        pmap_dbl(default_param_tbl, scen_1b_fun) - pmap_dbl(default_param_tbl, scen_2_fun),
        pmap_dbl(default_param_tbl, scen_1b_fun) - pmap_dbl(default_param_tbl, scen_3_fun)
      )
    ),
    by = "result"
  )

library(patchwork)
library(cowplot)


# saveRDS(resampled_pars, file = "resampled_pars_5.rds")
#
# resampled_pars <- readRDS(file = "resampled_pars_5.rds")

HCC_Time_Range_Plot_side_labs <- ggplot(
  mean_years_tbl %>%
    group_by(name) %>%
    summarise(mean025 = quantile(mean_years, 0.025), mean975 = quantile(mean_years, 0.975)),
  aes(name, ymin = mean025, ymax = mean975)
) +
  geom_errorbar() +
  geom_point(mapping = aes(x = rate, y = 1 / Point / 12), data = rates_tbl, inherit.aes = F) +
  scale_x_discrete(labels = c(
    p_0 = "BCLC Stage 0",
    p_A = "BCLC Stage A",
    p_B = "BCLC Stage B",
    p_C = "BCLC Stage C",
    p_D = "BCLC Stage D"
  )) +
  labs(
    x = NULL,
    y = "Mean Number of\nYears In BCLC Stage"
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_y_continuous(limits = c(0, 3), expand = expansion(0)) +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1
    )
  )

HCC_Time_Range_Plot_side_labs

plot_cor_base_logit_side_labs <- ggplot(pcor_base_logit, aes(x = Variable, y = value, fill = value)) +
  geom_col(colour = "black") +
  scale_fill_gradient2(
    limit = c(-1, 1),
    low = "red",
    high = "#075AFF"
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  # theme(
  #   axis.text.x = element_text(
  #     angle = 30,
  #     hjust = 1
  #   )
  # ) +
  labs(
    y = "Partial Correlation Coefficient",
    x = NULL
  ) +
  guides(fill = guide_colourbar(title = NULL)) +
  scale_y_continuous(limits = c(-1, 1), expand = expansion(0)) +
  scale_x_discrete(labels = TeX_labels) +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1
    ),
    legend.position = "bottom"
  )

US_sensitivity_Cal_side_labs <- ggplot(resampled_pars, aes(x = Sens_0, y = Sens_A)) +
  # geom_point() +
  # geom_density2d_filled(alpha = 0.75) +
  # geom_raster(stat = "density") +
  geom_area(aes(x = seq(0, 1, length.out = 10000), y = seq(0, 1, length.out = 10000)), alpha = 0.25) +
  geom_hex(alpha = 0.75, bins = c(100, 80)) +
  # geom_density2d() +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_fill_viridis(
    # discrete = T,
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
  geom_point(aes(x = Sens_0_def, y = Sens_A_def), shape = 4, stroke = 1, size = 1, colour = "black") +
  labs(
    # title = "Ultrasound Sensitivities Calibration",
    x = "BCLC Stage 0 HCC\nUltrasound Sensitivity",
    y = "BCLC Stage A HCC\nUltrasound Sensitivity",
    fill = "Number of\nSamples"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

multipanel_cals <- plot_grid(HCC_Time_Range_Plot_side_labs, US_sensitivity_Cal_side_labs, plot_cor_base_logit_side_labs, plot_univar_mod, labels = "AUTO", rel_widths = c(1, 2))

multipanel_cals

ggsave(here("figures", "multipanel_cals.png"), width = 6, height = 8, dpi = 600)
ggsave(here("figures", "multipanel_cals.pdf"), width = 6, height = 8, dpi = 600)
