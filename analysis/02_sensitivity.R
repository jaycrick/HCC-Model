# analysis/02_sensitivity.R
#
# Univariate sensitivity analyses for the HCC surveillance model.
# Reads Parameters.csv and Parameters_25.csv from data/.
# Writes figures to figures/ and updates data-derived/param_names.rds.
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

dir.create(here("figures"), showWarnings = FALSE)

# ----- Sensitivity data -------------------------------------------------------
param_tbl <- life_tibble(here("data", "Parameters.csv"))

facet_labels <- c(
  L = "Treatment Quality",
  Pr = "Adherence\nto Screening",
  p = "Monthly Risk",
  Sym = "Monthly probability\nof symptoms",
  Sens = "Ultrasound\nSensitivity",
  years_low_comp = "Lower Rate of HCC",
  years_high_comp = "Higher Rate of HCC"
)

variable_labels <- c(
  "L_0" = "Stage 0",
  "L_A" = "Stage A",
  "L_B" = "Stage B",
  "L_C" = "Stage C",
  "L_D" = "Stage D",
  "p_CA" = "Completely adherent",
  "p_NA" = "Non-adherent",
  "p_Att" = "Attendence probability\nwithin partial adherents",
  "p_HCC" = "Developing HCC",
  "p_0" = "Stage 0 to A",
  "p_A" = "Stage A to B HCC",
  "p_B" = "Stage B to C HCC",
  "p_C" = "Stage C to D HCC",
  "p_D" = "Stage D to Death",
  "Sym_0" = "Stage 0",
  "Sym_A" = "Stage A",
  "Sym_B" = "Stage B",
  "Sym_C" = "Stage C",
  "Sym_D" = "Stage D",
  "Sens_0" = "Stage 0",
  "Sens_A" = "Stage A",
  "Sens_B" = "Stage B",
  "Sens_C" = "Stage C",
  "Sens_D" = "Stage D"
)

# ----- Univariate sensitivity tornado ----------------------------------------
plot_univariate <- ggplot(
  data = param_tbl %>%
    mutate(
      Parameter = fct_reorder(.f = Parameter, .x = years_mid_comp, abssum),
      Type = factor(
        word(Parameter, sep = fixed("_")),
        levels = c("Pr", "Sens", "L", "p", "Sym")
      )
    ) %>%
    dplyr::filter(years_mid_comp != 0),
  aes(y = Parameter, x = years_mid_comp)
) +
  geom_col(position = "identity", colour = "black") +
  labs(
    x = "Additional Life Years Gained From\nIntervention per Cohort of 100",
    y = "Variable",
    title = "Univariate Sensitivity Analysis"
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_discrete(labels = variable_labels) +
  facet_wrap(
    ~Type,
    ncol = 1,
    scales = "free_y",
    labeller = labeller(Type = facet_labels)
  )

plot_univariate

ggsave(
  plot = plot_univariate,
  here("figures", "sens_anal_univariate.jpg"),
  width = 5.5,
  height = 8
)

# ----- f4_0 upper graph -------------------------------------------------------
plot_upper <- ggplot(
  data = param_tbl %>%
    select(!(value:years_high)) %>%
    mutate(
      Parameter = fct_reorder(.f = Parameter, .x = years_low_comp, abssum),
      Type = factor(case_when(
        str_detect(Parameter, "Pr") ~ "pr",
        str_detect(Parameter, "US") ~ "US"
      ))
    ) %>%
    dplyr::filter(years_high_comp != 0, Parameter != "p_HCC"),
  aes(y = Parameter, x = years_mid_comp)
) +
  geom_col(position = "identity", colour = "black") +
  labs(
    x = "Additional Life Years Gained From\nIntervention per Cohort of 100",
    y = "Variable",
    title = TeX(r'(Upper $\textit{p}_{HCC}$ Sensitivity Analysis)')
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_discrete(labels = TeX_labels)

plot_upper

ggsave(
  plot = plot_upper,
  here("figures", "sens_anal_f4_0_upper.jpg"),
  width = 5.5,
  height = 9
)

# ----- f4_0 lower graph -------------------------------------------------------
plot_lower <- ggplot(
  data = param_tbl %>%
    mutate(
      Parameter = fct_reorder(.f = Parameter, .x = years_low_comp, abssum)
    ) %>%
    dplyr::filter(years_low_comp != 0, Parameter != "p_HCC"),
  aes(y = Parameter, x = years_low_comp)
) +
  geom_col(position = "identity", colour = "black") +
  labs(
    x = "Additional Life Years Gained From\nIntervention per Cohort of 100",
    y = "Variable",
    title = TeX(r"(Lower $\textit{p}_{HCC}$ Sensitivity Analysis)")
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_discrete(labels = TeX_labels)

plot_lower

ggsave(
  plot = plot_lower,
  here("figures", "sens_anal_f4_0_lower.jpg"),
  width = 5.5,
  height = 9
)

# ----- Combined upper and lower graph -----------------------------------------
facet_sens <- ggplot(
  data = param_tbl %>%
    mutate(
      Parameter = fct_reorder(.f = Parameter, .x = years_low_comp, abssum),
      Type = factor(
        word(Parameter, sep = fixed("_")),
        levels = c("Pr", "Sens", "L", "p", "Sym")
      )
    ) %>%
    dplyr::filter(years_low_comp != 0, Parameter != "p_HCC") %>%
    dplyr::select(c(
      "Parameter",
      "Bound",
      "years_low_comp",
      "years_high_comp",
      "Type"
    )) %>%
    pivot_longer(
      c("years_low_comp", "years_high_comp"),
      names_to = "years_comp"
    ),
  aes(y = Parameter, x = value)
) +
  geom_col(position = "identity", colour = "black") +
  labs(
    x = "Additional Life Years Gained From\nIntervention per Cohort of 100",
    y = "Variable",
    title = "Sensitivity Analysis"
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  scale_y_discrete(labels = variable_labels) +
  facet_grid(
    vars(Type),
    vars(years_comp),
    scales = "free",
    labeller = labeller(Type = facet_labels, years_comp = facet_labels)
  )
facet_sens

ggsave(
  plot = facet_sens,
  here("figures", "sens_anal_f4_0_both.jpg"),
  width = 5.5,
  height = 8
)

# ----- 25%-variation univariate sensitivity -----------------------------------
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

merged_names <- map2_chr(
  param_chang_tbl$Parameter,
  param_chang_tbl$Bound,
  ~ paste(.x, .y, sep = "_")
)

param_names <- readRDS(here("data-derived", "param_names.rds"))
param_names[merged_names] <- param_chang_tbl$years_mid_comp
saveRDS(param_names, file = here("data-derived", "param_names.rds"))

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

plot_univar_mod

ggsave(
  plot = plot_univar_mod,
  here("figures", "sens_anal_mod.jpg"),
  width = 5,
  height = 4
)
ggsave(
  plot = plot_univar_mod,
  here("figures", "sens_anal_mod.png"),
  width = 5,
  height = 4,
  dpi = 600
)
ggsave(
  plot = plot_univar_mod,
  here("figures", "sens_anal_mod.pdf"),
  width = 5,
  height = 4,
  dpi = 600
)

# ----- Full 25% univariate (all parameters) ----------------------------------
plot_univariate_25 <- ggplot(
  data = param_tbl_25 %>%
    mutate(
      Parameter = fct_reorder(.f = Parameter, .x = years_mid_comp, abssum),
      Type = factor(
        word(Parameter, sep = fixed("_")),
        levels = c("Pr", "Sens", "L", "p", "Sym")
      )
    ),
  aes(y = Parameter, x = years_mid_comp)
) +
  geom_col(position = "identity", colour = "black") +
  labs(x = "Additional Life Years Gained", y = NULL, title = NULL) +
  theme_bw(base_family = "serif", base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_discrete(labels = variable_labels) +
  facet_wrap(
    ~Type,
    ncol = 1,
    scales = "free_y",
    labeller = labeller(Type = facet_labels)
  )

plot_univariate_25

ggsave(
  plot = plot_univariate_25,
  here("figures", "sens_anal_25_univariate.jpg"),
  width = 5,
  height = 8
)

# ----- 25% f4_0 upper graph --------------------------------------------------
plot_upper_25 <- ggplot(
  data = param_tbl_25 %>%
    dplyr::filter(
      years_high_comp != 0,
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
      Parameter = fct_reorder(.f = Parameter, .x = years_high_comp, abssum)
    ),
  aes(y = Parameter, x = years_high_comp)
) +
  geom_col(position = "identity", colour = "black") +
  labs(
    x = "Additional Life Years Gained From\nIntervention per Cohort of 100",
    y = "Variable",
    title = TeX(
      r"(Upper $\textit{p}_{HCC}$ Estimated change in mortality from 25% change in each modifiable parameter)"
    )
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_discrete(labels = TeX_labels)

plot_upper_25

ggsave(
  plot = plot_upper_25,
  here("figures", "sens_anal_25_f4_0_upper.jpg"),
  width = 5.5,
  height = 9
)

# ----- 25% f4_0 lower graph --------------------------------------------------
plot_lower_25 <- ggplot(
  data = param_tbl_25 %>%
    dplyr::filter(
      years_low_comp != 0,
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
      Parameter = fct_reorder(.f = Parameter, .x = years_low_comp, abssum)
    ),
  aes(x = years_low_comp, y = Parameter)
) +
  geom_col(position = "identity", colour = "black") +
  labs(
    x = "Additional Life Years Gained From\nIntervention per Cohort of 100",
    y = "Variable",
    title = TeX(
      r"(Lower $\textit{p}_{HCC}$ Estimated change in mortality from 25% change\nin each modifiable parameter)"
    )
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_discrete(labels = TeX_labels)

plot_lower_25

ggsave(
  plot = plot_lower_25,
  here("figures", "sens_anal_25_f4_0_lower.jpg"),
  width = 5.5,
  height = 9
)

# ----- Regression anchors (verified against current Parameters.csv) ----------
# life_add() = 68.02326 with current default parameters
stopifnot(isTRUE(all.equal(life_add(), 68.02326, tolerance = 0.01)))
message("Regression check passed: life_add() = ", round(life_add(), 4))
