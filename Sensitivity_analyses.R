# Initialise ----------------------------
library(expm)
library(tidyverse)
library(here)
library(viridis)
library(latex2exp)
library(dplyr)
# detach("package:MASS", unload = TRUE)
set.seed(4567385) # generated from random.org between 1 and 10000000

stage_list <- c(
  "Cirrhosis",
  "BCLC 0",
  "BCLC A",
  "BCLC B",
  "BCLC C",
  "BCLC D",
  "Mortality"
)

param_def <- read.csv(file = "Parameters.csv", row.names = 1)[1:3] # the default parameters for all functions

for (i in row.names(param_def)) {
  # assigning ***_def to variable's default (point) estimates
  # assign(i, param_def[i, "POINT"])
  assign(paste(i, "def", sep = "_"), param_def[i, "POINT"])
}

# FUNCTIONS ----------
# BASE MODEL FUNCTION ------------
fromf4 <- function(
  n,
  t, # Given n people, how many will progress to what after t months
  p_HCC = p_HCC_def, # NEEDS TO BE ANNUAL
  p_0 = p_0_def,
  p_A = p_A_def,
  p_B = p_B_def,
  p_C = p_C_def,
  p_D = p_D_def
) {
  p_HCC <- 1 - (1 - p_HCC)^(1 / 12)
  p_0 <- 1 - (1 - p_0)^(1 / 12)
  p_A <- 1 - (1 - p_A)^(1 / 12)
  p_B <- 1 - (1 - p_B)^(1 / 12)
  p_C <- 1 - (1 - p_C)^(1 / 12)
  p_D <- 1 - (1 - p_D)^(1 / 12)
  tmat <- t(matrix(
    c(
      1 - p_HCC,
      p_HCC,
      rep(0, 6),
      1 - p_0,
      p_0,
      rep(0, 6),
      1 - p_A,
      p_A,
      rep(0, 6),
      1 - p_B,
      p_B,
      rep(0, 6),
      1 - p_C,
      p_C,
      rep(0, 6),
      1 - p_D,
      p_D,
      rep(0, 6),
      1
    ),
    7,
    7
  ))
  stage <- c(n, rep(0, 6)) %*% (tmat %^% t)
  tbl <- as_tibble(stage, .name_repair = ~stage_list)
  return(tbl)
}

fromf4(100, 120)


# USs --------------
us_theor_pa <- function(
  n = 100,
  t = 120,
  m = 6,
  Pr_Att, # Pr_Att is 1 for those who are comp. adh. and 0 for non
  p_HCC = p_HCC_def, # NEEDS TO BE ANNUAL
  p_0 = p_0_def,
  p_A = p_A_def,
  p_B = p_B_def,
  p_C = p_C_def,
  p_D = p_D_def,
  Sens_0 = Sens_0_def,
  Sens_A = Sens_A_def,
  Sens_B = Sens_B_def,
  Sens_C = Sens_C_def,
  Sens_D = Sens_D_def,
  Sym_0 = Sym_0_def,
  Sym_A = Sym_A_def,
  Sym_B = Sym_B_def,
  Sym_C = Sym_C_def,
  Sym_D = Sym_D_def
) {
  # theoretical version of us_sim_pa
  p_HCC <- 1 - (1 - p_HCC)^(1 / 12)
  tmat <- t(matrix(
    c(
      1 - p_HCC,
      p_HCC,
      rep(0, 6),
      1 - p_0,
      p_0,
      rep(0, 6),
      1 - p_A,
      p_A,
      rep(0, 6),
      1 - p_B,
      p_B,
      rep(0, 6),
      1 - p_C,
      p_C,
      rep(0, 6),
      1 - p_D,
      p_D,
      rep(0, 6),
      1
    ),
    7,
    7
  )) # Creates a transition matrix,
  vec <- c(n, rep(0, 6))
  Dx <- rep(0, 7) # diagnosed people
  sym_prob <- c(
    0,
    Sym_0,
    Sym_A,
    Sym_B,
    Sym_C,
    Sym_D,
    1
  )
  US_probs <- c(
    0,
    Sens_0,
    Sens_A,
    Sens_B,
    Sens_C,
    Sens_D,
    1
  )
  for (i in 1:t) {
    vec <- vec %*% tmat
    sym <- vec * sym_prob
    Dx <- Dx + sym
    vec <- vec - sym
    if (i %% m == 0) {
      go <- c(
        0,
        vec[2] * Pr_Att,
        vec[3] * Pr_Att,
        vec[4] * Pr_Att,
        vec[5] * Pr_Att,
        vec[6] * Pr_Att,
        vec[7]
      ) # how many actually go to the ultrasound
      us <- go * US_probs
      Dx <- Dx + us
      vec <- vec - us
    }
  }
  return(Dx)
}

# Adherence model ----------

adh_theor <- function(
  n = 100,
  t = 120,
  m = 6,
  p_NA = p_NA_def,
  p_CA = p_CA_def,
  p_Att = p_Att_def,
  p_HCC = p_HCC_def, # NEEDS TO BE ANNUAL
  p_0 = p_0_def,
  p_A = p_A_def,
  p_B = p_B_def,
  p_C = p_C_def,
  p_D = p_D_def,
  Sens_0 = Sens_0_def,
  Sens_A = Sens_A_def,
  Sens_B = Sens_B_def,
  Sens_C = Sens_C_def,
  Sens_D = Sens_D_def,
  Sym_0 = Sym_0_def,
  Sym_A = Sym_A_def,
  Sym_B = Sym_B_def,
  Sym_C = Sym_C_def,
  Sym_D = Sym_D_def,
  ...
) {
  # number of people, total time tracking, interval of US, uses us_theor_pa

  adh_group <- n * c(p_NA, 1 - p_NA - p_CA, p_CA) # number of people in each adherence group
  custom_us_theor <- function(n, Pr_Att) {
    # FUNCTION THAT JUST VARIES NUMBERS AND THEIR LIKELIHOOD TO ATTEND
    us_theor_pa(
      n,
      Pr_Att,
      t = t,
      m = m,
      p_HCC = p_HCC,
      p_0 = p_0,
      p_A = p_A,
      p_B = p_B,
      p_C = p_C,
      p_D = p_D,
      Sens_0 = Sens_0,
      Sens_A = Sens_A,
      Sens_B = Sens_B,
      Sens_C = Sens_C,
      Sens_D = Sens_D,
      Sym_0 = Sym_0,
      Sym_A = Sym_A,
      Sym_B = Sym_B,
      Sym_C = Sym_C,
      Sym_D = Sym_D
    ) %>%
      as_tibble(.name_repair = ~stage_list)
  }

  stage_caught <-
    map2_dfr(
      adh_group,
      c(0, p_Att, 1),
      .f = custom_us_theor
    ) %>%
    summarise_all(sum)
  # %>%
  #  as.vector(x = ., mode = "double")

  return(stage_caught)
}

# YEARS LIFE ADDED BY MEDICAL TREATMENT ----------

life_add <- function(
  n = 100,
  t = 120,
  m = 6,
  p_NA = p_NA_def,
  p_CA = p_CA_def,
  p_Att = p_Att_def,
  p_HCC = p_HCC_def, # NEEDS TO BE ANNUAL
  p_0 = p_0_def,
  p_A = p_A_def,
  p_B = p_B_def,
  p_C = p_C_def,
  p_D = p_D_def,
  Sens_0 = Sens_0_def,
  Sens_A = Sens_A_def,
  Sens_B = Sens_B_def,
  Sens_C = Sens_C_def,
  Sens_D = Sens_D_def,
  Sym_0 = Sym_0_def,
  Sym_A = Sym_A_def,
  Sym_B = Sym_B_def,
  Sym_C = Sym_C_def,
  Sym_D = Sym_D_def,
  L_0 = L_0_def,
  L_A = L_A_def,
  L_B = L_B_def,
  L_C = L_C_def,
  L_D = L_D_def,
  ...
) {
  stage_caught <- adh_theor(
    n = n,
    m = m,
    t = t,
    p_NA = p_NA,
    p_CA = p_CA,
    p_Att = p_Att,
    p_HCC = p_HCC, # NEEDS TO BE ANNUAL
    p_0 = p_0,
    p_A = p_A,
    p_B = p_B,
    p_C = p_C,
    p_D = p_D,
    Sens_0 = Sens_0,
    Sens_A = Sens_A,
    Sens_B = Sens_B,
    Sens_C = Sens_C,
    Sens_D = Sens_D,
    Sym_0 = Sym_0,
    Sym_A = Sym_A,
    Sym_B = Sym_B,
    Sym_C = Sym_C,
    Sym_D = Sym_D
  )
  life_added <- sum(stage_caught * c(0, L_0, L_A, L_B, L_C, L_D, 0))
  return(life_added)
}

life_add() # 68.02326

dfs <- function(Parameter, value, ...) {
  new_df <- tibble(Parameter, value) %>%
    pivot_wider(names_from = Parameter, values_from = value)
  return(new_df)
}

years <- function(df) pmap_dbl(df, life_add)

# tibble creation for graphs ---------

life_tibble <- function(csv) {
  param_csv <- read.csv(file = csv, row.names = 1)
  param_tidy <- param_csv[1:24, 2:3] %>%
    as.data.frame() %>%
    rownames_to_column(var = "Parameter") %>%
    as_tibble() %>%
    pivot_longer(-Parameter, names_to = "Bound")
  param_tbl <- param_tidy %>%
    mutate(
      dfs_new = pmap(., dfs),
      Bound = as.factor(str_to_title(Bound)),
      Parameter = as.factor(Parameter),
      dfs_low = map(
        dfs_new,
        ~ if (names(.x) == "p_HCC") .x else
          tibble(.x, p_HCC = param_csv["p_HCC", "LOWER"])
      ),
      dfs_high = map(
        dfs_new,
        ~ if (names(.x) == "p_HCC") .x else
          tibble(.x, p_HCC = param_csv["p_HCC", "UPPER"])
      ),
      years_mid = map_dbl(dfs_new, years),
      years_low = map_dbl(dfs_low, years),
      years_high = map_dbl(dfs_high, years),
      years_mid_comp = years_mid - life_add(),
      years_low_comp = years_low -
        life_add(p_HCC = param_csv["p_HCC", "LOWER"]),
      years_high_comp = years_high -
        life_add(p_HCC = param_csv["p_HCC", "UPPER"])
    )
  return(param_tbl)
}


# GRAPHS ----------
# UNIVARIATE SENSITIVITY -------

param_tbl <- life_tibble("Parameters.csv")

abssum <- function(x) sum(abs(x))

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

TeX_labels <- c(
  "L_0" = TeX(r"($\textit{L}_{0}$)"),
  "L_A" = TeX(r"($\textit{L}_{A}$)"),
  "L_B" = TeX(r"($\textit{L}_{B}$)"),
  "L_C" = TeX(r"($\textit{L}_{C}$)"),
  "L_D" = TeX(r"($\textit{L}_{D}$)"),
  "p_CA" = TeX(r"($\textit{p}_{CA}$)"),
  "p_NA" = expression(italic(p)[`NA`]),
  "p_Att" = TeX(r"($\textit{p}_{Att}$)"),
  "p_HCC" = TeX(r"($\textit{p}_{HCC}$)"),
  "p_0" = TeX(r"($\textit{p}_{0}$)"),
  "p_A" = TeX(r"($\textit{p}_{A}$)"),
  "p_B" = TeX(r"($\textit{p}_{C}$)"),
  "p_C" = TeX(r"($\textit{p}_{C}$)"),
  "p_D" = TeX(r"($\textit{p}_{D}$)"),
  "Sym_0" = TeX(r"($Sym_{0}$)"),
  "Sym_A" = TeX(r"($Sym_{A}$)"),
  "Sym_B" = TeX(r"($Sym_{B}$)"),
  "Sym_C" = TeX(r"($Sym_{C}$)"),
  "Sym_D" = TeX(r"($Sym_{D}$)"),
  "Sens_0" = TeX(r"($Sens_{0}$)"),
  "Sens_A" = TeX(r"($Sens_{A}$)"),
  "Sens_B" = TeX(r"($Sens_{B}$)"),
  "Sens_C" = TeX(r"($Sens_{C}$)"),
  "Sens_D" = TeX(r"($Sens_{D}$)")
)

Boot_labels <- c(
  "m_0" = TeX(r"($\textit{m}^{*}_{0}$)"),
  "m_A" = TeX(r"($\textit{m}^{*}_{A}$)"),
  "m_B" = TeX(r"($\textit{m}^{*}_{B}$)"),
  "m_C" = TeX(r"($\textit{m}^{*}_{C}$)"),
  "m_D" = TeX(r"($\textit{m}^{*}_{D}$)"),
  "D_0" = TeX(r"($\textit{D}^{*}_{0}$)"),
  "D_A" = TeX(r"($\textit{D}^{*}_{A}$)"),
  "D_alpha" = TeX(r"($\textit{D}^{*}_{\alpha}$)"),
  "p_0" = TeX(r"($\textit{p}^{*}_{0}$)"),
  "p_A" = TeX(r"($\textit{p}^{*}_{A}$)"),
  "p_B" = TeX(r"($\textit{p}^{*}_{B}$)"),
  "p_C" = TeX(r"($\textit{p}^{*}_{C}$)"),
  "p_D" = TeX(r"($\textit{p}^{*}_{D}$)")
)

plot_univariate <- ggplot(
  data = param_tbl %>%
    mutate(
      Parameter = fct_reorder(.f = Parameter, .x = years_mid_comp, abssum),
      Type = factor(
        word(Parameter, sep = fixed("_")),
        levels = c("Pr", "Sens", "L", "p", "Sym")
      )
    ) %>%
    filter(years_mid_comp != 0),
  aes(
    y = Parameter,
    x = years_mid_comp,
    # fill = Bound
  )
) +
  geom_col(
    position = "identity",
    colour = "black"
  ) +
  # scale_fill_viridis(discrete = T) +
  labs(
    x = "Additional Life Years Gained From\nIntervention per Cohort of 100",
    y = "Variable",
    title = "Univariate Sensitivity Analysis",
    # caption = "Variation from the point estimate\nin 100 people over 10 years"
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
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


# f4_0 UPPER GRAPH ------------
# library(tikzDevice)
# tikz('math.tex', standAlone = TRUE, width = 4, height = 4)

plot_upper <- ggplot(
  data = param_tbl %>%
    select(!(value:years_high)) %>%
    mutate(
      Parameter = fct_reorder(.f = Parameter, .x = years_low_comp, abssum), # filtering by low_comp in order to keep both graphs same order
      Type = factor(case_when(
        str_detect(Parameter, "Pr") ~ "pr",
        str_detect(Parameter, "US") ~ "US"
      ))
    ) %>%
    filter(years_high_comp != 0, Parameter != "p_HCC"),
  aes(
    y = Parameter,
    x = years_mid_comp,
    # fill = Bound
  )
) +
  geom_col(
    position = "identity",
    colour = "black"
  ) +
  # scale_fill_viridis(discrete = T) +
  labs(
    x = "Additional Life Years Gained From\nIntervention per Cohort of 100",
    y = "Variable",
    title = TeX(r'(Upper $\textit{p}_{HCC}$ Sensitivity Analysis)'),
    # caption = "Variation from the point estimate\nin 100 people over 10 years"
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_y_discrete(labels = TeX_labels)

plot_upper

ggsave(
  plot = plot_upper,
  here("figures", "sens_anal_f4_0_upper.jpg"),
  width = 5.5,
  height = 9
)


# f4_0 LOWER GRAPH ------------

plot_lower <- ggplot(
  data = param_tbl %>%
    mutate(
      Parameter = fct_reorder(.f = Parameter, .x = years_low_comp, abssum)
    ) %>%
    filter(years_low_comp != 0, Parameter != "p_HCC"),
  aes(
    y = Parameter,
    x = years_low_comp,
    # fill = Bound
  )
) +
  geom_col(
    position = "identity",
    colour = "black"
  ) +
  # scale_fill_viridis(discrete = T) +
  labs(
    x = "Additional Life Years Gained From\nIntervention per Cohort of 100",
    y = "Variable",
    title = TeX(r"(Lower $\textit{p}_{HCC}$ Sensitivity Analysis)"),
    # caption = "Variation from the point estimate\nin 100 people over 10 years"
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_y_discrete(labels = TeX_labels)

plot_lower

ggsave(
  plot = plot_lower,
  here("figures", "sens_anal_f4_0_lower.jpg"),
  width = 5.5,
  height = 9
)

# combined upper and lower graphs --------

facet_sens <- ggplot(
  data = param_tbl %>%
    mutate(
      Parameter = fct_reorder(.f = Parameter, .x = years_low_comp, abssum),
      Type = factor(
        word(Parameter, sep = fixed("_")),
        levels = c("Pr", "Sens", "L", "p", "Sym")
      )
    ) %>%
    filter(years_low_comp != 0, Parameter != "p_HCC") %>%
    select(c(
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
  aes(
    y = Parameter,
    x = value,
    # fill = Bound
  )
) +
  geom_col(
    position = "identity",
    colour = "black"
  ) +
  # scale_fill_viridis(discrete = T) +
  labs(
    x = "Additional Life Years Gained From\nIntervention per Cohort of 100",
    y = "Variable",
    title = "Sensitivity Analysis",
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  ) +
  scale_y_discrete(labels = variable_labels) +
  facet_grid(
    vars(Type),
    vars(years_comp),
    scales = "free",
    labeller = labeller(
      Type = facet_labels,
      years_comp = facet_labels
    )
  )
facet_sens

ggsave(
  plot = facet_sens,
  here("figures", "sens_anal_f4_0_both.jpg"),
  width = 5.5,
  height = 8
)

# PM 25% version univariate --------------

param_tbl_25 <- life_tibble("Parameters_25.csv")

param_chang_tbl <- param_tbl_25 %>%
  filter(
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
  select(Parameter, Bound, years_mid_comp, Type)

merged_names = map2_chr(
  param_chang_tbl$Parameter,
  param_chang_tbl$Bound,
  ~ paste(.x, .y, sep = "_")
)

param_names <- readRDS(file = "param_names.rds")
param_names[merged_names] <- param_chang_tbl$years_mid_comp
saveRDS(param_names, file = "param_names.rds")

plot_univar_mod <- ggplot(
  data = param_chang_tbl,
  aes(
    y = Parameter,
    x = years_mid_comp,
    # fill = Bound
  )
) +
  geom_col(
    position = "identity",
    colour = "black"
  ) +
  # scale_fill_viridis(discrete = T) +
  labs(
    x = "Additional Life Years Gained",
    y = NULL,
    title = NULL,
    # caption = "Variation from the point estimate\nin 100 people over 10 years"
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  scale_y_discrete(labels = TeX_labels) + #or variable_labels for plain english
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


plot_univariate_25 <- ggplot(
  data = param_tbl_25 %>%
    mutate(
      Parameter = fct_reorder(.f = Parameter, .x = years_mid_comp, abssum),
      Type = factor(
        word(Parameter, sep = fixed("_")),
        levels = c("Pr", "Sens", "L", "p", "Sym")
      )
    ),
  aes(
    y = Parameter,
    x = years_mid_comp,
    # fill = Bound
  )
) +
  geom_col(
    position = "identity",
    colour = "black"
  ) +
  # scale_fill_viridis(discrete = T) +
  labs(
    x = "Additional Life Years Gained",
    y = NULL,
    title = NULL,
    # caption = "Variation from the point estimate\nin 100 people over 10 years"
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
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


# f4_0 pm 25% UPPER GRAPH ------------

plot_upper_25 <- ggplot(
  data = param_tbl_25 %>%
    filter(
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
  aes(
    y = Parameter,
    x = years_high_comp,
    # fill = Bound
  )
) +
  geom_col(
    position = "identity",
    colour = "black"
  ) +
  # scale_fill_viridis(discrete = T) +
  labs(
    x = "Additional Life Years Gained From\nIntervention per Cohort of 100",
    y = "Variable",
    title = TeX(
      r"(Upper $\textit{p}_{HCC}$ Estimated change in mortality from 25% change in each modifiable parameter)"
    ),
    # caption = "Variation from the point estimate\nin 100 people over 10 years"
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_y_discrete(labels = TeX_labels)

plot_upper_25

ggsave(
  plot = plot_upper_25,
  here("figures", "sens_anal_25_f4_0_upper.jpg"),
  width = 5.5,
  height = 9
)


# f4_0 pm 25% LOWER GRAPH ------------

plot_lower_25 <- ggplot(
  data = param_tbl_25 %>%
    filter(
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
  aes(
    x = years_low_comp,
    y = Parameter,
    # fill = Bound
  )
) +
  geom_col(
    position = "identity",
    colour = "black"
  ) +
  # scale_fill_viridis(discrete = T) +
  labs(
    x = "Additional Life Years Gained From\nIntervention per Cohort of 100",
    y = "Variable",
    title = TeX(
      r"(Lower $\textit{p}_{HCC}$ Estimated change in mortality from 25% change\nin each modifiable parameter)"
    ),
    # caption = "Variation from the point estimate\nin 100 people over 10 years"
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_y_discrete(labels = TeX_labels)

plot_lower_25

ggsave(
  plot = plot_lower_25,
  here("figures", "sens_anal_25_f4_0_lower.jpg"),
  width = 5.5,
  height = 9
)


save.image(file = "Sensitivity_analyses.RData")
