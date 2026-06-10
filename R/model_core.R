# R/model_core.R
# Core Markov model functions: progression matrix, ultrasound detection,
# adherence model, life-years calculation, and sensitivity-analysis helpers.
# All functions are pure (no file I/O, no set.seed) — they depend only on
# objects initialised in R/setup.R.

# ----- Base progression model -------------------------------------------------

#' Markov progression given n people over t months
#'
#' @param n  Starting cohort size (all in Cirrhosis at t = 0).
#' @param t  Number of months to run.
#' @param p_HCC  Annual probability of developing HCC from cirrhosis.
#' @param p_0,p_A,p_B,p_C,p_D  Annual stage-to-stage transition probabilities.
#' @return Tibble with columns matching `stage_list`.
fromf4 <- function(
  n,
  t, # Given n people, how many will progress to what after t months
  p_HCC = p_HCC_def, # NEEDS TO BE ANNUAL
  p_0   = p_0_def,
  p_A   = p_A_def,
  p_B   = p_B_def,
  p_C   = p_C_def,
  p_D   = p_D_def
) {
  p_HCC <- 1 - (1 - p_HCC)^(1 / 12)
  p_0   <- 1 - (1 - p_0)^(1 / 12)
  p_A   <- 1 - (1 - p_A)^(1 / 12)
  p_B   <- 1 - (1 - p_B)^(1 / 12)
  p_C   <- 1 - (1 - p_C)^(1 / 12)
  p_D   <- 1 - (1 - p_D)^(1 / 12)
  tmat <- t(matrix(
    c(
      1 - p_HCC, p_HCC, rep(0, 6),
      1 - p_0,   p_0,   rep(0, 6),
      1 - p_A,   p_A,   rep(0, 6),
      1 - p_B,   p_B,   rep(0, 6),
      1 - p_C,   p_C,   rep(0, 6),
      1 - p_D,   p_D,   rep(0, 6),
      1
    ),
    7, 7
  ))
  stage <- c(n, rep(0, 6)) %*% (tmat %^% t)
  tbl <- as_tibble(stage, .name_repair = ~stage_list)
  return(tbl)
}


# ----- Ultrasound detection model (by attendance probability) -----------------

#' Theoretical ultrasound detection for a cohort with a given attendance prob.
#'
#' @param n       Cohort size.
#' @param t       Total follow-up in months.
#' @param m       Ultrasound interval in months.
#' @param Pr_Att  Probability of attending a scheduled ultrasound (0 = none,
#'                1 = complete adherence).
#' @param p_HCC … p_D,Sens_*,Sym_*  Model parameters (default to *_def values).
#' @return Numeric vector of length 7 (people diagnosed per stage).
us_theor_pa <- function(
  n      = 100,
  t      = 120,
  m      = 6,
  Pr_Att, # Pr_Att is 1 for those who are comp. adh. and 0 for non
  p_HCC  = p_HCC_def, # NEEDS TO BE ANNUAL
  p_0    = p_0_def,
  p_A    = p_A_def,
  p_B    = p_B_def,
  p_C    = p_C_def,
  p_D    = p_D_def,
  Sens_0 = Sens_0_def,
  Sens_A = Sens_A_def,
  Sens_B = Sens_B_def,
  Sens_C = Sens_C_def,
  Sens_D = Sens_D_def,
  Sym_0  = Sym_0_def,
  Sym_A  = Sym_A_def,
  Sym_B  = Sym_B_def,
  Sym_C  = Sym_C_def,
  Sym_D  = Sym_D_def
) {
  # theoretical version of us_sim_pa
  p_HCC <- 1 - (1 - p_HCC)^(1 / 12)
  tmat <- t(matrix(
    c(
      1 - p_HCC, p_HCC, rep(0, 6),
      1 - p_0,   p_0,   rep(0, 6),
      1 - p_A,   p_A,   rep(0, 6),
      1 - p_B,   p_B,   rep(0, 6),
      1 - p_C,   p_C,   rep(0, 6),
      1 - p_D,   p_D,   rep(0, 6),
      1
    ),
    7, 7
  )) # Creates a transition matrix
  vec     <- c(n, rep(0, 6))
  Dx      <- rep(0, 7) # diagnosed people
  sym_prob <- c(0, Sym_0, Sym_A, Sym_B, Sym_C, Sym_D, 1)
  US_probs <- c(0, Sens_0, Sens_A, Sens_B, Sens_C, Sens_D, 1)
  for (i in 1:t) {
    vec <- vec %*% tmat
    sym <- vec * sym_prob
    Dx  <- Dx + sym
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
      us  <- go * US_probs
      Dx  <- Dx + us
      vec <- vec - us
    }
  }
  return(Dx)
}


# ----- Adherence-stratified detection model -----------------------------------

#' Combine detection across three adherence groups (non-, partial, complete).
#'
#' @param n     Cohort size.
#' @param t     Follow-up months.
#' @param m     Ultrasound interval months.
#' @param p_NA  Proportion non-adherent.
#' @param p_CA  Proportion completely adherent.
#' @param p_Att Attendance probability within partial-adherence group.
#' @param ...   Passed through to us_theor_pa.
#' @return Tibble with columns matching `stage_list`.
adh_theor <- function(
  n     = 100,
  t     = 120,
  m     = 6,
  p_NA  = p_NA_def,
  p_CA  = p_CA_def,
  p_Att = p_Att_def,
  p_HCC = p_HCC_def, # NEEDS TO BE ANNUAL
  p_0   = p_0_def,
  p_A   = p_A_def,
  p_B   = p_B_def,
  p_C   = p_C_def,
  p_D   = p_D_def,
  Sens_0 = Sens_0_def,
  Sens_A = Sens_A_def,
  Sens_B = Sens_B_def,
  Sens_C = Sens_C_def,
  Sens_D = Sens_D_def,
  Sym_0  = Sym_0_def,
  Sym_A  = Sym_A_def,
  Sym_B  = Sym_B_def,
  Sym_C  = Sym_C_def,
  Sym_D  = Sym_D_def,
  ...
) {
  # number of people, total time tracking, interval of US, uses us_theor_pa
  adh_group <- n * c(p_NA, 1 - p_NA - p_CA, p_CA) # people in each adherence group

  custom_us_theor <- function(n, Pr_Att) {
    # FUNCTION THAT JUST VARIES NUMBERS AND THEIR LIKELIHOOD TO ATTEND
    us_theor_pa(
      n, Pr_Att,
      t = t, m = m,
      p_HCC = p_HCC,
      p_0 = p_0, p_A = p_A, p_B = p_B, p_C = p_C, p_D = p_D,
      Sens_0 = Sens_0, Sens_A = Sens_A, Sens_B = Sens_B,
      Sens_C = Sens_C, Sens_D = Sens_D,
      Sym_0 = Sym_0, Sym_A = Sym_A, Sym_B = Sym_B,
      Sym_C = Sym_C, Sym_D = Sym_D
    ) %>%
      as_tibble(.name_repair = ~stage_list)
  }

  stage_caught <- map2_dfr(
    adh_group,
    c(0, p_Att, 1),
    .f = custom_us_theor
  ) %>%
    summarise_all(sum)

  return(stage_caught)
}


# ----- Life-years added by surveillance and treatment -------------------------

#' Total life years added by the HCC surveillance programme.
#'
#' Sums the product of stage-specific diagnoses and stage-specific life-year
#' gains (L_*) across all adherence groups.
#'
#' @param n,t,m,p_NA,p_CA,p_Att,p_HCC,p_*,Sens_*,Sym_*,L_*
#'   All default to *_def point estimates.
#' @return Scalar life-years gained per cohort of n.
life_add <- function(
  n     = 100,
  t     = 120,
  m     = 6,
  p_NA  = p_NA_def,
  p_CA  = p_CA_def,
  p_Att = p_Att_def,
  p_HCC = p_HCC_def, # NEEDS TO BE ANNUAL
  p_0   = p_0_def,
  p_A   = p_A_def,
  p_B   = p_B_def,
  p_C   = p_C_def,
  p_D   = p_D_def,
  Sens_0 = Sens_0_def,
  Sens_A = Sens_A_def,
  Sens_B = Sens_B_def,
  Sens_C = Sens_C_def,
  Sens_D = Sens_D_def,
  Sym_0  = Sym_0_def,
  Sym_A  = Sym_A_def,
  Sym_B  = Sym_B_def,
  Sym_C  = Sym_C_def,
  Sym_D  = Sym_D_def,
  L_0    = L_0_def,
  L_A    = L_A_def,
  L_B    = L_B_def,
  L_C    = L_C_def,
  L_D    = L_D_def,
  ...
) {
  stage_caught <- adh_theor(
    n = n, m = m, t = t,
    p_NA = p_NA, p_CA = p_CA, p_Att = p_Att,
    p_HCC = p_HCC,
    p_0 = p_0, p_A = p_A, p_B = p_B, p_C = p_C, p_D = p_D,
    Sens_0 = Sens_0, Sens_A = Sens_A, Sens_B = Sens_B,
    Sens_C = Sens_C, Sens_D = Sens_D,
    Sym_0 = Sym_0, Sym_A = Sym_A, Sym_B = Sym_B,
    Sym_C = Sym_C, Sym_D = Sym_D
  )
  life_added <- sum(stage_caught * c(0, L_0, L_A, L_B, L_C, L_D, 0))
  return(life_added)
}


# ----- Sensitivity-analysis helpers -------------------------------------------

#' Build a one-row tibble of named parameter values for use with pmap.
dfs <- function(Parameter, value, ...) {
  new_df <- tibble(Parameter, value) %>%
    pivot_wider(names_from = Parameter, values_from = value)
  return(new_df)
}

#' Apply life_add across each row of a data-frame of parameter draws.
years <- function(df) pmap_dbl(df, life_add)

#' Sum of absolute values.
abssum <- function(x) sum(abs(x))

#' Build the parameter tibble used for univariate sensitivity tornado plots.
#'
#' @param csv  Full path to a Parameters CSV (e.g. here("data","Parameters.csv")).
#' @return Tibble with one row per parameter × bound combination plus computed
#'   life-year columns.
life_tibble <- function(csv) {
  param_csv <- read.csv(file = csv, row.names = 1)
  param_tidy <- param_csv[1:24, 2:3] %>%
    as.data.frame() %>%
    rownames_to_column(var = "Parameter") %>%
    as_tibble() %>%
    pivot_longer(-Parameter, names_to = "Bound")
  param_tbl <- param_tidy %>%
    mutate(
      dfs_new  = pmap(., dfs),
      Bound    = as.factor(str_to_title(Bound)),
      Parameter = as.factor(Parameter),
      dfs_low  = map(
        dfs_new,
        ~ if (names(.x) == "p_HCC") .x else
          tibble(.x, p_HCC = param_csv["p_HCC", "LOWER"])
      ),
      dfs_high = map(
        dfs_new,
        ~ if (names(.x) == "p_HCC") .x else
          tibble(.x, p_HCC = param_csv["p_HCC", "UPPER"])
      ),
      years_mid        = map_dbl(dfs_new, years),
      years_low        = map_dbl(dfs_low, years),
      years_high       = map_dbl(dfs_high, years),
      years_mid_comp   = years_mid  - life_add(),
      years_low_comp   = years_low  - life_add(p_HCC = param_csv["p_HCC", "LOWER"]),
      years_high_comp  = years_high - life_add(p_HCC = param_csv["p_HCC", "UPPER"])
    )
  return(param_tbl)
}
