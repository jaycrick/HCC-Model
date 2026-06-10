# R/simulation.R
# Stochastic (Monte Carlo) simulation functions.
# NOTE: These are exploratory functions not used in the published main analysis
# (which uses the theoretical model in R/model_core.R).  They reference a
# global `tmat` transition matrix that must be constructed by the caller; see
# analysis/03_progression_model.R for context.

# ----- Basic stochastic progression ------------------------------------------

#' Stochastically progress n people over m months (no screening).
#'
#' @param n  Cohort size.
#' @param m  Number of months.
#' @return Numeric vector of length 7: counts per stage.
hccsim <- function(n, m) {
  # n = number of people, m = number of months
  vec <- c(n, rep(0, 6)) # initialise with all people in cirrhosis (F4)
  for (i in 1:m) {
    for (j in (length(vec) - 1):1) {
      worse <- rbinom(1, vec[j], tmat[j, j + 1]) # people who progress
      vec[j + 1] <- vec[j + 1] + worse
      vec[j] <- vec[j] - worse
    }
  }
  return(vec)
}

#' Compare two screening intervals using hccsim.
simcomp <- function(n, l, s) hccsim(n, l) - hccsim(n, s)


# ----- Stochastic progression with ultrasound screening ----------------------

#' Stochastic model with 6-monthly (or m-monthly) ultrasound.
#'
#' References globals: tmat, Sym_0–Sym_D, Sens_0–Sens_D.
#'
#' @param n  Cohort size.
#' @param t  Total follow-up months.
#' @param m  Ultrasound interval months.
#' @return Numeric vector of length 7.
us_sim <- function(n, t, m) {
  vec <- c(n, rep(0, 6))
  fin <- rep(0, 7)
  for (i in 1:t) {
    for (j in (length(vec) - 1):1) {
      worse <- rbinom(1, vec[j], tmat[j, j + 1])
      vec[j + 1] <- vec[j + 1] + worse
      vec[j] <- vec[j] - worse
    }
    sym <- c(
      0,
      rbinom(1, vec[2], Sym_0),
      rbinom(1, vec[3], Sym_A),
      rbinom(1, vec[4], Sym_B),
      rbinom(1, vec[5], Sym_C),
      rbinom(1, vec[6], Sym_D),
      vec[7]
    )
    fin <- fin + sym
    vec <- vec - sym
    if (i %% m == 0) {
      us <- c(
        0,
        rbinom(1, vec[2], Sens_0),
        rbinom(1, vec[3], Sens_A),
        rbinom(1, vec[4], Sens_B),
        rbinom(1, vec[5], Sens_C),
        rbinom(1, vec[6], Sens_D),
        vec[7]
      ) # how many are picked up at the ultrasound
      fin <- fin + us
      vec <- vec - us
    }
  }
  return(fin)
}


# ----- Stochastic model for the partially-adherent group ---------------------

#' Stochastic model for the partial-adherence cohort (pr_att attendance rate).
#'
#' References globals: tmat, pr_att, Sym_0–Sym_D, Sens_0–Sens_D.
#'
#' @param n  Cohort size.
#' @param t  Total follow-up months.
#' @param m  Ultrasound interval months.
#' @return Numeric vector of length 7.
us_sim_pa <- function(n, t, m) {
  vec <- c(n, rep(0, 6))
  fin <- rep(0, 7)
  for (i in 1:t) {
    for (j in (length(vec) - 1):1) {
      worse <- rbinom(1, vec[j], tmat[j, j + 1])
      vec[j + 1] <- vec[j + 1] + worse
      vec[j] <- vec[j] - worse
      if (i %% m == 0) {
        sym <- c(
          0,
          rbinom(1, vec[2], Sym_0),
          rbinom(1, vec[3], Sym_A),
          rbinom(1, vec[4], Sym_B),
          rbinom(1, vec[5], Sym_C),
          rbinom(1, vec[6], Sym_D),
          vec[7]
        )
        fin <- fin + sym
        vec <- vec - sym
        go <- c(
          0,
          rbinom(1, vec[2], pr_att),
          rbinom(1, vec[3], pr_att),
          rbinom(1, vec[4], pr_att),
          rbinom(1, vec[5], pr_att),
          rbinom(1, vec[6], pr_att),
          vec[7]
        ) # how many actually go to the ultrasound
        us <- c(
          0,
          rbinom(1, go[2], Sens_0),
          rbinom(1, go[3], Sens_A),
          rbinom(1, go[4], Sens_B),
          rbinom(1, go[5], Sens_C),
          rbinom(1, go[6], Sens_D),
          vec[7]
        )
        fin <- fin + us
        vec <- vec - us
      }
    }
  }
  return(fin)
}


# ----- Adherence-stratified stochastic model ---------------------------------

#' Stochastic adherence-stratified model.
#'
#' References globals: pr_na, pr_pa (adherence proportions) plus those used
#' by us_sim and us_sim_pa.
#'
#' @param n  Cohort size.
#' @param t  Total follow-up months.
#' @param m  Ultrasound interval months.
#' @return Numeric vector of length 7.
adh_sim <- function(n, t, m) {
  adh_group <- rmultinom(1, n, c(pr_na, pr_pa, 1 - pr_na - pr_pa))
  stage_caught <- us_sim(adh_group[3], t, m) + us_sim_pa(adh_group[2], t, m)
  stage_caught[7] <- stage_caught[7] + adh_group[1]
  return(stage_caught)
}
