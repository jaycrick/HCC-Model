# R/bootstrap.R
# Functions used in the bootstrap parameter estimation and uncertainty analysis.
# Requires: R/setup.R, R/model_core.R, and the hypogeom package.

# ----- Negative log-likelihood for parametric bootstrap ----------------------

#' Negative log-likelihood for joint parameter estimation.
#'
#' Combines hypergeometric-order-statistics likelihoods for median survival
#' times with binomial likelihoods for diagnosis proportions.
#'
#' @param hypo_probs  Length-8 numeric vector of candidate probabilities.
#' @param meds        Observed median survival times per stage.
#' @param Dx_0,Dx_A,Dx_sym_early  Observed diagnosis counts.
#' @return Scalar NLL, or NaN if any probability is outside [0,1].
neg_log_lik <- function(
  hypo_probs,
  meds         = c(med_0, med_a, med_b, med_c, med_d),
  Dx_0         = 36,
  Dx_A         = 51,
  Dx_sym_early = 52
) {
  while (prod((0 <= hypo_probs) * (hypo_probs <= 1)) == 0) {
    return(NaN)
  }

  p_0a   <- hypo_probs[1]
  p_Ab   <- hypo_probs[2]
  p_Bc   <- hypo_probs[3]
  US_00  <- hypo_probs[6]
  US_Aa  <- hypo_probs[7]
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
      (1 - `p_0a`) * US_00, `p_0a` * US_Aa,
      (1 - `p_0a`) * (1 - US_00), rep(0, 5),
      `p_0a` * (1 - US_Aa), rep(0, 6),
      rep(0, 9), 1 - p_Ab, rep(0, 4), p_Ab,
      rep(0, 10), 1 - p_Ab, rep(0, 3), p_Ab,
      rep(0, 11), 1 - p_Ab, rep(0, 2), p_Ab,
      rep(0, 12), 1 - p_Ab, 0, p_Ab,
      rep(0, 13), 1 - p_Ab, p_Ab,
      0, (1 - p_Ab) * US_Aa, rep(0, 6),
      (1 - p_Ab) * (1 - US_Aa), rep(0, 5), p_Ab,
      rep(0, 14), 1
    ),
    byrow = TRUE,
    nrow  = 15
  )

  end_vector <- c(
    US_00 / 6, 0, (1 - US_00) / 6, rep(1 / 6, 5), rep(0, 7)
  ) %*%
    (matrix_of_trans %^% 100)

  NLL <- -(
    dOrder(meds[5] - 1, distn = "hypogeom", probs = hypo_probs[5],   mlen = 209, j = 210 / 2, log = TRUE) +
    dOrder(meds[4] - 2, distn = "hypogeom", probs = hypo_probs[4:5], mlen = 137, j = 138 / 2, log = TRUE) +
    dOrder(meds[3] - 3, distn = "hypogeom", probs = hypo_probs[3:5], mlen = 139, j = 140 / 2, log = TRUE) +
    dOrder(meds[2] - 4, distn = "hypogeom", probs = hypo_probs[2:5], mlen = 101, j = 102 / 2, log = TRUE) +
    dOrder(meds[1] - 5, distn = "hypogeom", probs = hypo_probs[1:5], mlen = 11,  j = 12  / 2, log = TRUE) +
    dbinom(Dx_0, 109, end_vector[1], log = TRUE) +
    dbinom(Dx_A, 109 - Dx_0, end_vector[2] / (1 - end_vector[1]), log = TRUE) +
    dbinom(
      Dx_sym_early, 158,
      Sym_Aa / (1 - (1 - p_Ab) * (1 - Sym_Aa)) +
        (1 - Sym_Aa / (1 - (1 - p_Ab) * (1 - Sym_Aa))) *
          Sym_Bb / (1 - (1 - p_Bc) * (1 - Sym_Bb)),
      log = TRUE
    )
  )
  return(NLL)
}


# ----- Scenario functions for bootstrap uncertainty analysis -----------------
# Each function wraps life_add with the parameter set from one scenario,
# accepting columns from the bootstrapped resampled_pars tibble via pmap.

#' Baseline: current practice.
base_scen_fun <- function(
  `p_0`, p_A, p_B, p_C, p_D, Sens_0, Sens_A,
  Sym_A = Sym_A_def, Sym_B = Sym_B_def,
  p_HCC = p_HCC_def, ...
) {
  life_add(
    p_HCC = p_HCC,
    `p_0` = `p_0`, p_A = p_A, p_B = p_B, p_C = p_C, p_D = p_D,
    Sens_0 = Sens_0, Sens_A = Sens_A,
    Sym_A = Sym_A, Sym_B = Sym_A * 8
  )
}

#' Scenario 1a: realistic adherence improvement.
scen_1a_fun <- function(
  `p_0`, p_A, p_B, p_C, p_D, Sens_0, Sens_A,
  Sym_A = Sym_A_def, Sym_B = Sym_B_def,
  p_HCC = p_HCC_def, ...
) {
  life_add(
    p_HCC = p_HCC,
    p_NA  = p_NA_def  - 0.05,
    p_CA  = p_CA_def  + 0.10,
    p_Att = p_Att_def + 0.10,
    `p_0` = `p_0`, p_A = p_A, p_B = p_B, p_C = p_C, p_D = p_D,
    Sens_0 = Sens_0, Sens_A = Sens_A,
    Sym_A = Sym_A, Sym_B = Sym_A * 8
  )
}

#' Scenario 1b: optimistic adherence improvement.
scen_1b_fun <- function(
  `p_0`, p_A, p_B, p_C, p_D, Sens_0, Sens_A,
  Sym_A = Sym_A_def, Sym_B = Sym_B_def,
  p_HCC = p_HCC_def, ...
) {
  life_add(
    p_HCC = p_HCC,
    p_NA  = p_NA_def  - 0.10,
    p_CA  = p_CA_def  + 0.15,
    p_Att = p_Att_def + 0.20,
    `p_0` = `p_0`, p_A = p_A, p_B = p_B, p_C = p_C, p_D = p_D,
    Sens_0 = Sens_0, Sens_A = Sens_A,
    Sym_A = Sym_A, Sym_B = Sym_A * 8
  )
}

#' Scenario 2: improved ultrasound sensitivity.
scen_2_fun <- function(
  `p_0`, p_A, p_B, p_C, p_D, Sens_0, Sens_A,
  Sym_A = Sym_A_def, Sym_B = Sym_B_def,
  p_HCC = p_HCC_def, ...
) {
  life_add(
    p_HCC = p_HCC,
    `p_0`  = `p_0`, p_A = p_A, p_B = p_B, p_C = p_C, p_D = p_D,
    Sens_0 = Sens_0 + 0.05,
    Sens_A = Sens_A + 0.05,
    Sens_B = Sens_B_def + 0.05,
    Sym_A  = Sym_A, Sym_B = Sym_A * 8
  )
}

#' Scenario 3: better treatments.
scen_3_fun <- function(
  `p_0`, p_A, p_B, p_C, p_D, Sens_0, Sens_A,
  Sym_A = Sym_A_def, Sym_B = Sym_B_def,
  p_HCC = p_HCC_def, ...
) {
  life_add(
    p_HCC = p_HCC,
    L_A   = L_A_def + 1,
    L_B   = L_B_def + 0.5,
    L_C   = L_C_def + 0.5,
    `p_0` = `p_0`, p_A = p_A, p_B = p_B, p_C = p_C, p_D = p_D,
    Sens_0 = Sens_0, Sens_A = Sens_A,
    Sym_A  = Sym_A, Sym_B = Sym_A * 8
  )
}


# ----- LaTeX label helper ----------------------------------------------------

#' Convert a variable name to its TeX expression label.
Tex_labeller <- function(value) {
  return(TeX_labels[value])
}
