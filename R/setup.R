# R/setup.R
# Shared initialisation: loaded at the top of every analysis script.
# Defines libraries, conflict preferences, global constants, default
# parameters, and shared visual objects.

# ----- Libraries --------------------------------------------------------------
library(conflicted) # must load first to intercept ambiguous calls
library(cli)
library(expm) # matrix exponentiation (%^%)
library(tidyverse)
library(here)
library(viridis)
library(scales)
library(latex2exp)

# ----- Conflict preferences ---------------------------------------------------
conflicts_prefer(
  dplyr::filter,
  dplyr::select,
  dplyr::lag
)

# ----- Reproducibility --------------------------------------------------------
set.seed(4567385)

# ----- Stage labels -----------------------------------------------------------
stage_list <- c(
  "Cirrhosis",
  "BCLC 0",
  "BCLC A",
  "BCLC B",
  "BCLC C",
  "BCLC D",
  "Mortality"
)

# ----- Default parameters (point estimates) -----------------------------------
# Reads Parameters.csv from data/ and assigns every row as  <name>_def
param_def <- read.csv(
  file = here("data", "Parameters.csv"),
  row.names = 1
)[1:3]

for (i in row.names(param_def)) {
  assign(paste(i, "def", sep = "_"), param_def[i, "POINT"])
}

# ----- Colour palettes and ggplot scales --------------------------------------
library(ggsci) # pal_igv
library(viridis) # turbo

my_pal <- pal_igv("default")(7)
names(my_pal) <- c(
  "Cirrhosis",
  "Mortality",
  "BCLC 0",
  "BCLC A",
  "BCLC D",
  "BCLC B",
  "BCLC C"
)

my_linet <- as.vector(c(
  "solid",
  "solid",
  "solid",
  "solid",
  "solid",
  "solid",
  "dotdash"
))
names(my_linet) <- stage_list

# Use turbo palette (overrides IGV-based scales)
turbo_pal <- turbo(7)
names(turbo_pal) <- stage_list

my_scale <- scale_colour_manual(
  name = "HCC_Stage",
  values = turbo_pal,
  limits = stage_list
)
my_fill <- scale_fill_manual(
  name = "HCC_Stage",
  values = turbo_pal,
  limits = stage_list
)
my_line <- scale_linetype_manual(
  name = "HCC_Stage",
  values = my_linet,
  limits = stage_list
)

# ----- Shared axis / legend labels --------------------------------------------
# TeX_labels: parameter names as LaTeX expressions (used in multiple scripts)
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

# Boot_labels: bootstrap parameter names as LaTeX expressions
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
