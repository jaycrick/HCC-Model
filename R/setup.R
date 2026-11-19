# R/setup.R
# Shared initialisation: loaded at the top of every analysis script.
# Defines libraries, conflict preferences, global constants, default
# parameters, and shared visual objects.

# ----- Libraries --------------------------------------------------------------
library(conflicted) # must load first to intercept ambiguous calls
library(expm) # matrix exponentiation (%^%)
library(tidyverse)
library(here)
library(viridis) # turbo palette, scale_fill_viridis
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

#' Read a Parameters CSV into a named list.
#'
#' Values are matched by the PARAMETER column, never by position, and the
#' UTF-8 byte-order mark in the CSV header is handled.
#'
#' @param csv     Path to a Parameters CSV (columns PARAMETER, POINT, LOWER,
#'                UPPER).
#' @param column  Which value column to read ("POINT", "LOWER" or "UPPER").
#' @return Named list: one numeric element per PARAMETER row.
load_params <- function(csv, column = "POINT") {
  tbl <- read.csv(file = csv, fileEncoding = "UTF-8-BOM")
  stopifnot(
    "PARAMETER" %in% names(tbl),
    column %in% names(tbl),
    !anyDuplicated(tbl$PARAMETER)
  )
  as.list(setNames(tbl[[column]], tbl$PARAMETER))
}

# Named list of POINT estimates by PARAMETER name, e.g. params_def$p_HCC.
params_def <- load_params(here("data", "Parameters.csv"))

# ----- Line types, colour palette and ggplot scales ---------------------------
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

#' Build one LaTeX label expression, e.g. \\textit{p}_{CA}.
#'
#' @param base    Symbol before the subscript (e.g. "p", "L", "Sym").
#' @param sub     Subscript text (LaTeX, e.g. "CA" or "\\alpha").
#' @param italic  Wrap `base` in \\textit{}? (p_*, L_*, m, D are italic.)
#' @param star    Add the bootstrap-replicate superscript `*`?
#' @return A plotmath expression from latex2exp::TeX().
tex_label <- function(base, sub, italic = TRUE, star = FALSE) {
  b <- if (italic) paste0("\\textit{", base, "}") else base
  sup <- if (star) "^{*}" else ""
  TeX(paste0("$", b, sup, "_{", sub, "}$"))
}

# TeX_labels: parameter names as LaTeX expressions (used in multiple scripts)
tex_param_names <- names(params_def)
TeX_labels <- lapply(tex_param_names, function(nm) {
  parts <- strsplit(nm, "_", fixed = TRUE)[[1]]
  if (nm == "p_NA") {
    # `NA` needs backticks in plotmath, so TeX() cannot be used here
    return(expression(italic(p)[`NA`]))
  }
  if (nm == "p_B") {
    # G2: preserved bug, p_B renders as p_C
    parts <- c("p", "C")
  }
  tex_label(parts[1], parts[2], italic = parts[1] %in% c("p", "L"))
})
TeX_labels <- do.call(c, TeX_labels) # expression vector, as before
names(TeX_labels) <- tex_param_names

# Boot_labels: bootstrap parameter names as LaTeX expressions
boot_specs <- list(
  m = c("0", "A", "B", "C", "D"),
  D = c("0", "A", "alpha"),
  p = c("0", "A", "B", "C", "D")
)
Boot_labels <- do.call(
  c,
  unlist(
    lapply(names(boot_specs), function(base) {
      lapply(
        setNames(boot_specs[[base]], paste0(base, "_", boot_specs[[base]])),
        function(sub) {
          tex_label(base, if (sub == "alpha") "\\alpha" else sub, star = TRUE)
        }
      )
    }),
    recursive = FALSE
  )
)
