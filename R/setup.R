# Shared setup (file R/setup.R)
# Shared initialisation: loaded at the top of every analysis script.
# Defines libraries, conflict preferences, global constants, default
# parameters, and shared visual objects.
# Parameter values: Cumming 2024, Tables 2-3; their sources are listed in
# the file Parameters_sources.csv (in papers/2024-cancers/data/).

# ----- Libraries --------------------------------------------------------------
library(conflicted) # must load first to intercept ambiguous calls
library(expm) # matrix exponentiation (%^%)
library(dplyr) # data manipulation (filter, mutate, summarise, ...)
library(tidyr) # pivot_longer
library(tibble) # tibble, as_tibble
library(ggplot2) # plotting
library(purrr) # map, pmap
library(forcats) # fct_inorder
library(stringr) # str_to_title
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

# ----- Paper directories ------------------------------------------------------
# R/ is the model core shared by every paper. Each paper keeps its own inputs
# and outputs in papers/<id>/ (data/, data-derived/, output/, figures/), and
# the core reads and writes them only through paper_path(). The paper is
# chosen with options(hcc.paper = "<id>"), which each paper's scripts set
# before sourcing this file. The default is the published model, which is
# also the regression reference for the test suite.
DEFAULT_PAPER <- "2024-cancers"

#' Path inside the active paper's directory.
#'
#' Ref: none (repository layout).
#'
#' @param ...    Path components below papers/<paper>/.
#' @param paper  Paper id (directory name under papers/).
#' @return An absolute path.
paper_path <- function(..., paper = getOption("hcc.paper", DEFAULT_PAPER)) {
  here::here("papers", paper, ...)
}

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

# Stage suffixes of the per-stage parameters (BCLC 0, A, B, C, D), and the
# parameter names built from them, in the order of the transition matrix rows.
STAGE_IDS <- c("0", "A", "B", "C", "D")
STAGE_P_NAMES <- paste0("p_", STAGE_IDS) # monthly progression probabilities
STAGE_SYM_NAMES <- paste0("Sym_", STAGE_IDS) # symptomatic presentation
STAGE_SENS_NAMES <- paste0("Sens_", STAGE_IDS) # ultrasound sensitivity
STAGE_L_NAMES <- paste0("L_", STAGE_IDS) # life expectancy at diagnosis

# ----- Default parameters (point estimates) -----------------------------------

#' Read a Parameters CSV into a named list.
#'
#' Values are matched by the PARAMETER column, never by position, and the
#' UTF-8 byte-order mark in the CSV header is handled.
#'
#' Ref: Cumming 2024, Tables 2 and 3 (values); sources per parameter in
#' data/Parameters_sources.csv.
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
params_def <- load_params(paper_path("data", "Parameters.csv"))

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
#' Ref: none (plot labels in paper notation, Tables 2-3).
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
  tex_label(parts[1], parts[2], italic = parts[1] %in% c("p", "L"))
})
TeX_labels <- do.call(c, TeX_labels) # expression vector, as before
names(TeX_labels) <- tex_param_names

# Boot_labels: bootstrap parameter names as LaTeX expressions
boot_specs <- list(
  m = STAGE_IDS,
  D = c("0", "A", "alpha"),
  p = STAGE_IDS
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

# ----- Parameter groups and labels (one lookup) -------------------------------

# The three adherence parameters (non-attenders, completely adherent, partial
# attendance probability).
ADHERENCE_PARAMS <- c("p_NA", "p_CA", "p_Att")

# Facet groups, in plotting order, with their facet labels. Adherence is its own
# group (p_NA, p_CA, p_Att), separate from the monthly transition probabilities.
group_labels <- c(
  Adherence = "Adherence\nto Screening",
  Sens = "Ultrasound\nSensitivity",
  L = "Treatment Quality",
  p = "Monthly Risk",
  Sym = "Monthly probability\nof symptoms"
)

#' Group of each parameter: the adherence parameters, else the name prefix
#' (Sens, L, p = transition probabilities incl. p_HCC, Sym).
#'
#' Ref: none (plot grouping; groups follow the Fig. 2D panels).
#'
#' @param nm  Parameter names.
#' @return Character vector of group names (names of `group_labels`).
param_group <- function(nm) {
  ifelse(
    nm %in% ADHERENCE_PARAMS,
    "Adherence",
    vapply(strsplit(nm, "_", fixed = TRUE), `[`, "", 1)
  )
}

# Axis labels of the parameters inside their group.
.variable_labels <- c(
  L_0 = "Stage 0",
  L_A = "Stage A",
  L_B = "Stage B",
  L_C = "Stage C",
  L_D = "Stage D",
  p_CA = "Completely adherent",
  p_NA = "Non-adherent",
  p_Att = "Attendence probability\nwithin partial adherents",
  p_HCC = "Developing HCC",
  p_0 = "Stage 0 to A",
  p_A = "Stage A to B HCC",
  p_B = "Stage B to C HCC",
  p_C = "Stage C to D HCC",
  p_D = "Stage D to Death",
  Sym_0 = "Stage 0",
  Sym_A = "Stage A",
  Sym_B = "Stage B",
  Sym_C = "Stage C",
  Sym_D = "Stage D",
  Sens_0 = "Stage 0",
  Sens_A = "Stage A",
  Sens_B = "Stage B",
  Sens_C = "Stage C",
  Sens_D = "Stage D"
)

# The single lookup: one row per parameter of `params_def`.
param_lookup <- tibble::tibble(
  parameter = names(params_def),
  group = param_group(names(params_def)),
  label = unname(.variable_labels[names(params_def)])
)

# Derived, never edited separately: axis labels per parameter, and facet labels
# (parameter groups, plus the low/high-HCC-rate columns of the sensitivity
# grid).
variable_labels <- setNames(param_lookup$label, param_lookup$parameter)
facet_labels <- c(
  group_labels,
  years_low_comp = "Lower Rate of HCC",
  years_high_comp = "Higher Rate of HCC"
)

# Parameters that no intervention can change (disease progression, symptom
# probabilities, HCC incidence): excluded from the "modifiable" sensitivity
# plots.
NON_MODIFIABLE_PARAMS <- c("p_HCC", STAGE_P_NAMES, STAGE_SYM_NAMES)

# ----- Manuscript-values ledger -----------------------------------------------
# Every number quoted in the manuscript is recorded here, in one place, by
# record_values(). The ledger is output/manuscript_values.csv (name, value),
# sorted by name and written at a precision that round-trips exactly, so it
# is text-diffable. data-derived/ is read-only input: nothing writes to it
# except the deliberate `--bootstrap` regeneration in analysis/01.

# Shortest decimal string (15-17 significant digits) that R's own parser reads
# back to exactly the same double. (R's parser is not always correctly rounded,
# so "%.17g" alone does not guarantee an exact round trip.)
#' @param x  A double.
#' @return A string that parses back to exactly `x`.
round_trip_string <- function(x) {
  for (digits in 15:17) {
    s <- sprintf(paste0("%.", digits, "g"), x)
    if (identical(as.numeric(s), x)) return(s)
  }
  sprintf("%.17g", x)
}

#' Merge named scalar values into the manuscript-values ledger.
#'
#' Re-recording a name overwrites it. Rows are sorted by name (C locale, so the
#' order does not depend on the user's locale).
#'
#' Ref: none (bookkeeping for the numbers quoted in the manuscript).
#'
#' @param values Named numeric vector or list of scalars.
#' @param path Ledger CSV path.
#' @return The full ledger as a named numeric vector, invisibly.
record_values <- function(
  values,
  path = paper_path("output", "manuscript_values.csv")
) {
  values <- unlist(values)
  stopifnot(
    is.numeric(values),
    !is.null(names(values)),
    !anyNA(names(values)),
    all(nzchar(names(values))),
    !anyDuplicated(names(values))
  )
  ledger <- if (file.exists(path)) {
    old <- utils::read.csv(path, stringsAsFactors = FALSE)
    stats::setNames(as.numeric(old$value), old$name)
  } else {
    numeric(0)
  }
  ledger[names(values)] <- values
  ledger <- ledger[order(names(ledger), method = "radix")]
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  out <- data.frame(
    name = names(ledger),
    value = vapply(ledger, round_trip_string, character(1), USE.NAMES = FALSE)
  )
  utils::write.csv(out, path, row.names = FALSE, quote = FALSE)
  invisible(ledger)
}
