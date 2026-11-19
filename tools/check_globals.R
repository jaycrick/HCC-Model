#!/usr/bin/env Rscript
# Undefined-name check for the project's own functions.
#
# lintr's object_usage_linter lints each file in isolation, so in this sourced
# layout (library() calls in R/setup.R, helpers in other R/ files) it reports
# almost only false positives. This script instead sources the project exactly
# as tests/testthat/helper-setup.R does, then runs codetools over the project
# functions only (never package functions) and reports every name that is
# neither defined in the project, nor in an attached package, nor local.
#
# Usage: Rscript tools/check_globals.R      exit status 1 if anything is found
#
# Tidy-evaluation column names (`mutate(x = col + 1)`) and a few historic
# globals look undefined to codetools. They are listed per function in
# `allowed` below, each with its reason. Anything not listed fails.

suppressPackageStartupMessages(library(codetools))

# ----- Load the project in the helper-setup.R order ---------------------------
r_files <- list.files("R", pattern = "\\.R$")
first <- c("setup.R", "model_core.R", "bootstrap.R", "scenarios.R")
load_order <- c(first, setdiff(r_files, first))
before <- ls(globalenv(), all.names = TRUE)
for (f in load_order) {
  source(file.path("R", f), local = FALSE)
}
project_fns <- setdiff(ls(globalenv(), all.names = TRUE), before)
project_fns <- project_fns[vapply(
  project_fns,
  function(n) is.function(get(n, globalenv())),
  logical(1)
)]

# ----- Allow-list: known false positives and tracked debt ---------------------
# Remove an entry as soon as its cause is fixed; stale entries fail the check.
allowed <- list(
  # tidy-eval column names inside dplyr/tidyr verbs
  life_tibble = c(
    ".",
    "Parameter",
    "Bound",
    "dfs_new",
    "dfs_low",
    "dfs_high",
    "years_mid",
    "years_low",
    "years_high"
  ),
  # tidy-eval / aesthetic column names in the analysis/01 figure builders
  plot_baseline_hist = "Baseline",
  plot_bootstrap_counts = c("name", "value"),
  plot_rates_hist = c("name", "Point", "rate", "value"),
  plot_rates_distribution = c("name", "Point", "rate", "value"),
  plot_pcor_baseline = c("value", "Variable"),
  plot_pcor_heatmap = "value",
  plot_sensitivity_hex = c("Sens_0", "Sens_A"),
  plot_time_range = c(
    "mean_years",
    "mean025",
    "mean975",
    "name",
    "Point",
    "rate",
    "value"
  ),
  bootstrap_estimates = "list_of_pars",
  # tracked debt, finding C10 (R/simulation.R depends on undefined globals;
  # fixed in the simulation step of the plan)
  adh_sim = c("pr_na", "pr_pa"),
  hccsim = "tmat",
  us_sim = c(
    "tmat",
    paste0("Sym_", c("0", "A", "B", "C", "D")),
    paste0("Sens_", c("0", "A", "B", "C", "D"))
  ),
  us_sim_pa = c(
    "tmat",
    "pr_att",
    paste0("Sym_", c("0", "A", "B", "C", "D")),
    paste0("Sens_", c("0", "A", "B", "C", "D"))
  )
)

# ----- Run codetools over the project functions -------------------------------
collector <- new.env()
collector$found <- character(0)
collect <- function(msg) {
  collector$found <- c(collector$found, trimws(msg))
}
for (fn in project_fns) {
  codetools::checkUsage(
    get(fn, globalenv()),
    name = fn,
    report = collect,
    suppressLocal = TRUE
  )
}

# Messages look like "pafunc : tempfun: no visible binding for global variable
# 'x'" (nested) or "pafunc: no visible global function definition for 'f'".
pattern <- paste0(
  "^([^ :]+).*: no visible ",
  "(binding for global variable|global function definition for) .(.+).$"
)
found <- collector$found
hits <- found[grepl("no visible", found)]
fn_of <- sub(pattern, "\\1", hits)
name_of <- sub(pattern, "\\3", hits)
what_of <- ifelse(
  grepl("global function", hits),
  "undefined function",
  "undefined variable"
)
is_allowed <- mapply(
  function(fn, nm) nm %in% allowed[[fn]],
  fn_of,
  name_of
)
bad <- unique(data.frame(
  fn = fn_of[!is_allowed],
  name = name_of[!is_allowed],
  what = what_of[!is_allowed]
))

# Stale allow-list entries (cause fixed but entry kept) are also an error.
stale <- unlist(lapply(names(allowed), function(fn) {
  gone <- setdiff(allowed[[fn]], name_of[fn_of == fn])
  if (length(gone)) paste0(fn, ": ", gone)
}))

if (nrow(bad) > 0) {
  cat("check_globals: undefined names in project functions:\n")
  cat(sprintf("  %s: %s '%s'\n", bad$fn, bad$what, bad$name), sep = "")
}
if (length(stale) > 0) {
  cat("check_globals: stale allow-list entries (remove them):\n")
  cat(sprintf("  %s\n", stale), sep = "")
}
if (nrow(bad) > 0 || length(stale) > 0) {
  quit(status = 1)
}
cat(sprintf(
  "check_globals: OK (%d project functions, no undefined names)\n",
  length(project_fns)
))
