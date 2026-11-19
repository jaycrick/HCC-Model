# ONE-OFF generator for tests/testthat/fixtures/ledger_baseline.rds.
#
# The baseline that test-ledger.R compares the 02 and 04 ledger values against:
# a named numeric vector of the manuscript-values ledger
# (output/manuscript_values.csv) for the names that also appear in the frozen
# data-derived/param_names.rds (the published reference, never rewritten).
#
# RULE: regenerate only under an approved gate (G1..G9) and record the gate in
# the commit message. Never regenerate to make a failing test pass.
# Current content: G9 (refit; also G8), from `Rscript run_all.R --bootstrap`.
#
# Usage (from repo root), after run_all.R --bootstrap:
#   Rscript tests/testthat/fixtures/make_ledger_baseline.R

frozen <- readRDS(here::here("data-derived", "param_names.rds"))
ledger <- utils::read.csv(here::here("output", "manuscript_values.csv"))
ledger <- stats::setNames(as.numeric(ledger$value), ledger$name)

keep <- intersect(names(frozen), names(ledger))
baseline <- ledger[keep]
saveRDS(
  baseline,
  here::here("tests", "testthat", "fixtures", "ledger_baseline.rds")
)
message("wrote ledger_baseline.rds (", length(baseline), " names)")
