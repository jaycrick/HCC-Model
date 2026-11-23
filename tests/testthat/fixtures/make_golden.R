# ONE-OFF generator for tests/testthat/fixtures/golden.rds.
#
# RULE: run ONCE against the unmodified code (step S1 of the review plan).
# NEVER regenerate to make a failing test pass. Regenerate only under an
# approved gate (G1..G9; G8 = exact absorption probabilities, G9 = optimiser and
# local bootstrap refit) and record the gate in the commit message.
#
# Usage (from repo root):  Rscript tests/testthat/fixtures/make_golden.R

source(here::here("tests", "testthat", "helper-setup.R"))

rows <- resampled_pars_fin[1:20, ]
scen_ly <- function(id) {
  long <- run_scenarios(rows, risk = risk_multipliers["base"])
  long$ly[long$scenario == id]
}

golden <- list(
  life_add_baseline = life_add(),
  life_add_no_screening = life_add(p_NA = 1, p_CA = 0),
  screen_by_adherence_default = unlist(screen_by_adherence()),
  # KNOWN BUG C1 (gate G1): pins the CURRENT (wrong) output
  progress_untreated_100_120 = unlist(progress_untreated(n = 100, t = 120)),
  # G7: at the CSV point estimates (pre-G7 hypo_probs: 12.587514779190872).
  # G8/G9: exact absorption probabilities and the refit CSV point
  neg_log_lik_csv = neg_log_lik(csv_par),
  scen_base = scen_ly("baseline"),
  scen_1a = scen_ly("real_adh"),
  scen_1b = scen_ly("opt_adh"),
  scen_2 = scen_ly("diag_sens"),
  scen_3 = scen_ly("treat"),
  life_tibble_main = life_tibble(csv_main),
  life_tibble_25 = life_tibble(csv_25),
  param_names = names(readRDS(paper_path("data-derived", "param_names.rds")))
)

saveRDS(golden, golden_path())
message("wrote golden.rds")
