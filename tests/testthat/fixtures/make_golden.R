# ONE-OFF generator for tests/testthat/fixtures/golden.rds.
#
# RULE: run ONCE against the unmodified code (step S1 of the review plan).
# NEVER regenerate to make a failing test pass. Regenerate only under an
# approved gate (G1..G7) and record the gate in the commit message.
#
# Usage (from repo root):  Rscript tests/testthat/fixtures/make_golden.R

source(here::here("tests", "testthat", "helper-setup.R"))

rows <- resampled_pars_fin[1:20, ]

golden <- list(
  life_add_baseline = life_add(),
  life_add_no_screening = life_add(p_NA = 1, p_CA = 0),
  screen_by_adherence_default = unlist(screen_by_adherence()),
  # KNOWN BUG C1 (gate G1): pins the CURRENT (wrong) output
  progress_untreated_100_120 = unlist(progress_untreated(n = 100, t = 120)),
  neg_log_lik_hypo_probs = neg_log_lik(hypo_probs, meds = meds),
  scen_base = purrr::pmap_dbl(rows, base_scen_fun),
  scen_1a = purrr::pmap_dbl(rows, scen_1a_fun),
  scen_1b = purrr::pmap_dbl(rows, scen_1b_fun),
  scen_2 = purrr::pmap_dbl(rows, scen_2_fun),
  scen_3 = purrr::pmap_dbl(rows, scen_3_fun),
  life_tibble_main = life_tibble(csv_main),
  life_tibble_25 = life_tibble(csv_25),
  param_names = names(readRDS(here::here("data-derived", "param_names.rds")))
)

saveRDS(golden, golden_path())
message("wrote golden.rds")
