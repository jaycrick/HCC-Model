# ONE-OFF generator for tests/testthat/fixtures/adherence_golden.rds.
#
# Pins the free 3-class adherence fit and the Att_Low = 0 fit on the committed
# aggregate table papers/stochastic/data/adherence_counts.csv (no raw data).
# RULE: run ONCE. NEVER regenerate to make a failing test pass. Regenerate only
# under an approved gate and record the gate in the commit message.
#
# REGENERATED 2026-11-24: the input table changed intentionally (data-rule
# change per plan amendment: scans and end converted to whole months before the
# merge/push rule). The fit code did not change.
#
# REGENERATED 2026-11-25 (approved): the base table switched to the pre-COVID
# cohort (follow-up capped at 2020-03-01) and gained the exposure column T. The
# fit code did not change (the likelihood ignores T).
#
# Usage (from repo root):  Rscript tests/testthat/fixtures/make_adherence_golden.R

source(here::here("tests", "testthat", "helper-setup.R"))

tab <- read.csv(here::here("papers", "stochastic", "data", "adherence_counts.csv"))
stopifnot(identical(names(tab), c("A", "F", "T", "n")))

f3 <- fit_adherence(tab)
f3_0 <- fit_adherence(tab, fix = list(Att_Low = 0))
golden <- list(
  free = list(par = f3$par, loglik = f3$loglik),
  att_low_0 = list(par = f3_0$par, loglik = f3_0$loglik),
  provenance = list(r_version = R.version.string)
)

saveRDS(golden, here::here("tests", "testthat", "fixtures", "adherence_golden.rds"))
message("wrote adherence_golden.rds")
