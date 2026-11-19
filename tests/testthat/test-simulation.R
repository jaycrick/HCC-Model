# The stochastic simulation must agree, on average, with the deterministic
# cohort model (R/model_core.R): same operations, same order, same
# probabilities. Counts are compared per stage at n = 1000 (the deterministic
# model is run with the same n, so no rescaling is needed).

source(here::here("R", "simulation.R"))

test_that("mean of 500 simulate_by_adherence() runs matches the model", {
  n <- 1000
  n_sim <- 500
  set.seed(20240801)
  sims <- t(replicate(n_sim, simulate_by_adherence(n = n)))
  expect_equal(dim(sims), c(n_sim, 7))
  expect_true(all(rowSums(sims) <= n)) # the rest stay undiagnosed

  det <- as.numeric(screen_by_adherence(n = n))
  se <- apply(sims, 2, sd) / sqrt(n_sim)
  mn <- unname(colMeans(sims))
  for (k in seq_along(det)) {
    if (se[k] == 0) {
      expect_equal(mn[k], det[k], info = stage_list[k])
    } else {
      expect_lt(abs(mn[k] - det[k]) / se[k], 4, label = stage_list[k])
    }
  }
})

test_that("simulate_cohort() rejects unknown parameter names", {
  expect_error(simulate_cohort(0.5, p_att = 1), "unknown parameter")
})
