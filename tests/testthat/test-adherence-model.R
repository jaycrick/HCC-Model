# "due-until-attend" adherence engine (screen_cohort_due_mat, plan 2026-11-24
# S3): consistency with the scheduled engine, conservation, person-level scan
# counts and the stochastic analogue.

source(here::here("R", "simulation.R"))

due_params <- function(...) {
  modify_params(
    c(
      params_def,
      list(
        Pi_Low = 0.3,
        Pi_Med = 0.45,
        Pi_High = 0.25,
        Att_Low = 0.02,
        Att_Med = 0.2,
        Att_High = 0.6
      )
    ),
    list(...)
  )
}

test_that("adherence = 'scheduled' is identical to the default", {
  expect_identical(screen_by_adherence(), screen_by_adherence(adherence = "scheduled"))
  expect_identical(
    screen_by_adherence(adherence = "scheduled", init = "gap"),
    screen_by_adherence()
  )
  expect_identical(life_add(), life_add(adherence = "scheduled"))
  expect_identical(
    life_add_vec(params_def, 100, 120, 6),
    life_add_vec(params_def, 100, 120, 6, adherence = "scheduled")
  )
})

test_that("due engine, init = 'gap', matches the scheduled engine at q = 0, 1", {
  for (m in c(3, 6, 12)) {
    for (t in c(m, m + 1, 2 * m - 1, 3 * m, 40, 61)) {
      for (qq in c(0, 1)) {
        a <- screen_cohort_due_mat(qq, params_def, 1000, t, m, "gap")
        b <- screen_cohort_mat(qq, params_def, 1000, t, m)
        expect_equal(a$dx, b$dx, tolerance = 1e-12, info = paste(m, t, qq))
        expect_equal(a$undiagnosed, b$undiagnosed, tolerance = 1e-12)
      }
    }
  }
})

test_that("people are conserved", {
  for (init in c("due", "gap")) {
    r <- screen_cohort_due_mat(0.3, params_def, 500, 77, 6, init)
    expect_equal(sum(r$dx) + sum(r$undiagnosed), 500, tolerance = 1e-9)
  }
  dx <- screen_by_adherence_mat(due_params(), 500, 77, 6, "due")
  r <- screen_cohort_due_mat(c(0.02, 0.2, 0.6), params_def, 500, 77, 6)
  expect_equal(sum(r$dx) + sum(r$undiagnosed), 3 * 500, tolerance = 1e-9)
  expect_lte(sum(dx), 500)
})

test_that("expected scans match the person-level model (exposure = t)", {
  # No progression: everyone stays in Cirrhosis, so only the clock acts.
  p0 <- modify_params(params_def, list(p_HCC = 0))
  t <- 60
  m <- 6
  q <- c(0.05, 0.3, 0.8)
  sc <- screen_cohort_due_mat(q, p0, 1, t, m, "due")$scans
  set.seed(1)
  par <- c(Pi_Low = 0, Pi_High = 0, Att_Low = q[1], Att_Med = q[2], Att_High = q[3])
  for (k in 1:3) {
    pk <- par
    pk[c("Pi_Low", "Pi_High")] <- if (k == 1) {
      c(1, 0)
    } else if (k == 3) {
      c(0, 1)
    } else {
      c(0, 0)
    }
    sim <- simulate_adherence_people(pk, rep(t, 20000), m)
    se <- sd(sim$A) / sqrt(20000)
    expect_lt(abs(mean(sim$A) - sc[k]), 4 * se, label = paste("class", k))
  }
})

test_that("stochastic due engine matches the deterministic one", {
  n <- 1000
  n_sim <- 200
  prm <- due_params()
  set.seed(20241124)
  sims <- t(replicate(
    n_sim,
    simulate_by_adherence(n = n, t = 60, params = prm, adherence = "due")
  ))
  det <- as.numeric(screen_by_adherence(n = n, t = 60, params = prm, adherence = "due"))
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

test_that("'due' requires Pi_* and Att_*; 'scheduled' does not", {
  expect_error(screen_by_adherence(adherence = "due"), "Att_Low, Att_Med, Att_High")
  expect_error(life_add(adherence = "due"), "requires parameter")
  expect_error(simulate_by_adherence(adherence = "due"), "requires parameter")
  expect_false("Att_Low" %in% names(params_def))
  expect_no_error(screen_by_adherence())
})

test_that("'due' checks the class shares sum to 1 and ignores p_NA, p_CA, p_Att", {
  bad <- due_params(Pi_Med = 0.5)
  expect_error(screen_by_adherence_mat(bad, 100, 60, 6, "due"), "must equal 1")
  expect_error(simulate_by_adherence(params = bad, adherence = "due"), "must equal 1")
  a <- screen_by_adherence_mat(due_params(), 100, 60, 6, "due")
  b <- screen_by_adherence_mat(
    due_params(p_NA = 0.9, p_CA = 0.05, p_Att = 0.99),
    100,
    60,
    6,
    "due"
  )
  expect_identical(a, b)
})
