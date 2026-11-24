# Tests for R/adherence.R (likelihood, EM / optim fit, profile CI, simulator)

adh_tab <- function() {
  read.csv(here::here("papers", "stochastic", "data", "adherence_counts.csv"))
}
adh_par <- c(
  Pi_Low = 0.3,
  Pi_High = 0.3,
  Att_Low = 0.02,
  Att_Med = 0.25,
  Att_High = 0.8
)

# Brute force: month-by-month path likelihood from explicit scan months, not
# using due_exposure() for the probability itself.
path_lik <- function(par, scans, horizon, m = 6) {
  u <- adh_unpack(par)
  sum(vapply(
    1:3,
    function(k) {
      d <- 0
      lik <- 1
      for (t in seq_len(horizon) - 1) {
        if (t < d) {
          next
        }
        if (t %in% scans) {
          lik <- lik * u$q[k]
          d <- t + m
        } else {
          lik <- lik * (1 - u$q[k])
        }
      }
      u$pi[k] * lik
    },
    numeric(1)
  ))
}

test_that("loglik equals brute-force path likelihood", {
  people <- list(
    list(scans = integer(0), hz = 12),
    list(scans = c(0), hz = 5),
    list(scans = c(3, 10), hz = 20),
    list(scans = c(0, 6, 12), hz = 30),
    list(scans = c(2), hz = 7),
    list(scans = integer(0), hz = 12),
    list(scans = c(4, 11), hz = 18)
  )
  ex <- t(vapply(people, function(p) due_exposure(p$scans, p$hz), numeric(2)))
  tab <- adherence_counts(ex[, "A"], ex[, "F"])
  brute <- sum(vapply(
    people,
    function(p) log(path_lik(adh_par, p$scans, p$hz)),
    numeric(1)
  ))
  expect_equal(adherence_loglik(adh_par, tab), brute, tolerance = 1e-12)
  expect_equal(adherence_loglik(as.list(adh_par), tab), brute, tolerance = 1e-12)
})

test_that("posterior rows sum to one", {
  r <- adherence_posterior(adh_par, adh_tab())
  expect_equal(unname(rowSums(r)), rep(1, nrow(r)), tolerance = 1e-12)
})

test_that("invalid pi gives -Inf", {
  bad <- replace(adh_par, "Pi_Low", 0.8)
  expect_equal(adherence_loglik(bad, adh_tab()), -Inf)
})

test_that("EM and optim agree on the committed table", {
  tab <- adh_tab()
  e <- fit_adherence(tab, "em", n_starts = 6)
  o <- fit_adherence(tab, "optim", n_starts = 6)
  expect_equal(e$loglik, o$loglik, tolerance = 1e-6)
  expect_equal(e$par, o$par, tolerance = 1e-3)
  expect_equal(e$n_par_free, 5)
  expect_identical(names(e$par), ADH_PAR_NAMES)
})

test_that("label ordering holds", {
  e <- fit_adherence(adh_tab(), "em", n_starts = 6)
  expect_true(!is.unsorted(unname(e$par[c("Att_Low", "Att_Med", "Att_High")])))
  p <- c(Pi_Low = 0.2, Pi_High = 0.5, Att_Low = 0.7, Att_Med = 0.4, Att_High = 0.1)
  s <- adh_order_classes(p)
  expect_equal(unname(s), c(0.5, 0.2, 0.1, 0.4, 0.7))
  s2 <- adh_order_classes(
    c(Pi_Low = 0.2, Pi_High = 0.5, Att_Low = 0, Att_Med = 0.7, Att_High = 0.1),
    "Att_Low"
  )
  expect_equal(s2[["Att_Low"]], 0)
  expect_equal(s2[["Att_Med"]], 0.1)
})

test_that("recovery on a simulated table, EM agrees with optim", {
  set.seed(20261124)
  ex <- sample(20:160, 10000, replace = TRUE)
  sim <- simulate_adherence_people(adh_par, ex)
  tab <- adherence_counts(sim$A, sim$F)
  e <- fit_adherence(tab, "em", n_starts = 4)
  o <- fit_adherence(tab, "optim", n_starts = 4)
  expect_equal(e$loglik, o$loglik, tolerance = 1e-6)
  expect_lt(max(abs(e$par[names(adh_par)] - adh_par)), 0.03)
})

test_that("simulator reproduces due_exposure of its own scan months", {
  set.seed(1)
  # q = 1: attend as soon as due, so scans fall at 0, m, 2m, ...
  s <- simulate_adherence_people(
    c(Pi_Low = 0.1, Pi_High = 0.1, Att_Low = 1, Att_Med = 1, Att_High = 1),
    c(1, 6, 7, 13, 30)
  )
  expected <- t(vapply(
    c(1, 6, 7, 13, 30),
    function(hz) {
      due_exposure(seq(0, hz - 1, by = 6), hz)
    },
    numeric(2)
  ))
  expect_equal(s$A, unname(expected[, "A"]))
  expect_equal(s$F, unname(expected[, "F"]))
  # q = 0: never attends, F = exposure
  z <- simulate_adherence_people(
    c(Pi_Low = 0.1, Pi_High = 0.1, Att_Low = 0, Att_Med = 0, Att_High = 0),
    c(5, 50)
  )
  expect_equal(z$A, c(0, 0))
  expect_equal(z$F, c(5, 50))
})

test_that("fix = list(Att_Low = 0) works", {
  tab <- adh_tab()
  free <- fit_adherence(tab, "em", n_starts = 6)
  f <- fit_adherence(tab, "em", fix = list(Att_Low = 0), n_starts = 6)
  o <- fit_adherence(tab, "optim", fix = list(Att_Low = 0), n_starts = 6)
  expect_equal(f$par[["Att_Low"]], 0)
  expect_equal(f$n_par_free, 4)
  expect_lte(f$loglik, free$loglik + 1e-9)
  expect_equal(f$loglik, o$loglik, tolerance = 1e-6)
  expect_error(fit_adherence(tab, "em", fix = list(Pi_Low = 0.3)), "optim")
  fo <- fit_adherence(tab, "optim", fix = list(Pi_Low = 0.3), n_starts = 3)
  expect_equal(fo$par[["Pi_Low"]], 0.3)
  expect_lte(fo$loglik, free$loglik + 1e-9)
})

test_that("boundary q = 0 / 1 give finite loglik", {
  tab <- data.frame(A = c(0, 0, 2), F = c(10, 0, 0), n = c(3, 1, 2))
  p <- c(Pi_Low = 0.2, Pi_High = 0.3, Att_Low = 0, Att_Med = 0.2, Att_High = 1)
  expect_true(is.finite(adherence_loglik(p, tab)))
  # A = F = 0 person (zero exposure) has likelihood 1 whatever q
  tab0 <- data.frame(A = 0, F = 0, n = 4)
  expect_equal(adherence_loglik(p, tab0), 0)
  # impossible under the parameters: attended with q = 0 for every class
  q0 <- c(Pi_Low = 0.2, Pi_High = 0.3, Att_Low = 0, Att_Med = 0, Att_High = 0)
  expect_equal(adherence_loglik(q0, data.frame(A = 1, F = 0, n = 1)), -Inf)
})

test_that("profile CI brackets the estimate and respects boundaries", {
  tab <- adh_tab()
  e <- fit_adherence(tab, "em", n_starts = 4)
  ci <- adherence_profile_ci(e, tab, "Att_Med")
  expect_lt(ci$lower, ci$estimate)
  expect_gt(ci$upper, ci$estimate)
  prof <- adh_fit_ll(tab, c(Att_Med = ci$upper), list(e$par))
  expect_equal(prof, ci$cutoff, tolerance = 1e-4)
  # tiny table: the interval can reach the [0, 1] bounds
  small <- data.frame(A = c(0, 3), F = c(12, 2), n = c(5, 5))
  es <- fit_adherence(small, "em", n_starts = 4)
  cs <- adherence_profile_ci(es, small, "Pi_Low")
  expect_gte(cs$lower, 0)
  expect_lte(cs$upper, 1)
})

test_that("bootstrap helper draws counts that sum to N and fits run", {
  tab <- adh_tab()
  set.seed(3)
  cnt <- draw_bootstrap_counts(tab, 3)
  expect_equal(unname(rowSums(cnt)), rep(sum(tab$n), 3))
  b <- adherence_bootstrap(tab, 2, list(n_starts = 3), cnt)
  expect_equal(dim(b), c(2, 7))
  expect_true(all(is.finite(b)))
})
