# Tests for the vectorised engine (screen_cohort_mat, screen_by_adherence_mat,
# life_add_vec) against loops of the scalar public functions (plan EFF-1, S2).

boot50 <- prepare_draws(head(resampled_pars_fin, 50))
# Full parameter lists (one element per parameter, length 50), merged over the
# defaults for any parameter the draws do not carry.
pars50 <- modifyList(params_def, as.list(boot50))
row_pars <- function(i) modifyList(params_def, as.list(boot50[i, ]))

for (m in c(1, 6, 24)) {
  test_that(paste0("screen_cohort_mat matches scalar loop, m = ", m), {
    for (att in c(0, 0.37, 1)) {
      mat <- screen_cohort_mat(att, pars50, n = 100, t = 120, m = m)
      ref <- lapply(1:50, function(i) {
        screen_cohort(att, n = 100, t = 120, m = m, params = row_pars(i))
      })
      dx <- do.call(rbind, lapply(ref, `[[`, "dx"))
      und <- do.call(rbind, lapply(ref, `[[`, "undiagnosed"))
      # Both sides use the same element-wise code path, so identical holds.
      expect_identical(unname(mat$dx), unname(dx))
      expect_identical(unname(mat$undiagnosed), unname(und))
    }
  })

  test_that(paste0("screen_by_adherence_mat matches scalar loop, m = ", m), {
    mat <- screen_by_adherence_mat(pars50, n = 100, t = 120, m = m)
    ref <- do.call(
      rbind,
      lapply(1:50, function(i) {
        unlist(screen_by_adherence(m = m, params = row_pars(i)))
      })
    )
    expect_identical(unname(mat), unname(ref))
  })

  test_that(paste0("life_add_vec matches scalar loop, m = ", m), {
    vec <- life_add_vec(pars50, n = 100, t = 120, m = m)
    ref <- vapply(
      1:50,
      function(i) life_add(m = m, params = row_pars(i)),
      numeric(1)
    )
    expect_identical(vec, ref)
  })
}

test_that("vectorised engine matches the previous matrix-product model", {
  # Reference: the original `vec %*% tmat` loop. An element-wise update is not
  # bit-identical to BLAS (FMA on arm64; measured ~1e-15 relative), hence
  # tolerance 1e-14 rather than expect_identical.
  ref_loop <- function(Pr_Att, params, n = 100, t = 120, m = 6) {
    tmat <- transition_matrix(
      annual_to_monthly(params$p_HCC),
      unlist(params[STAGE_P_NAMES])
    )
    sp <- c(0, unlist(params[STAGE_SYM_NAMES]), 1)
    us_p <- c(0, unlist(params[STAGE_SENS_NAMES]), 1)
    att <- c(0, rep(Pr_Att, 5), 1)
    vec <- c(n, rep(0, 6))
    Dx <- rep(0, 7)
    for (i in 1:t) {
      vec <- vec %*% tmat
      sym <- vec * sp
      Dx <- Dx + sym
      vec <- vec - sym
      if (i %% m == 0) {
        us <- (vec * att) * us_p
        Dx <- Dx + us
        vec <- vec - us
      }
    }
    list(dx = as.numeric(Dx), undiagnosed = as.numeric(vec))
  }
  for (m in c(1, 6, 24)) {
    for (att in c(0, 0.37, 1)) {
      new <- screen_cohort(att, m = m)
      old <- ref_loop(att, params_def, m = m)
      expect_equal(unname(new$dx), old$dx, tolerance = 1e-14)
      expect_equal(unname(new$undiagnosed), old$undiagnosed, tolerance = 1e-14)
    }
  }
})

test_that("length-1 params give N = 1 and mass is conserved", {
  for (m in c(1, 6, 24)) {
    for (att in c(0, 0.37, 1)) {
      res <- screen_cohort_mat(att, params_def, n = 100, t = 120, m = m)
      expect_identical(nrow(res$dx), 1L)
      expect_identical(ncol(res$dx), 7L)
      expect_identical(colnames(res$dx), stage_list)
      expect_lt(abs(rowSums(res$dx) + rowSums(res$undiagnosed) - 100), 1e-12)
    }
  }
  expect_identical(nrow(screen_by_adherence_mat(params_def, 100, 120, 6)), 1L)
  expect_length(life_add_vec(params_def, 100, 120, 6), 1L)
  res <- screen_cohort_mat(0.37, pars50, n = 100, t = 120, m = 6)
  expect_equal(
    rowSums(res$dx) + rowSums(res$undiagnosed),
    rep(100, 50),
    tolerance = 1e-12
  )
})

test_that("sweep_life_years on 5 draws equals a loop of scalar life_add", {
  b5 <- resampled_pars_fin[1:5, ]
  scalar_ly <- function(i, ...) {
    p <- as.list(draws(b5)[i, ])
    do.call(life_add, c(p, list(Sym_B = 8 * b5$Sym_A[i]), list(...)))
  }
  q <- c(0.025, 0.975)
  # Same element-wise engine on both sides, so tolerance 0 (identical values).
  sw <- sweep_life_years(c(1, 6, 24), "m", b5, t = 60)
  for (k in 1:3) {
    ly <- vapply(1:5, scalar_ly, 1, m = c(1, 6, 24)[k], t = 60)
    expect_equal(
      c(sw$lower[k], sw$upper[k]),
      unname(quantile(ly, q)),
      tolerance = 0
    )
  }
  sw <- sweep_life_years(
    c(0, 0.37, 1),
    "p_Att",
    b5,
    fixed = list(p_NA = 0, p_CA = 0),
    t = 60
  )
  for (k in 1:3) {
    ly <- vapply(
      1:5,
      scalar_ly,
      1,
      p_Att = c(0, 0.37, 1)[k],
      p_NA = 0,
      p_CA = 0,
      t = 60
    )
    expect_equal(
      c(sw$lower[k], sw$upper[k]),
      unname(quantile(ly, q)),
      tolerance = 0
    )
  }
})

test_that("sweep_model rejects an unsupported fn", {
  expect_error(
    sweep_model(sum, c(1, 6), "m", resampled_pars_fin[1:2, ]),
    "must be life_add or screen_by_adherence"
  )
})

test_that("zero-row boot gives NA bands and an unchanged point estimate", {
  b0 <- resampled_pars_fin[0, ]
  sw <- sweep_life_years(c(1, 6), "m", b0, t = 60)
  expect_identical(sw$point, vapply(c(1, 6), \(m) life_add(m = m, t = 60), 1))
  expect_true(all(is.na(sw$lower)) && all(is.na(sw$upper)))
  st <- sweep_stage_diagnoses(c(0, 1), "p_Att", b0, t = 60)
  expect_true(all(is.na(st$lower)) && all(is.na(st$upper)))
  expect_false(anyNA(st$number))
})

test_that("scenario modifiers on vector params equal elementwise scalar use", {
  pars <- modifyList(
    params_def,
    as.list(prepare_draws(head(resampled_pars_fin, 8)))
  )
  for (nm in names(scenario_defs)) {
    v <- scenario_defs[[nm]](pars)
    for (i in 1:8) {
      sc <- scenario_defs[[nm]](lapply(
        pars,
        function(x) x[[min(i, length(x))]]
      ))
      got <- lapply(v, function(x) rep_len(x, 8)[[i]])
      expect_identical(got, sc, label = paste(nm, "draw", i))
    }
  }
})
