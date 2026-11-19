# Tests for R/plot_*.R and the label tables in R/setup.R.
# Gates G2-G5 (docs in CHANGES.md) are each covered by a labelled test here.

# ----- G2 (finding C5): TeX_labels subscripts match the parameter names -------

# The plotmath string a label for parameter `nm` must deparse to: base symbol
# (italic for p_* and L_*), subscript = everything after the first underscore
# (`NA` needs backticks in plotmath).
expected_label <- function(nm) {
  parts <- strsplit(nm, "_", fixed = TRUE)[[1]]
  base <- if (parts[1] %in% c("p", "L")) {
    paste0("italic(", parts[1], ")")
  } else {
    parts[1]
  }
  sub <- if (parts[2] == "NA") "`NA`" else parts[2]
  paste0(base, "[", sub, "]")
}

test_that("G2: TeX_labels has one entry per parameter, in parameter order", {
  expect_identical(names(TeX_labels), names(params_def))
})

test_that("G2: every TeX_labels entry renders its own parameter name", {
  for (nm in names(params_def)) {
    expect_identical(
      paste(deparse(TeX_labels[[nm]]), collapse = ""),
      expected_label(nm),
      label = paste("rendered label of", nm)
    )
  }
})

test_that("G2: p_B renders as p_B, not p_C (finding C5)", {
  expect_identical(deparse(TeX_labels[["p_B"]]), "italic(p)[B]")
  expect_false(identical(TeX_labels[["p_B"]], TeX_labels[["p_C"]]))
})

# ----- G5(b) (finding C9): one parameter -> group lookup ----------------------

test_that("G5b: every parameter has exactly one group, label and row", {
  expect_setequal(param_lookup$parameter, names(params_def))
  expect_false(anyDuplicated(param_lookup$parameter) > 0)
  expect_true(all(!is.na(param_lookup$group)))
  expect_true(all(nzchar(param_lookup$label)))
  expect_true(all(param_lookup$group %in% names(group_labels)))
})

test_that("G5b: adherence parameters are in the Adherence group", {
  grp <- function(nm) param_lookup$group[param_lookup$parameter == nm]
  for (nm in ADHERENCE_PARAMS) {
    expect_identical(grp(nm), "Adherence", label = paste("group of", nm))
  }
  expect_identical(grp("p_HCC"), "p")
  expect_identical(grp("p_B"), "p")
  expect_identical(grp("Sym_A"), "Sym")
  expect_identical(grp("Sens_0"), "Sens")
  expect_identical(grp("L_C"), "L")
})

test_that("G5b: facet and variable labels are derived from the lookup", {
  expect_identical(
    variable_labels,
    setNames(param_lookup$label, param_lookup$parameter)
  )
  expect_identical(
    facet_labels[names(group_labels)],
    group_labels
  )
  # the Adherence facet exists and the old, never-occurring "Pr" one is gone
  expect_true("Adherence" %in% names(facet_labels))
  expect_false("Pr" %in% names(facet_labels))
})

test_that("G5b: the non-modifiable parameters are all real parameters", {
  expect_true(all(NON_MODIFIABLE_PARAMS %in% names(params_def)))
  expect_false(any(ADHERENCE_PARAMS %in% NON_MODIFIABLE_PARAMS))
})

# ----- G3 (finding C6): the upper-p_HCC tornado plots years_high_comp ---------

param_tbl <- life_tibble(csv_main)

test_that("G3: the 'high' bound tornado plots and sorts by years_high_comp", {
  p <- plot_bound_tornado(param_tbl, "high")
  d <- p$data
  # bars: x of the plotted layer is exactly years_high_comp
  expect_equal(layer_data(p)$x, d$years_high_comp)
  expect_false(isTRUE(all.equal(d$years_high_comp, d$years_mid_comp)))
  # order: levels sorted by the summed |years_high_comp|, as for the low bound
  sorted <- levels(order_parameters(param_tbl, "years_high_comp")$Parameter)
  shown <- levels(droplevels(d$Parameter))
  expect_identical(shown, sorted[sorted %in% shown])
})

test_that("G3: the 'low' bound tornado plots and sorts by years_low_comp", {
  p <- plot_bound_tornado(param_tbl, "low")
  expect_equal(layer_data(p)$x, p$data$years_low_comp)
  sorted <- levels(order_parameters(param_tbl, "years_low_comp")$Parameter)
  shown <- levels(droplevels(p$data$Parameter))
  expect_identical(shown, sorted[sorted %in% shown])
})

test_that("G3: p_HCC and parameters without effect are not plotted", {
  p <- plot_bound_tornado(param_tbl, "high")
  expect_false("p_HCC" %in% as.character(p$data$Parameter))
  expect_true(all(p$data$years_high_comp != 0))
})

# ----- G4 (finding C7): Sym_B = SYM_B_RATIO * Sym_A in every sweep ------------

boot5 <- resampled_pars_fin[1:5, ]

test_that("G4: prepare_draws sets Sym_B = SYM_B_RATIO * Sym_A per draw", {
  pd <- prepare_draws(boot5)
  expect_identical(pd$Sym_B, SYM_B_RATIO * boot5$Sym_A)
  expect_identical(SYM_B_RATIO, 8)
  # only parameter columns remain, and the other parameters are untouched
  expect_true(all(names(pd) %in% names(params_def)))
  expect_identical(pd$p_0, boot5$p_0)
  expect_false("scen_1a" %in% names(pd))
})

test_that("G4: prepare_draws overwrites Sym_B and handles a missing Sym_A", {
  df <- tibble::tibble(Sym_A = c(0.1, 0.2), Sym_B = c(9, 9))
  expect_identical(prepare_draws(df)$Sym_B, c(0.8, 1.6))
  expect_identical(
    prepare_draws(tibble::tibble(.rows = 1))$Sym_B,
    SYM_B_RATIO * params_def$Sym_A
  )
})

# life_add on draw i with an explicit Sym_B, the sweep's expected value
explicit_ly <- function(i, ...) {
  p <- as.list(draws(boot5)[i, ])
  do.call(life_add, c(p, list(Sym_B = 8 * boot5$Sym_A[i]), list(...)))
}

test_that("G4: sweep over m on 5 draws equals life_add with explicit Sym_B", {
  grid <- c(1, 6, 24)
  sw <- sweep_life_years(grid, "m", boot5, t = 60)
  for (k in seq_along(grid)) {
    ly <- vapply(1:5, explicit_ly, 1, m = grid[k], t = 60)
    expect_equal(
      c(sw$lower[k], sw$upper[k]),
      unname(quantile(ly, c(0.025, 0.975))),
      tolerance = 1e-12,
      label = paste("band at m =", grid[k])
    )
  }
})

test_that("G4: sweeps over p_CA and p_Att (fixed p_NA, p_CA) use Sym_B too", {
  sw <- sweep_life_years(c(0.2, 0.6), "p_CA", boot5, t = 60)
  ly <- vapply(1:5, explicit_ly, 1, p_CA = 0.6, t = 60)
  expect_equal(sw$upper[2], unname(quantile(ly, 0.975)), tolerance = 1e-12)
  sw <- sweep_life_years(
    c(0, 0.5),
    "p_Att",
    boot5,
    fixed = list(p_NA = 0, p_CA = 0),
    t = 60
  )
  ly <- vapply(1:5, explicit_ly, 1, p_Att = 0.5, p_NA = 0, p_CA = 0, t = 60)
  expect_equal(sw$lower[2], unname(quantile(ly, 0.025)), tolerance = 1e-12)
})

test_that("G4: the point-estimate line is unchanged (default parameters)", {
  sw <- sweep_life_years(c(1, 6, 12), "m", boot5, t = 120)
  expect_identical(sw$point, vapply(c(1, 6, 12), \(m) life_add(m = m), 1))
  expect_identical(sw$point[2], life_add())
})

test_that("G4: stage-diagnosis sweep uses Sym_B and gives IQR", {
  sw <- sweep_stage_diagnoses(
    c(0, 1),
    "p_Att",
    boot5,
    fixed = list(p_NA = 0, p_CA = 0),
    t = 60
  )
  expect_identical(attr(sw, "probs"), c(0.25, 0.75))
  ref <- vapply(
    1:5,
    function(i) {
      p <- as.list(draws(boot5)[i, ])
      unlist(
        do.call(
          screen_by_adherence,
          c(
            p,
            list(
              Sym_B = 8 * boot5$Sym_A[i],
              p_Att = 1,
              p_NA = 0,
              p_CA = 0,
              t = 60
            )
          )
        )[["BCLC A"]]
      )
    },
    1
  )
  row <- sw[sw$x == 1 & sw$stage == "BCLC A", ]
  expect_equal(
    c(row$lower, row$upper),
    unname(quantile(ref, c(0.25, 0.75))),
    tolerance = 1e-12
  )
  expect_identical(
    row$number,
    screen_by_adherence(p_Att = 1, p_NA = 0, p_CA = 0, t = 60)[["BCLC A"]]
  )
})

test_that("G4: run_scenarios shares prepare_draws (Sym_B derived per draw)", {
  one <- run_scenarios(boot5[1, ], risk = risk_multipliers["base"])
  expect_equal(
    one$ly[one$scenario == "baseline"],
    explicit_ly(1),
    tolerance = 1e-12
  )
})

# ----- G5(a) (finding C8): honest band labels ---------------------------------

test_that("G5a: band_label names the real quantile range", {
  expect_identical(band_label(c(0.025, 0.975)), "95%")
  expect_identical(band_label(c(0.25, 0.75)), "IQR")
  expect_identical(band_label(c(0.05, 0.95)), "90%")
})

test_that("G5a: the sweep ribbon legend says 95%, never 50%", {
  sw <- sweep_life_years(c(1, 6), "m", boot5, t = 24)
  p <- plot_sweep_band(sw, "x")
  fill <- rlang::eval_tidy(p$layers[[2]]$mapping$fill)
  expect_identical(fill, "95%")
})

test_that("G5a: the stage ribbon plot names the IQR", {
  sw <- sweep_stage_diagnoses(c(0, 1), "p_Att", boot5, t = 24)
  p <- plot_stage_diagnoses(sw, "title", "caption", ribbon = TRUE)
  expect_match(p$labels$caption, "IQR")
  expect_no_match(
    plot_stage_diagnoses(sw, "t", "caption")$labels$caption,
    "IQR"
  )
})

# ----- save_figure, themes and the remaining helpers --------------------------

test_that("save_figure writes exactly the requested formats at the given dpi", {
  dir <- withr::local_tempdir()
  p <- ggplot(data.frame(x = 1:2, y = 1:2), aes(x, y)) + geom_line()
  paths <- save_figure(
    p,
    "fig",
    2,
    1,
    formats = c("jpg", "png"),
    dpi = c(jpg = 100, png = 200),
    dir = dir
  )
  expect_identical(sort(list.files(dir)), c("fig.jpg", "fig.png"))
  expect_identical(basename(paths), c("fig.jpg", "fig.png"))
  # PNG header: width and height are big-endian integers at bytes 17-24
  hdr <- readBin(file.path(dir, "fig.png"), "raw", 24)
  dims <- vapply(
    list(17:20, 21:24),
    function(b) sum(as.integer(hdr[b]) * 256^(3:0)),
    1
  )
  expect_identical(dims, c(400, 200)) # 2 x 1 inches at 200 dpi
})

test_that("save_figure defaults to jpg+png+pdf at 600 dpi and checks inputs", {
  dir <- withr::local_tempdir()
  p <- ggplot(data.frame(x = 1:2, y = 1:2), aes(x, y)) + geom_line()
  save_figure(p, "a", 1, 1, dir = dir)
  expect_setequal(list.files(dir), c("a.jpg", "a.png", "a.pdf"))
  expect_error(save_figure(p, "b", 1, 1, formats = "gif", dir = dir))
  expect_error(
    save_figure(p, "b", 1, 1, formats = "png", dpi = c(jpg = 300), dir = dir),
    "no value for format"
  )
})

test_that("dfs() builds the one-row tibble pivot_wider used to build", {
  expect_identical(dfs("p_A", 0.1), tibble::tibble(p_A = 0.1))
  expect_identical(
    dfs("p_A", 0.1, Bound = "LOWER"),
    tibble::tibble(p_A = 0.1)
  )
})

test_that("summarise_ci returns the 2.5% and 97.5% quantiles", {
  x <- c(1:100, 1000)
  expect_identical(
    summarise_ci(x),
    c(
      lower = quantile(x, 0.025, names = FALSE),
      upper = quantile(x, 0.975, names = FALSE)
    )
  )
})

test_that("stage distributions are progress_untreated() rows", {
  by_risk <- stage_distribution_by_risk(2, risk = c(0.01, 0.05))
  expect_identical(by_risk$x, c(0.01, 0.05))
  expect_equal(
    unlist(by_risk[2, stage_list]),
    unlist(progress_untreated(t = 24, p_HCC = 0.05)),
    tolerance = 1e-12
  )
  over <- stage_distribution_over_time(1)
  expect_identical(nrow(over), 13L)
  expect_identical(over$x[13], 1)
  expect_equal(sum(over[13, stage_list]), 100, tolerance = 1e-12)
  p <- plot_stage_distribution(over, "x", "t", "c")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("response curves are life_add at shifted parameters", {
  cv <- response_curves(list(
    list(name = "NA", param = "p_NA", max = 2, sign = -1)
  ))
  expect_identical(cv$x, c(0, 0.01, 0.02))
  expect_identical(cv$value[3], life_add(p_NA = params_def$p_NA - 0.02))
})

test_that("tornado, side-label and scenario plots build", {
  p <- plot_tornado(
    add_group(order_parameters(param_tbl, "years_mid_comp")),
    "years_mid_comp",
    variable_labels,
    facet = "group"
  )
  expect_no_error(ggplot2::ggplot_build(p))
  expect_identical(
    levels(add_group(param_tbl)$Type),
    c("Adherence", "Sens", "L", "p", "Sym")
  )
  expect_identical(
    nrow(modifiable_sensitivity(param_tbl)[
      modifiable_sensitivity(param_tbl)$Parameter %in% NON_MODIFIABLE_PARAMS,
    ]),
    0L
  )
})
