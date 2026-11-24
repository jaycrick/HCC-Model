# Sweep figures (file R/plot_sweeps.R)
# Bootstrap sweeps and their band plots, and the untreated stage
# distributions (analysis/02, analysis/03).
# Requires R/setup.R, R/model_core.R and R/plot_theme.R.

# ----- Bootstrap sweeps -------------------------------------------------------

#' One model evaluation: `fn` at default parameters overridden by `overrides`,
#' with `arg` set to `x`. `arg` is a model parameter, or a keyword argument of
#' `fn` such as the screening interval `m`.
eval_at <- function(fn, x, arg, overrides, t) {
  swept <- setNames(list(x), arg)
  if (arg %in% names(params_def)) {
    overrides <- utils::modifyList(overrides, swept)
    swept <- list()
  }
  do.call(fn, c(list(params = modify_params(params_def, overrides), t = t), swept))
}

#' Evaluate a model function over a grid, at default parameters and for every
#' bootstrap draw.
#'
#' Draws go through `prepare_draws()`, so Sym_B = SYM_B_RATIO * Sym_A in each
#' (finding C7). The point estimate uses the default parameters and the scalar
#' path; the draws are evaluated by the vectorised engine, one call per grid
#' value, on a column-parameter list built once.
#'
#' @param fn     `life_add` or `screen_by_adherence`.
#' @param grid   Values of `arg`.
#' @param arg    Name of the swept parameter (or "m").
#' @param boot   Tibble of bootstrap parameter draws.
#' @param fixed  Named list of parameter overrides held fixed.
#' @param t      Follow-up months.
#' @return List with `point` (one result per grid value) and `draws` (per grid
#'   value, a numeric vector over draws for `life_add`, or a draws x 7 matrix
#'   for `screen_by_adherence`).
sweep_model <- function(fn, grid, arg, boot, fixed = list(), t = 120) {
  vec_fn <- if (identical(fn, life_add)) {
    life_add_vec
  } else if (identical(fn, screen_by_adherence)) {
    screen_by_adherence_mat
  } else {
    stop("`fn` must be life_add or screen_by_adherence")
  }
  pd <- prepare_draws(boot)
  # defaults, then each draw's columns, then the fixed overrides
  cols <- modify_params(params_def, as.list(pd))
  cols <- modify_params(cols, fixed)
  m <- 6 # default ultrasound interval of life_add()
  eval_grid <- function(x) {
    p <- cols
    if (arg == "m") {
      m <- x
    } else {
      p <- modify_params(p, setNames(list(x), arg))
    }
    if (nrow(pd) == 0) {
      r <- vec_fn(params_def, 100, t, m) # no draws: keep the result's shape
      return(if (is.matrix(r)) r[0, , drop = FALSE] else r[0])
    }
    vec_fn(p, 100, t, m)
  }
  list(
    point = lapply(grid, eval_at, fn = fn, arg = arg, overrides = fixed, t = t),
    draws = lapply(grid, eval_grid)
  )
}

#' Life-years gained over a grid of one parameter, with a bootstrap band.
#'
#' @param grid   Values of the swept parameter.
#' @param arg    Its name: a model parameter ("p_CA", "p_Att") or "m".
#' @param boot   Tibble of bootstrap draws.
#' @param probs  Band quantiles (default 2.5 % and 97.5 %).
#' @param fixed  Named list of parameters held fixed (e.g. p_NA = 0).
#' @param t      Follow-up months.
#' @return Tibble `x`, `point`, `lower`, `upper`; attribute `probs`.
sweep_life_years <- function(
  grid,
  arg,
  boot,
  probs = c(0.025, 0.975),
  fixed = list(),
  t = 120
) {
  res <- sweep_model(life_add, grid, arg, boot, fixed, t)
  ci <- vapply(res$draws, summarise_ci, numeric(2), probs = probs)
  out <- tibble::tibble(
    x = grid,
    point = unlist(res$point),
    lower = ci["lower", ],
    upper = ci["upper", ]
  )
  attr(out, "probs") <- probs
  out
}

#' Number diagnosed per BCLC stage over a grid, with a bootstrap band.
#'
#' @inheritParams sweep_life_years
#' @return Long tibble `x`, `stage` (BCLC 0 .. D), `number` (point estimate),
#'   `lower`, `upper`; attribute `probs`. Default band is the IQR.
sweep_stage_diagnoses <- function(
  grid,
  arg,
  boot,
  probs = c(0.25, 0.75),
  fixed = list(),
  t = 120
) {
  stages <- stage_list[2:6]
  res <- sweep_model(screen_by_adherence, grid, arg, boot, fixed, t)
  rows <- lapply(seq_along(grid), function(k) {
    mat <- res$draws[[k]]
    ci <- vapply(stages, function(s) summarise_ci(mat[, s], probs), numeric(2))
    tibble::tibble(
      x = grid[k],
      stage = stages,
      number = unname(unlist(res$point[[k]][stages])),
      lower = unname(ci["lower", ]),
      upper = unname(ci["upper", ])
    )
  })
  out <- dplyr::bind_rows(rows)
  attr(out, "probs") <- probs
  out
}

#' Line plus bootstrap band for a `sweep_life_years()` result.
#'
#' The legend names the band's real quantile range ("95%" or "IQR").
#'
#' @param sweep       Output of `sweep_life_years()`.
#' @param x_lab       x-axis label.
#' @param title,caption  Optional title and caption.
#' @param x_scale     x scale (default: continuous, flush with the axis).
#' @param y_scale     y scale (default: flush with the axis; NULL for ggplot's).
#' @param linewidth   Width of the point-estimate line.
#' @param fill_title  Legend title.
#' @return A ggplot.
plot_sweep_band <- function(
  sweep,
  x_lab,
  title = NULL,
  caption = NULL,
  x_scale = scale_x_continuous(expand = expansion(0)),
  y_scale = scale_y_continuous(expand = expansion(0)),
  linewidth = 1.25,
  fill_title = "Confidence\nBound"
) {
  ggplot(sweep, aes(.data$x, .data$point)) +
    geom_line(linewidth = linewidth) +
    geom_ribbon(
      aes(
        ymin = .data$lower,
        ymax = .data$upper,
        fill = band_label(attr(sweep, "probs"))
      ),
      alpha = 0.2
    ) +
    scale_fill_manual(values = "black") +
    x_scale +
    y_scale +
    labs(
      title = title,
      x = x_lab,
      y = "Additional Life Years Gained From\nIntervention per Cohort of 100",
      caption = caption,
      fill = fill_title
    ) +
    theme_hcc(legend_title_centre = TRUE)
}

#' Diagnoses per stage against attendance probability, optionally with the
#' bootstrap band as ribbons.
#'
#' @param sweep    Output of `sweep_stage_diagnoses()`.
#' @param title    Plot title.
#' @param caption  Caption.
#' @param ribbon   Draw the band? (Its quantile range is named in the caption.)
#' @return A ggplot.
plot_stage_diagnoses <- function(sweep, title, caption, ribbon = FALSE) {
  p <- ggplot(
    sweep,
    aes(
      .data$x,
      .data$number,
      group = .data$stage,
      colour = .data$stage,
      linetype = .data$stage,
      fill = .data$stage
    )
  ) +
    geom_line(linewidth = 1.5) +
    ggtitle(title) +
    scale_x_continuous(labels = label_percent(), expand = expansion(0))
  if (ribbon) {
    p <- p +
      my_fill +
      geom_ribbon(aes(ymin = .data$lower, ymax = .data$upper), alpha = 0.2)
    caption <- paste0(
      caption,
      "\nRibbon: ",
      band_label(attr(sweep, "probs")),
      " across bootstrap draws"
    )
  }
  p +
    my_line +
    my_scale +
    labs(
      fill = "Stage of\nDiagnosis",
      caption = caption,
      x = "Probability of Attending\nScheduled Ultrasound Screenings",
      y = "Number Diagnosed by Stage"
    ) +
    theme_hcc(11, legend_title_centre = TRUE)
}

# ----- Untreated stage distributions ------------------------------------------

#' Untreated stage distribution after `years` for a range of annual HCC risk.
#'
#' @param years  Horizon in years.
#' @param risk   Annual HCC risks to evaluate.
#' @param n      Cohort size.
#' @return Wide tibble: `x` (risk) and one column per stage.
stage_distribution_by_risk <- function(
  years,
  risk = seq(0, 0.08, length.out = 1000),
  n = 100
) {
  rows <- purrr::map(risk, ~ progress_untreated(n = n, t = years * 12, p_HCC = .x))
  dplyr::bind_cols(tibble::tibble(x = risk), dplyr::bind_rows(rows))
}

#' Untreated stage distribution month by month for `years` at default risk.
#'
#' @param years  Horizon in years.
#' @param n      Cohort size.
#' @return Wide tibble: `x` (years) and one column per stage.
stage_distribution_over_time <- function(years, n = 100) {
  months <- 0:(years * 12)
  rows <- purrr::map(months, ~ progress_untreated(n = n, t = .x))
  # months to years: Ref: none (unit conversion)
  dplyr::bind_cols(tibble::tibble(x = months / 12), dplyr::bind_rows(rows))
}

#' Stacked-panel line plot of the untreated stage distribution (Cirrhosis and
#' Mortality on top, the BCLC stages below).
#'
#' @param tbl      Output of `stage_distribution_by_risk()` or
#'                 `stage_distribution_over_time()`.
#' @param x_lab,title,caption  Labels.
#' @param x_scale  x scale.
#' @return A ggplot.
plot_stage_distribution <- function(
  tbl,
  x_lab,
  title,
  caption,
  x_scale = scale_x_continuous(expand = expansion(0))
) {
  long <- tidyr::pivot_longer(
    tbl,
    dplyr::all_of(stage_list),
    names_to = "HCC_stage",
    values_to = "percent"
  )
  long$HCC_stage <- factor(long$HCC_stage, levels = stage_list)
  long$graph_location <- factor(
    ifelse(long$HCC_stage %in% c("Cirrhosis", "Mortality"), "top", "bottom"),
    levels = c("top", "bottom")
  )
  ggplot(
    long,
    aes(
      x = .data$x,
      y = .data$percent,
      group = .data$HCC_stage,
      colour = .data$HCC_stage,
      linetype = .data$HCC_stage
    )
  ) +
    geom_line(linewidth = 1.5) +
    my_scale +
    my_line +
    x_scale +
    scale_y_continuous(
      labels = label_percent(scale = 1),
      expand = expansion(c(0, 0.01))
    ) +
    guides(colour = guide_legend("HCC Stage"), linetype = guide_legend("HCC Stage")) +
    facet_wrap(~graph_location, dir = "v", scales = "free_y") +
    labs(
      x = x_lab,
      y = "Percentage of Population per Stage",
      title = title,
      caption = caption
    ) +
    theme_hcc(11, legend_title_centre = TRUE) +
    theme(legend.position = "bottom", strip.text.x = element_blank())
}
