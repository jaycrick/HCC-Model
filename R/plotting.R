# Figure helpers (file R/plotting.R)
# Figure helpers: theme, saving, bootstrap sweeps and their band plots, stage
# distributions, tornado (sensitivity) plots and scenario figures. Replaces the
# old R/progression.R: its functions built a data table, a plot and saved it
# in one go, so each is now split into a data builder, a plot function and the
# caller's `save_figure()`. Nothing in here reads files or sets seeds; the
# analysis scripts pass in the data (e.g. bootstrap draws) explicitly.
# Requires R/setup.R, R/model_core.R and R/scenarios.R.

# ----- Theme and saving -------------------------------------------------------

#' The project's ggplot theme: black-and-white, serif, centred title.
#'
#' @param base_size           Base font size (12 for most figures; the stage
#'                            distribution and sweep figures use 11).
#' @param legend_title_centre Also centre the legend title?
#' @return A ggplot2 theme.
theme_hcc <- function(base_size = 12, legend_title_centre = FALSE) {
  th <- theme_bw(base_family = "serif", base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5))
  if (legend_title_centre) {
    th <- th + theme(legend.title = element_text(hjust = 0.5))
  }
  th
}

#' Save one figure in one or more formats.
#'
#' Every figure keeps exactly the files and resolution it always had, so the
#' formats and dpi are chosen per call: jpg only at 300 dpi, or all three.
#'
#' @param plot     A ggplot (default: the last plot drawn).
#' @param name     File name without extension, e.g. "scenarios_new".
#' @param width,height  Size in inches.
#' @param formats  Any of "jpg", "png", "pdf".
#' @param dpi      One number for all formats, or a vector named by format,
#'                 e.g. `c(jpg = 300, png = 600, pdf = 600)`.
#' @param dir      Output directory (default figures/ in the project).
#' @return Invisibly, the paths written.
save_figure <- function(
  plot = ggplot2::last_plot(),
  name,
  width,
  height,
  formats = c("jpg", "png", "pdf"),
  dpi = 600,
  dir = here::here("figures")
) {
  formats <- match.arg(formats, c("jpg", "png", "pdf"), several.ok = TRUE)
  dpis <- if (is.null(names(dpi))) rep(dpi, length(formats)) else dpi[formats]
  if (anyNA(dpis)) {
    stop("dpi has no value for format(s): ", paste(formats[is.na(dpis)]))
  }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  paths <- file.path(dir, paste0(name, ".", formats))
  for (i in seq_along(formats)) {
    ggplot2::ggsave(
      paths[i],
      plot,
      width = width,
      height = height,
      dpi = dpis[[i]]
    )
  }
  invisible(paths)
}

# ----- Intervals --------------------------------------------------------------

#' Lower and upper limits of a sample.
#'
#' @param x      Numeric sample (e.g. bootstrap results).
#' @param probs  Two quantile probabilities; 2.5 % and 97.5 % by default.
#' @return Named numeric vector `c(lower, upper)`.
summarise_ci <- function(x, probs = c(0.025, 0.975)) {
  q <- quantile(x, probs = probs, names = FALSE)
  c(lower = q[[1]], upper = q[[2]])
}

#' Legend text for a band between two quantiles: "95%" for 2.5-97.5 %,
#' "IQR" for 25-75 %.
#'
#' @param probs  Two quantile probabilities.
#' @return A string.
band_label <- function(probs) {
  if (isTRUE(all.equal(probs, c(0.25, 0.75)))) {
    return("IQR")
  }
  paste0(round(100 * (probs[2] - probs[1])), "%")
}

# ----- Bootstrap sweeps -------------------------------------------------------

# One model evaluation: `fn` at default parameters overridden by `overrides`,
# with `arg` set to `x`. `arg` is a model parameter, or a keyword argument of
# `fn` such as the screening interval `m`.
eval_at <- function(fn, x, arg, overrides, t) {
  swept <- setNames(list(x), arg)
  if (arg %in% names(params_def)) {
    overrides <- utils::modifyList(overrides, swept)
    swept <- list()
  }
  do.call(
    fn,
    c(list(params = modify_params(params_def, overrides), t = t), swept)
  )
}

#' Evaluate a model function over a grid, at default parameters and for every
#' bootstrap draw.
#'
#' Draws go through `prepare_draws()`, so Sym_B = SYM_B_RATIO * Sym_A in each
#' (finding C7). The point estimate uses the default parameters.
#'
#' @param fn     `life_add` or `screen_by_adherence`.
#' @param grid   Values of `arg`.
#' @param arg    Name of the swept parameter (or "m").
#' @param boot   Tibble of bootstrap parameter draws.
#' @param fixed  Named list of parameter overrides held fixed.
#' @param t      Follow-up months.
#' @return List with `point` (one result per grid value) and `draws` (per grid
#'   value, a list with one result per draw).
sweep_model <- function(fn, grid, arg, boot, fixed = list(), t = 120) {
  pd <- prepare_draws(boot)
  rows <- lapply(seq_len(nrow(pd)), function(i) lapply(pd, `[[`, i))
  list(
    point = lapply(grid, eval_at, fn = fn, arg = arg, overrides = fixed, t = t),
    draws = lapply(grid, function(x) {
      lapply(rows, function(r) {
        eval_at(fn, x, arg, utils::modifyList(r, fixed), t)
      })
    })
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
  ci <- vapply(
    res$draws,
    function(v) summarise_ci(unlist(v), probs),
    numeric(2)
  )
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
    mat <- do.call(rbind, lapply(res$draws[[k]], function(d) unlist(d[stages])))
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
  rows <- purrr::map(
    risk,
    ~ progress_untreated(n = n, t = years * 12, p_HCC = .x)
  )
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
    guides(
      colour = guide_legend("HCC Stage"),
      linetype = guide_legend("HCC Stage")
    ) +
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

# ----- Sensitivity (tornado) plots --------------------------------------------

#' Order the Parameter factor by the sum of absolute values of `by`.
#'
#' @param df  Table with a `Parameter` column.
#' @param by  Name of the column to order by.
#' @return `df` with Parameter reordered (largest effects at the top).
order_parameters <- function(df, by) {
  df$Parameter <- forcats::fct_reorder(df$Parameter, df[[by]], abssum)
  df
}

#' Add the facet group (`Type`) of each parameter, from `param_lookup`.
#'
#' @param df  Table with a `Parameter` column.
#' @return `df` with a factor column `Type` (levels = `names(group_labels)`).
add_group <- function(df) {
  df$Type <- factor(
    param_group(as.character(df$Parameter)),
    levels = names(group_labels)
  )
  df
}

#' Sensitivity table of the modifiable parameters only (those an intervention
#' could change), ordered and grouped for plotting.
#'
#' @param param_tbl  Output of `life_tibble()`.
#' @return Tibble `Parameter`, `Bound`, `years_mid_comp`, `Type`.
modifiable_sensitivity <- function(param_tbl) {
  keep <- param_tbl$years_mid_comp != 0 &
    !param_tbl$Parameter %in% NON_MODIFIABLE_PARAMS
  df <- add_group(order_parameters(param_tbl[keep, ], "years_mid_comp"))
  df[c("Parameter", "Bound", "years_mid_comp", "Type")]
}

#' Horizontal bar ("tornado") plot of life-years gained per parameter.
#'
#' @param data      Table with `Parameter` (ordered factor) and the column
#'                  named by `value`; plus `Type` for faceting.
#' @param value     Name of the column plotted on x.
#' @param y_labels  Axis labels of the parameters (`TeX_labels` or
#'                  `variable_labels`).
#' @param title,x_lab,y_lab  Labels.
#' @param facet     "none", "group" (one column of panels per `Type`) or "grid"
#'                  (Type by the `years_comp` column).
#' @return A ggplot.
plot_tornado <- function(
  data,
  value,
  y_labels = TeX_labels,
  title = NULL,
  x_lab = "Additional Life Years Gained From\nIntervention per Cohort of 100",
  y_lab = "Variable",
  facet = c("none", "group", "grid")
) {
  facet <- match.arg(facet)
  p <- ggplot(data, aes(y = .data$Parameter, x = .data[[value]])) +
    geom_col(position = "identity", colour = "black") +
    labs(x = x_lab, y = y_lab, title = title) +
    theme_hcc() +
    scale_y_discrete(labels = y_labels)
  switch(
    facet,
    none = p,
    group = p +
      facet_wrap(
        ~Type,
        ncol = 1,
        scales = "free_y",
        labeller = labeller(Type = facet_labels)
      ),
    grid = p +
      theme(legend.position = "bottom") +
      facet_grid(
        Type ~ years_comp,
        scales = "free",
        labeller = labeller(Type = facet_labels, years_comp = facet_labels)
      )
  )
}

#' Tornado plot for the lower or upper HCC-incidence bound of a sensitivity
#' table (p_HCC itself and parameters with no effect are dropped).
#'
#' @param param_tbl  Output of `life_tibble()`.
#' @param bound      "low" or "high" (column `years_<bound>_comp`).
#' @param title      Plot title.
#' @param exclude    Parameters to drop (default only p_HCC).
#' @param ...        Passed to `plot_tornado()`.
#' @return A ggplot.
plot_bound_tornado <- function(
  param_tbl,
  bound = c("low", "high"),
  title = NULL,
  exclude = "p_HCC",
  ...
) {
  bound <- match.arg(bound)
  col <- paste0("years_", bound, "_comp")
  df <- order_parameters(param_tbl, col)
  df <- df[df[[col]] != 0 & !df$Parameter %in% exclude, ]
  plot_tornado(df, col, title = title, ...)
}


#' Point estimates, 95% intervals and bootstrap draws of the life-year gain of
#' each improvement scenario over baseline, for one risk level.
#'
#' @param point  `run_scenarios()` output at default parameters.
#' @param boot   Wide bootstrap results (legacy columns) with the `_diff`
#'               columns of `scenario_columns()`.
#' @param risk   Risk level, a name of `risk_multipliers`.
#' @return List: `points` (Scenario, `Life Years per Cohort of 100`, type,
#'   lowCI, highCI, ledger_key) and `violin` (Scenario, `Life Years per Cohort
#'   of 100` for every bootstrap draw), both in the order treatment, diagnosis,
#'   optimistic adherence, realistic adherence.
scenario_diff_data <- function(point, boot, risk) {
  sfx <- ifelse(risk == "base", "", paste0("_", risk))
  tbl <- scenario_table[scenario_table$id != "baseline", ]
  tbl <- tbl[rev(seq_len(nrow(tbl))), ]
  boot_col <- paste0(tbl$column, "_diff", sfx)
  y <- "Life Years per Cohort of 100"
  points <- tibble::tibble(
    Scenario = tbl$label,
    gain = point$diff[point$risk == risk][match(
      tbl$id,
      point$scenario[point$risk == risk]
    )],
    type = "Point\nEstimate",
    lowCI = vapply(boot_col, \(cl) summarise_ci(boot[[cl]])[["lower"]], 1),
    highCI = vapply(boot_col, \(cl) summarise_ci(boot[[cl]])[["upper"]], 1),
    ledger_key = tbl$ledger_key
  )
  names(points)[names(points) == "gain"] <- y
  violin <- tibble::tibble(
    Scenario = rep(tbl$label, each = nrow(boot)),
    value = unlist(boot[boot_col], use.names = FALSE)
  )
  names(violin)[2] <- y
  list(points = points, violin = violin)
}

#' Violin plot of scenario gains with point estimates (Figs. in
#' figures/scenarios_new*.{jpg,png,pdf}).
#'
#' @param d     Output of `scenario_diff_data()`.
#' @param xlab  x-axis label.
#' @param ymax  Upper y limit.
plot_scenario_violins <- function(d, xlab, ymax) {
  y <- "Life Years per Cohort of 100"
  ggplot(d$violin, aes(.data[["Scenario"]], .data[[y]])) +
    geom_violin(trim = TRUE) +
    theme_hcc() +
    scale_x_discrete(
      limits = scenario_table$label[scenario_table$id != "baseline"]
    ) +
    labs(x = xlab) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)
    ) +
    geom_point(
      data = d$points,
      aes(.data[["Scenario"]], .data[[y]], shape = .data[["type"]]),
      size = 3,
      show.legend = TRUE
    ) +
    scale_y_continuous(limits = c(0, ymax), expand = expansion(0))
}

#' Absolute life-years of every scenario: bars at the point estimate and
#' violins of the bootstrap draws (figures/scenarios_all).
#'
#' @param point  `run_scenarios()` output at default parameters.
#' @param boot   Wide bootstrap results with the scenario columns of
#'               `scenario_table$column`.
#' @return A ggplot. Scenarios are ordered treatment, diagnosis, optimistic,
#'   realistic adherence, baseline.
plot_scenarios_all <- function(point, boot) {
  tbl <- scenario_table[rev(seq_len(nrow(scenario_table))), ]
  y <- "Life Years per Cohort of 100"
  at_base <- point[point$risk == "base", ]
  bars <- tibble::tibble(
    Scenario = tbl$label,
    gain = at_base$ly[match(tbl$id, at_base$scenario)]
  )
  names(bars)[2] <- y
  violin <- tibble::tibble(
    Scenario = rep(tbl$label, each = nrow(boot)),
    value = unlist(boot[tbl$column], use.names = FALSE)
  )
  names(violin)[2] <- y
  ggplot(bars, aes(forcats::fct_inorder(.data[["Scenario"]]), .data[[y]])) +
    geom_col() +
    geom_violin(
      data = violin,
      mapping = aes(.data[["Scenario"]], .data[[y]])
    ) +
    theme_hcc() +
    labs(x = "Scenario") +
    theme(legend.position = "none")
}

#' Life-years against the size of an improvement, one curve per parameter.
#'
#' @param specs  List of `list(name, param, max, sign)`: curve `name` is
#'               `life_add(<param> = default + sign * x)` for x = 0, 0.01, ..,
#'               max / 100.
#' @return Long tibble `x`, `value` (life-years), `group` (the curve name).
response_curves <- function(specs) {
  dplyr::bind_rows(lapply(specs, function(s) {
    x <- 0:s$max / 100
    vals <- vapply(
      x,
      function(d) {
        life_add(
          params = modify_params(
            params_def,
            setNames(list(params_def[[s$param]] + s$sign * d), s$param)
          )
        )
      },
      1
    )
    tibble::tibble(x = x, value = vals, group = s$name)
  }))
}

#' Line plot of `response_curves()`.
#'
#' @param curves     Output of `response_curves()`.
#' @param title      Plot title.
#' @param x_lab      x-axis label.
#' @param group_lab  Legend title.
#' @return A ggplot.
plot_response_curves <- function(curves, title, x_lab, group_lab) {
  ggplot(
    curves,
    aes(.data$x, .data$value, colour = .data$group)
  ) +
    geom_line(linewidth = 1.5) +
    ggtitle(title) +
    scale_x_continuous(labels = label_percent(), expand = expansion(0)) +
    labs(x = x_lab, y = "Life Extended", colour = group_lab) +
    theme_hcc(legend_title_centre = TRUE)
}

# ----- Bootstrap calibration figures (analysis/01) ----------------------------
# The five monthly progression probabilities p_0..p_D in long form.
rates_long <- function(pars) {
  pars[STAGE_P_NAMES] |> pivot_longer(everything())
}

#' Bootstrap point estimates and 95% limits for p_0..p_D.
#' @param pars   Bootstrap draws (columns p_0..p_D).
#' @param point  Named vector of point estimates.
rates_table <- function(pars, point) {
  ci <- vapply(pars[STAGE_P_NAMES], summarise_ci, numeric(2))
  tibble(
    rate = factor(STAGE_P_NAMES),
    Point = unname(point[STAGE_P_NAMES]),
    Lower = unname(ci["lower", ]),
    Upper = unname(ci["upper", ])
  )
}

#' Counts of each bootstrapped calibration input (m_*, D_*), faceted.
plot_bootstrap_counts <- function(long) {
  long$name <- vapply(long$name, \(n) paste(Boot_labels[n]), character(1))
  ggplot(long, aes(x = value)) +
    geom_bar() +
    facet_wrap(vars(name), scales = "free", labeller = label_parsed, nrow = 3) +
    theme_hcc() +
    scale_y_continuous(expand = expansion(c(0, 0.05)))
}

#' Box plot or violin plot of the bootstrap rates, red dot = point estimate.
plot_rates_distribution <- function(
  pars,
  rates_tbl,
  type = c("box", "violin")
) {
  type <- match.arg(type)
  geom <- if (type == "box") {
    geom_boxplot(outlier.shape = 1, outlier.alpha = 0.2)
  } else {
    geom_violin(scale = "width")
  }
  labs_x <- if (type == "box") {
    c(
      p_0 = "BCLC Stage\n0 to A HCC",
      p_A = "BCLC Stage\nA to B HCC",
      p_B = "BCLC Stage\nB to C HCC",
      p_C = "BCLC Stage\nC to D HCC",
      p_D = "BCLC Stage\nD HCC to Death"
    )
  } else {
    Boot_labels
  }
  ggplot(rates_long(pars), aes(fct_inorder(name), value)) +
    geom +
    geom_point(
      aes(x = rate, y = Point),
      data = rates_tbl,
      inherit.aes = FALSE,
      colour = "red"
    ) +
    scale_x_discrete(labels = labs_x) +
    labs(
      x = NULL,
      y = if (type == "box") "Probability of\nProgressing per Month" else NULL
    ) +
    theme_hcc() +
    scale_y_continuous(limits = c(0, 1), expand = expansion(0))
}

#' Histograms of the bootstrap rates, red line = point estimate.
plot_rates_hist <- function(pars, rates_tbl) {
  lab <- \(x) vapply(as.character(x), \(n) paste(Boot_labels[n]), character(1))
  long <- dplyr::rename(rates_long(pars), rate = name)
  long$rate <- lab(long$rate)
  rates_tbl$rate <- lab(rates_tbl$rate)
  ggplot(long, aes(value)) +
    geom_histogram(bins = 20) +
    geom_vline(aes(xintercept = Point), data = rates_tbl, colour = "red") +
    facet_wrap(vars(rate), nrow = 3, scales = "free", labeller = label_parsed) +
    theme_hcc() +
    scale_y_continuous(expand = expansion(c(0, 0.05)))
}

#' 95% interval of the mean years in each BCLC stage (1 / p / 12).
plot_time_range <- function(pars, rates_tbl) {
  ci <- rates_long(pars) |>
    mutate(mean_years = 1 / value / 12, name = fct_inorder(name)) |>
    group_by(name) |>
    summarise(
      mean025 = quantile(mean_years, 0.025),
      mean975 = quantile(mean_years, 0.975)
    )
  ggplot(ci, aes(name, ymin = mean025, ymax = mean975)) +
    geom_errorbar() +
    geom_point(
      aes(x = rate, y = 1 / Point / 12),
      data = rates_tbl,
      inherit.aes = FALSE
    ) +
    scale_x_discrete(
      labels = setNames(
        paste("BCLC Stage", c(0, "A", "B", "C", "D")),
        STAGE_P_NAMES
      )
    ) +
    labs(x = NULL, y = "Mean Number of\nYears In BCLC Stage") +
    theme_hcc() +
    scale_y_continuous(limits = c(0, 3), expand = expansion(0))
}

#' Hexbin of the ultrasound sensitivities in BCLC 0 and A; X = point estimate.
plot_sensitivity_hex <- function(pars, point) {
  pct <- scale_x_continuous(
    labels = label_percent(),
    expand = expansion(0),
    limits = c(0, 1)
  )
  ggplot(pars, aes(x = Sens_0, y = Sens_A)) +
    geom_area(
      aes(
        x = seq(0, 1, length.out = nrow(pars)),
        y = seq(0, 1, length.out = nrow(pars))
      ),
      alpha = 0.25
    ) +
    geom_hex(alpha = 0.75, bins = c(100, 80)) +
    theme_hcc(legend_title_centre = TRUE) +
    scale_fill_viridis(
      option = "C",
      trans = "log",
      breaks = c(1, 2, 5, 10, 20, 35)
    ) +
    pct +
    scale_y_continuous(
      labels = label_percent(),
      expand = expansion(0),
      limits = c(0, 1)
    ) +
    geom_point(
      aes(x = point[["Sens_0"]], y = point[["Sens_A"]]),
      shape = 4,
      stroke = 1,
      size = 1,
      colour = "black"
    ) +
    labs(
      x = "BCLC Stage 0 HCC\nUltrasound Sensitivity",
      y = "BCLC Stage A HCC\nUltrasound Sensitivity",
      fill = "Number of\nSamples"
    )
}

#' Histogram of bootstrap life years gained with a fitted normal; red line =
#' the point estimate `life_add()`.
plot_baseline_hist <- function(pars) {
  ggplot(pars, aes(x = Baseline)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.5) +
    geom_vline(xintercept = life_add(), colour = "red") +
    theme_hcc() +
    scale_y_continuous(expand = expansion(0), limits = c(0, 0.08)) +
    geom_function(
      fun = dnorm,
      args = list(mean = mean(pars$Baseline), sd = sd(pars$Baseline)),
      colour = "blue",
      linewidth = 1
    ) +
    labs(x = "Life Years Gained", y = "Density")
}

#' Partial-correlation heat map (square matrix `pcor`) and the bar chart of
#' the Baseline row against the other variables.
plot_pcor_heatmap <- function(pcor) {
  vars <- colnames(pcor)
  long <- as_tibble(pcor) |>
    mutate(row = fct_inorder(vars)) |>
    pivot_longer(-row, names_to = "col") |>
    mutate(col = fct_inorder(col))
  ggplot(long, aes(x = col, y = row, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(limit = c(-1, 1), low = "red", high = "#075AFF") +
    scale_y_discrete(expand = expansion(0), limits = rev) +
    scale_x_discrete(expand = expansion(0), position = "top") +
    theme_hcc() +
    geom_text(aes(label = round(value, 2)), color = "black") +
    theme(axis.text.x = element_text(angle = 45, hjust = 0)) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_colourbar(title = "Partial Rank\nCoefficient"))
}

plot_pcor_baseline <- function(pcor) {
  d <- tibble(Variable = fct_inorder(setdiff(colnames(pcor), "Baseline")))
  d$value <- unname(pcor["Baseline", as.character(d$Variable)])
  ggplot(d, aes(x = Variable, y = value, fill = value)) +
    geom_col(colour = "black") +
    scale_fill_gradient2(limit = c(-1, 1), low = "red", high = "#075AFF") +
    theme_hcc() +
    labs(y = "Partial Correlation Coefficient", x = NULL) +
    guides(fill = guide_colourbar(title = NULL)) +
    scale_y_continuous(limits = c(-1, 1), expand = expansion(0)) +
    scale_x_discrete(labels = TeX_labels)
}
