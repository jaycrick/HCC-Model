# Scenario figures (file R/plot_scenarios.R)
# Sensitivity (tornado) plots, clinical-scenario figures and response curves
# (analysis/02, analysis/04).
# Requires R/setup.R, R/model_core.R, R/scenarios.R and R/plot_theme.R.

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
    scale_x_discrete(limits = scenario_table$label[scenario_table$id != "baseline"]) +
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
    geom_violin(data = violin, mapping = aes(.data[["Scenario"]], .data[[y]])) +
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
    # percentage points to a proportion; Ref: Cumming 2024, Sec. 2.4
    # (absolute percentage-point changes)
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
  ggplot(curves, aes(.data$x, .data$value, colour = .data$group)) +
    geom_line(linewidth = 1.5) +
    ggtitle(title) +
    scale_x_continuous(labels = label_percent(), expand = expansion(0)) +
    labs(x = x_lab, y = "Life Extended", colour = group_lab) +
    theme_hcc(legend_title_centre = TRUE)
}
