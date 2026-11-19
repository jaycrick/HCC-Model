# Calibration figures (file R/plot_calibration.R)
# Bootstrap calibration figures (analysis/01).
# Requires R/setup.R, R/bootstrap.R and R/plot_theme.R.

# ----- Bootstrap calibration figures (analysis/01) ----------------------------
#' The five monthly progression probabilities p_0..p_D in long form.
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
    # Ref: Cumming 2024, Sec. 3.1 / Fig. 2A (mean time in stage = 1 / p
    # months, shown in years)
    mutate(mean_years = 1 / value / 12, name = fct_inorder(name)) |>
    group_by(name) |>
    summarise(
      mean025 = quantile(mean_years, 0.025),
      mean975 = quantile(mean_years, 0.975)
    )
  ggplot(ci, aes(name, ymin = mean025, ymax = mean975)) +
    geom_errorbar() +
    geom_point(
      # Ref: Cumming 2024, Fig. 2A (point estimate, 1 / p / 12 years)
      aes(x = rate, y = 1 / Point / 12),
      data = rates_tbl,
      inherit.aes = FALSE
    ) +
    scale_x_discrete(
      labels = setNames(
        paste("BCLC Stage", STAGE_IDS),
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
    # Finding P4: Supp. Fig. S3 caption says N(mu = 68.0, sd = 5.7). Here the
    # normal uses the bootstrap mean (about 67.6); 68.0 is life_add(), the
    # red line. Ref: Supplement Fig. S3 (p. 10).
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

#' Bar chart of the partial correlations of the Baseline row.
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
