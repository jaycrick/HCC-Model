# Theme and saving (file R/plot_theme.R)
# Figure helpers shared by every figure: the ggplot theme, save_figure(), and
# the interval helpers (summarise_ci(), band_label()).  Nothing in the plot_*.R
# files reads files or sets seeds; the analysis scripts pass the data (e.g.
# bootstrap draws) in explicitly.  Requires R/setup.R.

# ----- Theme and saving -------------------------------------------------------

# Figure formats and the dpi each is saved at when a figure is written in all
# three (jpg 300, png and pdf 600).
FIGURE_FORMATS <- c("jpg", "png", "pdf")
FIGURE_DPI_ALL <- c(jpg = 300, png = 600, pdf = 600)

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
  formats = FIGURE_FORMATS,
  dpi = 600,
  dir = here::here("figures")
) {
  formats <- match.arg(formats, FIGURE_FORMATS, several.ok = TRUE)
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
