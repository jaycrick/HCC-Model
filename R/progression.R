# R/progression.R
# Analysis-and-figure functions for the progression and screening-frequency
# components of the model.  Each function computes a summary tibble, builds
# one or more ggplots, saves them to figures/, and returns the primary plot.
# Requires: R/setup.R, R/model_core.R, and resampled_pars in the environment.

# ----- Randomised PSA parameter draws ----------------------------------------

#' Draw n randomised parameter sets using PERT distributions.
#'
#' @param x       Number of draws.
#' @param ignore  Character vector of parameter names to drop (e.g. "p_HCC").
#' @return List of x one-row tibbles, each a full parameter set.
ran_vars <- function(x = 1, ignore = NULL) {
  rtbls <- vector("list", length = x)

  for (i in 1:x) {
    rparam_tbl <- as_tibble(param_def) %>%
      rowwise() %>%
      mutate(new = rpert(1, LOWER, POINT, UPPER)) %>%
      ungroup() %>%
      bind_cols(
        Stage = c(
          "p_HCC",
          "p_0",
          "p_A",
          "p_B",
          "p_C",
          "p_D",
          "Sens_0",
          "Sens_A",
          "Sens_B",
          "Sens_C",
          "Sens_D",
          "p_NA",
          "p_CA",
          "p_Att",
          "L_0",
          "L_A",
          "L_B",
          "L_C",
          "L_D",
          "Sym_0",
          "Sym_A",
          "Sym_B",
          "Sym_C",
          "Sym_D"
        )
      )

    rtbls[[i]] <- dfs(
      Parameter = rparam_tbl$Stage,
      value = rparam_tbl$new
    ) %>%
      select(!all_of(ignore))
  }
  return(rtbls)
}


# ----- Distribution of BCLC stages across a range of HCC risk ----------------

#' Plot HCC stage distribution over a range of annual HCC risk levels.
#'
#' @param y  Horizon in years.
#' @return ggplot object; also saves to figures/.
risk_no_int <- function(n, y) {
  risk <- seq(0, 0.08, length.out = 1000)

  my_tbl <- map_dfr(risk, ~ fromf4(n = 100, t = y * 12, p_HCC = .x)) %>%
    tibble(`Annual HCC Risk Profile` = risk, .) %>%
    pivot_longer(
      c(Cirrhosis:Mortality),
      names_to = "HCC_stage",
      values_to = "percent"
    ) %>%
    mutate(
      HCC_stage = factor(HCC_stage, levels = stage_list),
      graph_location = factor(
        ifelse(HCC_stage %in% c("Cirrhosis", "Mortality"), "top", "bottom"),
        levels = c("top", "bottom")
      )
    )

  my_plot <- ggplot(
    my_tbl,
    aes(
      x = `Annual HCC Risk Profile`,
      y = percent,
      group = HCC_stage,
      colour = HCC_stage,
      linetype = HCC_stage
    )
  ) +
    geom_line(linewidth = 1.5) +
    theme_bw(base_family = "serif") +
    my_scale +
    my_line +
    scale_x_continuous(
      labels = label_percent(),
      expand = expansion(0)
    ) +
    scale_y_continuous(
      labels = label_percent(scale = 1),
      expand = expansion(c(0, 0.01))
    ) +
    guides(
      colour = guide_legend("HCC Stage"),
      linetype = guide_legend("HCC Stage")
    ) +
    facet_wrap(vars(graph_location), dir = "v", scales = "free_y") +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      strip.text.x = element_blank()
    ) +
    labs(
      x = "Annual HCC Risk Profile",
      y = "Percentage of Population per Stage",
      title = paste(
        "Model-estimated Distribution of BCLC\nStages From Untreated HCC After",
        y,
        "Years"
      ),
      caption = "With the assumption no-one has\nHCC at the start of the model"
    )

  ggsave(
    filename = here(
      "figures",
      paste0("Untreated_HCp_Cistribution_by_Risk_", y, "_Year.jpg")
    ),
    plot = my_plot,
    width = 5.5,
    height = 6
  )

  return(my_plot)
}


# ----- BCLC stage distribution over time -------------------------------------

#' Plot HCC stage distribution over time at default HCC risk.
#'
#' @param y  Horizon in years.
#' @return ggplot object; also saves to figures/.
years_no_int <- function(y) {
  months <- 0:(y * 12)

  my_tbl <- map_dfr(months, ~ fromf4(n = 100, t = .x)) %>%
    tibble(Years = months / 12, .) %>%
    pivot_longer(
      c(Cirrhosis:Mortality),
      names_to = "HCC_stage",
      values_to = "percent"
    ) %>%
    mutate(
      HCC_stage = factor(HCC_stage, levels = stage_list),
      graph_location = factor(
        ifelse(HCC_stage %in% c("Cirrhosis", "Mortality"), "top", "bottom"),
        levels = c("top", "bottom")
      )
    )

  my_plot <- ggplot(
    my_tbl,
    aes(
      x = Years,
      y = percent,
      group = HCC_stage,
      colour = HCC_stage,
      linetype = HCC_stage
    )
  ) +
    geom_line(linewidth = 1.5) +
    theme_bw(base_family = "serif") +
    my_scale +
    my_line +
    scale_x_continuous(expand = expansion(0)) +
    scale_y_continuous(
      labels = label_percent(scale = 1),
      expand = expansion(c(0, 0.01))
    ) +
    guides(
      colour = guide_legend("HCC Stage"),
      linetype = guide_legend("HCC Stage")
    ) +
    facet_wrap(vars(graph_location), scales = "free_y", dir = "v") +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      strip.text.x = element_blank()
    ) +
    labs(
      x = "Years Since Model Start",
      y = "Percentage of Population per Stage",
      title = paste(
        "Model-estimated Distribution of BCLC\nStages from Untreated HCC Over",
        y,
        "Years"
      ),
      caption = paste(
        "With the assumption no-one has\nHCC at the start of the model and a",
        p_HCC_def,
        "annual HCC risk"
      )
    )

  ggsave(
    filename = here(
      "figures",
      paste0("Untreated_HCC_over_time_", y, "_Year.jpg")
    ),
    plot = my_plot,
    width = 5.5,
    height = 6
  )

  return(my_plot)
}


# ----- Effect of ultrasound screening frequency ------------------------------

#' Life years gained vs ultrasound frequency with bootstrap uncertainty band.
#'
#' Requires `resampled_pars` in the calling environment.
#'
#' @param y  Horizon in years.
#' @return ggplot object (linear x-axis); also saves log- and linear-scale
#'   versions to figures/.
usfreq_graph <- function(y) {
  usfreq <- 1:24
  tempfun <- function(m) pmap_dbl(resampled_pars, m = m, t = y * 12, life_add)
  testyll <- map_dbl(usfreq, ~ life_add(m = .x, t = y * 12))
  rand_list <- map(usfreq, tempfun)
  saveddf <- tibble(usfreq, testyll, rand_list) %>%
    mutate(
      min = map_dbl(rand_list, ~ quantile(.x, probs = 0.025)),
      max = map_dbl(rand_list, ~ quantile(.x, probs = 0.975))
    )

  yat_base <- ggplot(saveddf, aes(usfreq, testyll)) +
    geom_line(linewidth = 1.25) +
    geom_ribbon(aes(ymin = min, ymax = max, fill = "50%"), alpha = 0.2) +
    theme_bw(base_family = "serif") +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_text(hjust = 0.5)
    ) +
    scale_y_continuous(expand = expansion(0)) +
    labs(
      title = paste(
        "Effect of Ultrasound Screening Frequency\non Life Years over",
        y,
        "Years"
      ),
      y = "Additional Life Years Gained From\nIntervention per Cohort of 100",
      fill = "Confidence\nBound"
    ) +
    scale_fill_manual(values = "black")

  yat_log <- yat_base +
    scale_x_log10(
      breaks = round(10^(seq(0, log10(24), length.out = 7))),
      expand = expansion(0)
    ) +
    xlab("Ultrasound Frequency (log scale)")

  yat <- yat_base +
    xlab("Ultrasound Frequency") +
    scale_x_continuous(expand = expansion(0))

  ggsave(
    filename = here(
      "figures",
      paste0("HCC_Years_of_Life_Saved_by_Testing_Frequency_", y, "_Year.jpg")
    ),
    plot = yat,
    width = 5.5,
    height = 6
  )
  ggsave(
    filename = here(
      "figures",
      paste0(
        "HCC_Years_of_Life_Saved_by_Testing_Frequency_",
        y,
        "_log_Year.jpg"
      )
    ),
    plot = yat_log,
    width = 5.5,
    height = 6
  )

  return(yat)
}


# ----- Effect of complete-adherence proportion --------------------------------

#' Life years gained vs proportion in completely-adherent group.
#'
#' Requires `resampled_pars` in the calling environment.
#'
#' @param y  Horizon in years.
#' @return ggplot object; also saves to figures/.
pafunc <- function(y) {
  probs <- seq(0, 1 - p_NA_def, length.out = 50)
  tempfun <- function(p)
    pmap_dbl(resampled_pars, p_CA = p, t = y * 12, life_add)
  testyll <- map_dbl(probs, ~ life_add(p_CA = .x, t = y * 12))
  rand_list <- map(probs, tempfun)
  saveddf <- tibble(probs, testyll, rand_list) %>%
    mutate(
      min = map_dbl(rand_list, ~ quantile(.x, probs = 0.025)),
      max = map_dbl(rand_list, ~ quantile(.x, probs = 0.975))
    )

  yatpa <- ggplot(saveddf, aes(probs, testyll)) +
    geom_line(linewidth = 1.25) +
    geom_ribbon(aes(ymin = min, ymax = max, fill = "50%"), alpha = 0.2) +
    theme_bw(base_family = "serif") +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_text(hjust = 0.5)
    ) +
    scale_fill_manual(values = "black") +
    labs(
      title = paste(
        "Effect of Complete Adherence\non Life Years over",
        y,
        "Years"
      ),
      x = "Population in Completely Adherent Cohort",
      y = "Additional Life Years Gained From\nIntervention per Cohort of 100",
      caption = "People are taken from the partially\nadherent cohort, whilst the non-adherent\ncohort is assumed to be stationary",
      fill = "Confidence\nBound"
    ) +
    scale_x_continuous(labels = label_percent(), expand = expansion(0)) +
    scale_y_continuous(expand = expansion(0))

  ggsave(
    filename = here(
      "figures",
      paste0("HCC_Years_of_Life_Saved_by_p_CA_", y, "_Year.jpg")
    ),
    plot = yatpa,
    width = 5.5,
    height = 9
  )

  return(yatpa)
}


# ----- Stage-specific diagnoses vs attendance probability --------------------

#' Diagnoses per stage vs screening attendance probability in the PA group.
#'
#' Requires `resampled_pars` in the calling environment.
#'
#' @param y  Horizon in years.
#' @return ggplot object; also saves ribbon and non-ribbon versions to figures/.
p_Att_Stages <- function(y) {
  p_Att <- (0:20) / 20
  tempfun <- function(p)
    map_df(
      resampled_pars,
      ~ pmap_df(.x, p_NA = 0, p_CA = 0, p_Att = p, t = y * 12, adh_theor)
    )
  number <- map_df(
    p_Att,
    ~ adh_theor(p_Att = .x, p_NA = 0, p_CA = 0, t = y * 12)
  ) %>%
    select(-c(Cirrhosis, Mortality)) %>%
    pivot_longer(`BCLC 0`:`BCLC D`, values_to = "number") %>%
    select(number)
  rand_list <- map_df(p_Att, ~ tibble(p_Att = .x, tempfun(.x)))
  saveddf <- rand_list %>%
    select(-c(Cirrhosis, Mortality)) %>%
    pivot_longer(`BCLC 0`:`BCLC D`, names_to = "HCC_Stage") %>%
    group_by(p_Att, HCC_Stage) %>%
    summarise(
      min = quantile(value, probs = 0.25),
      max = quantile(value, probs = 0.75)
    ) %>%
    tibble(., number)

  plot_att_stages <- ggplot(
    saveddf,
    aes(
      p_Att,
      number,
      group = HCC_Stage,
      colour = HCC_Stage,
      linetype = HCC_Stage,
      fill = HCC_Stage
    )
  ) +
    geom_line(linewidth = 1.5) +
    theme_bw(base_family = "serif") +
    ggtitle(
      "Effect of Screening Attendance within the\nPartial Adherence Cohort on Diagnosis of HCC"
    ) +
    scale_x_continuous(labels = label_percent(), expand = expansion(0)) +
    my_line +
    my_scale +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_text(hjust = 0.5)
    ) +
    labs(
      fill = "Stage of\nDiagnosis",
      caption = paste("In 100 people over", y, "years"),
      x = "Probability of Attending\nScheduled Ultrasound Screenings",
      y = "Number Diagnosed by Stage"
    )

  ggsave(
    here("figures", paste0("attendance_plot_stages_", y, "_year.jpg")),
    plot_att_stages,
    width = 5,
    height = 8
  )

  plot_att_stages_ribbon <- ggplot(
    saveddf,
    aes(
      p_Att,
      number,
      group = HCC_Stage,
      colour = HCC_Stage,
      linetype = HCC_Stage,
      fill = HCC_Stage
    )
  ) +
    geom_line(linewidth = 1.5) +
    theme_bw(base_family = "serif") +
    ggtitle(
      "Effect of Screening Attendance\non Stage of Diagnosis within\nthe Partial Adherence Cohort"
    ) +
    scale_x_continuous(labels = label_percent(), expand = expansion(0)) +
    my_fill +
    geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.2) +
    my_line +
    my_scale +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_text(hjust = 0.5)
    ) +
    labs(
      fill = "Stage of\nDiagnosis",
      caption = paste("In 100 people over", y, "years"),
      x = "Probability of Attending\nScheduled Ultrasound Screenings",
      y = "Number Diagnosed by Stage"
    )

  ggsave(
    here("figures", paste0("attendance_plot_stages_ribbon_", y, "_year.jpg")),
    plot_att_stages_ribbon,
    width = 5,
    height = 8
  )

  return(plot_att_stages)
}
