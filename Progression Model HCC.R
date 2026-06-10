# Init. ----------------------------
load("Sensitivity_analyses.RData")

library(conflicted)
library(expm)
library(beepr)
library(tidyverse)
library(here)
library(RColorBrewer)
library(ggthemes)
library(ggsci)
library(gridExtra)
library(scales)
library(mc2d)
library(magrittr)
library(viridis)

set.seed(4567385)

my_pal <- pal_igv("default")(7)
names(my_pal) <- c(
  "Cirrhosis",
  "Mortality",
  "BCLC 0",
  "BCLC A",
  "BCLC D",
  "BCLC B",
  "BCLC C"
)
my_linet <- as.vector(c(
  "solid",
  "solid",
  "solid",
  "solid",
  "solid",
  "solid",
  "dotdash"
))
names(my_linet) <- stage_list

my_scale <- scale_colour_manual(
  name = "HCC_Stage",
  values = my_pal,
  limits = stage_list
)
my_fill <- scale_fill_manual(
  name = "HCC_Stage",
  values = my_pal,
  limits = stage_list
)
my_line <- scale_linetype_manual(
  name = "HCC_Stage",
  values = my_linet,
  limits = stage_list
)

turbo_pal <- turbo(7)
names(turbo_pal) <- stage_list
my_scale <- scale_colour_manual(
  name = "HCC_Stage",
  values = turbo_pal,
  limits = stage_list
)
my_fill <- scale_fill_manual(
  name = "HCC_Stage",
  values = turbo_pal,
  limits = stage_list
)


# Initialising some randomised variables ------------
ran_vars <- function(x = 1, ignore = NULL) {
  # x = number of randomisations, with a fixed p_HCC

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

    rtbls[[i]] <- dfs(Parameter = rparam_tbl$Stage, value = rparam_tbl$new) %>%
      select(!all_of(ignore))
  }
  return(rtbls)
}

random_vars <- ran_vars(x = 100)
resampled_pars <- readRDS(file = "resampled_pars_fin.rds")

# Setting parameters--------------

risk_no_int <- function(n, y) {
  # Graphs no intervention over y years
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
    # setting up aes etc.
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
      # strip.background = element_blank(),
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
      "Figures",
      paste0("Untreated_HCp_Cistribution_by_Risk_", y, "_Year.jpg")
    ),
    plot = my_plot,
    width = 5.5,
    height = 6
  )

  return(my_plot)
}

plot_untreated_risk_10_years <- risk_no_int(100, 10)

# plot_untreated_risk_20_years <- risk_no_int(100, 20)

# plot_untreated_risk_5_years <- risk_no_int(100, 5)

# GRAPHING PROGRESSION OVER MULTIPLE MONTHS ---------------------------------

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
    # setting up aes etc.
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
    scale_x_continuous(
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
    facet_wrap(vars(graph_location), scales = "free_y", dir = "v") +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      # strip.background = element_blank(),
      strip.text.x = element_blank()
      # aspect.ratio = 3
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
      "Figures",
      paste0("Untreated_HCC_over_time_", y, "_Year.jpg")
    ),
    plot = my_plot,
    width = 5.5,
    height = 6
  )

  return(my_plot)
}

# plot_untreated_prog_20_years <- years_no_int(20)

plot_untreated_prog_10_years <- years_no_int(10)

# plot_untreated_prog_5_years <- years_no_int(5)

# US FREQ GRAPHS ---------------------------

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
    scale_y_continuous(
      expand = expansion(0)
    ) +
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
      breaks = round(
        10^(seq(0, log10(24), length.out = 7))
      ),
      expand = expansion(0)
    ) +
    xlab("Ultrasound Frequency (log scale)")

  yat <- yat_base +
    xlab("Ultrasound Frequency") +
    scale_x_continuous(
      expand = expansion(0)
    )

  ggsave(
    filename = here(
      "Figures",
      paste0("HCC_Years_of_Life_Saved_by_Testing_Frequency_", y, "_Year.jpg")
    ),
    plot = yat,
    width = 5.5,
    height = 6
  )
  ggsave(
    filename = here(
      "Figures",
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

plot_usfreq_10_years <- usfreq_graph(10)
# plot_usfreq_20_years <- usfreq_graph(20)

# GRAPHS BY LIKELIHOOD OF ATTENDED IF IN PA -----------

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
    scale_x_continuous(
      labels = label_percent(),
      expand = expansion(0)
    ) +
    scale_y_continuous(
      expand = expansion(0)
    )

  ggsave(
    filename = here(
      "Figures",
      paste0("HCC_Years_of_Life_Saved_by_p_CA_", y, "_Year.jpg")
    ),
    plot = yatpa,
    width = 5.5,
    height = 9
  )

  return(yatpa)
}

# plot_partial_adherence_5_years <- pafunc(5)
plot_partial_adherence_10_years <- pafunc(10)
# plot_partial_adherence_20_years <- pafunc(20)

# Randomising some values and generating a range ----------------------------------------------------

life_tbl <- ran_vars(100, ignore = c("p_HCC")) %>%
  tibble(dfs = ., life = map_dbl(., years))

ggplot(life_tbl, aes(y = life)) +
  geom_boxplot(fill = "lightblue") +
  theme_bw(base_family = "serif") +
  geom_hline(yintercept = life_add(), colour = "red", size = 1)

plot_sim <- ggplot(life_tbl, aes(x = version, y = years_saved)) +
  geom_boxplot() +
  theme_bw(base_family = "serif")

ggsave(
  here("Figures", "simulation_range.jpg"),
  plot_sim,
  width = 5.5,
  height = 9
)


# Varying US Attendance for partial adh -----------

p_Att <- (0:20) / 20
tempfun <- function(p)
  pmap_dbl(resampled_pars, p_NA = 0, p_CA = 0, p_Att = p, life_add)
testyll <- map_dbl(p_Att, ~ life_add(p_Att = .x, p_NA = 0, p_CA = 0))
rand_list <- map(p_Att, tempfun)
saveddf <- tibble(p_Att, testyll, rand_list) %>%
  mutate(
    min = map_dbl(rand_list, ~ quantile(.x, probs = 0.025)),
    max = map_dbl(rand_list, ~ quantile(.x, probs = 0.975))
  )

plot_att <- ggplot(saveddf, aes(p_Att, testyll)) +
  geom_line(linewidth = 1.5) +
  theme_bw(base_family = "serif") +
  ggtitle(NULL) +
  scale_x_continuous(
    labels = label_percent(),
    expand = expansion(0)
  ) +
  geom_ribbon(aes(ymin = min, ymax = max, fill = "95%"), alpha = 0.2) +
  scale_fill_manual(values = "black") +
  labs(
    fill = "Confidence\nBand",
    x = "Probability of Attending\nScheduled Ultrasound Screenings",
    y = "Additional Life Years Gained From\nIntervention per Cohort of 100",
    # caption = "In 100 people over 10 years\nwith 100% partial adherence"
  ) +
  # scale_x_log10(
  #   labels = label_percent(),
  #   expand = expansion(0)
  # ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  )

plot_att

ggsave(
  here("Figures", "attendance_plot_10_year.jpg"),
  plot_att,
  width = 6,
  height = 5.5
)
ggsave(
  here("Figures", "attendance_plot_10_year.png"),
  plot_att,
  width = 6,
  height = 5.5,
  dpi = 600
)
ggsave(
  here("Figures", "attendance_plot_10_year.pdf"),
  plot_att,
  width = 6,
  height = 5.5,
  dpi = 600
)


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
    scale_x_continuous(
      labels = label_percent(),
      expand = expansion(0)
    ) +
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
    here("Figures", paste0("attendance_plot_stages_", y, "_year.jpg")),
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
    scale_x_continuous(
      labels = label_percent(),
      expand = expansion(0)
    ) +
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
    here("Figures", paste0("attendance_plot_stages_ribbon_", y, "_year.jpg")),
    plot_att_stages_ribbon,
    width = 5,
    height = 8
  )

  return(plot_att_stages)
}

# plot_att_stages_20_years <- p_Att_Stages(20)
plot_att_stages_10_years <- p_Att_Stages(10)
# plot_att_stages_5_years <- p_Att_Stages(5)

all_objects_in_workspace <- data.frame(list = ls())

plotList <- all_objects_in_workspace %>%
  filter(str_detect(list, "plot"))

save(list = plotList$list, file = "plots.RData")

# SIMULATION----

hccsim <- function(n, m) {
  # n=number of people, m=number of months
  vec <- c(n, rep(0, 6)) # initialise the vector, with all people in F4
  for (i in 1:m) {
    # for each month
    for (j in (length(vec) - 1):1) {
      worse <- rbinom(1, vec[j], tmat[j, j + 1]) # the number of people who progress
      vec[j + 1] <- vec[j + 1] + worse
      vec[j] <- vec[j] - worse
    }
  }
  return(vec)
}

# hccsim(100, 24)

simcomp <- function(n, l, s) hccsim(n, l) - hccsim(n, s) # compares testing increments with a simulation

# simcomp(13000, 6, 4)

# Now with ultrasound----------------------------

us_sim <- function(n, t, m) {
  # n=number of people, t=total time in months, m=number of months checkup
  vec <- c(n, rep(0, 6)) # initialise the vector, with all people in F4
  fin <- rep(0, 7)
  for (i in 1:t) {
    # for each month
    for (j in (length(vec) - 1):1) {
      worse <- rbinom(1, vec[j], tmat[j, j + 1]) # the number of people who progress
      vec[j + 1] <- vec[j + 1] + worse
      vec[j] <- vec[j] - worse
    }
    sym <- c(
      0,
      rbinom(1, vec[2], Sym_0),
      rbinom(1, vec[3], Sym_A),
      rbinom(1, vec[4], Sym_B),
      rbinom(1, vec[5], Sym_C),
      rbinom(1, vec[6], Sym_D),
      vec[7]
    )
    fin <- fin + sym
    vec <- vec - sym
    if (i %% m == 0) {
      us <- c(
        0,
        rbinom(1, vec[2], Sens_0),
        rbinom(1, vec[3], Sens_A),
        rbinom(1, vec[4], Sens_B),
        rbinom(1, vec[5], Sens_C),
        rbinom(1, vec[6], Sens_D),
        vec[7]
      ) # how many are picked up at the ultrasound
      fin <- fin + us
      vec <- vec - us
    }
  }
  return(fin)
}

# us_sim(13000, 36, 6)
us_sim(13000, 6, 6)

# WITH ADHERRENCE GROUPS -------------------------

us_sim_pa <- function(n, t, m) {
  # now just for the pa group, same as us_sim except for comments
  vec <- c(n, rep(0, 6))
  fin <- rep(0, 7)
  for (i in 1:t) {
    for (j in (length(vec) - 1):1) {
      worse <- rbinom(1, vec[j], tmat[j, j + 1])
      vec[j + 1] <- vec[j + 1] + worse
      vec[j] <- vec[j] - worse
      if (i %% m == 0) {
        sym <- c(
          0,
          rbinom(1, vec[2], Sym_0),
          rbinom(1, vec[3], Sym_A),
          rbinom(1, vec[4], Sym_B),
          rbinom(1, vec[5], Sym_C),
          rbinom(1, vec[6], Sym_D),
          vec[7]
        )
        fin <- fin + sym
        vec <- vec - sym
        go <- c(
          0,
          rbinom(1, vec[2], pr_att),
          rbinom(1, vec[3], pr_att),
          rbinom(1, vec[4], pr_att),
          rbinom(1, vec[5], pr_att),
          rbinom(1, vec[6], pr_att),
          vec[7]
        ) # how many actually go to the ultrasound
        us <- c(
          0,
          rbinom(1, go[2], Sens_0),
          rbinom(1, go[3], Sens_A),
          rbinom(1, go[4], Sens_B),
          rbinom(1, go[5], Sens_C),
          rbinom(1, go[6], Sens_D),
          vec[7]
        ) # how many are picked up at the ultrasound
        fin <- fin + us
        vec <- vec - us
      }
    }
  }
  return(fin)
}

adh_sim <- function(n, t, m) {
  # number of people, total time tracking, interval of US
  adh_group <- rmultinom(1, n, c(pr_na, pr_pa, 1 - pr_na - pr_pa)) # number of people in each adherence group
  stage_caught <- us_sim(adh_group[3], t, m) + us_sim_pa(adh_group[2], t, m)
  stage_caught[7] <- stage_caught[7] + adh_group[1]
  return(stage_caught)
}

adh_sim(13000, 60, 6)
