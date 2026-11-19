# analysis/03_progression_model.R
#
# Progression model figures: stage distribution over time, screening
# frequency effects, and attendance analyses.
# Requires data-derived/resampled_pars_fin.rds (from 01_bootstrap_params.R).
#
# Figures written to figures/:
#   Untreated_HCp_Cistribution_by_Risk_10_Year.jpg
#   Untreated_HCC_over_time_10_Year.jpg
#   HCC_Years_of_Life_Saved_by_Testing_Frequency_10_Year.jpg
#   HCC_Years_of_Life_Saved_by_Testing_Frequency_10_log_Year.jpg
#   HCC_Years_of_Life_Saved_by_p_CA_10_Year.jpg
#   attendance_plot_10_year.{jpg,png,pdf}
#   attendance_plot_stages_10_year.jpg
#   attendance_plot_stages_ribbon_10_year.jpg

# ----- Initialise -------------------------------------------------------------
source(here::here("R", "setup.R"))
source(here::here("R", "model_core.R"))
source(here::here("R", "progression.R"))

library(RColorBrewer)
library(ggthemes)
library(gridExtra)
library(magrittr)

dir.create(here("figures"), showWarnings = FALSE)

# ----- Bootstrap parameter draws ----------------------------------------------
resampled_pars <- readRDS(here("data-derived", "resampled_pars_fin.rds"))

# ----- Stage distribution across a range of HCC risk ------------------------
plot_untreated_risk_10_years <- risk_no_int(100, 10)

# ----- Stage distribution over time ------------------------------------------
plot_untreated_prog_10_years <- years_no_int(10)

# ----- Screening frequency vs life years -------------------------------------
plot_usfreq_10_years <- usfreq_graph(10)

# ----- Complete-adherence proportion vs life years ---------------------------
plot_partial_adherence_10_years <- pafunc(10)

# ----- Attendance probability vs life years (inline analysis) ----------------
p_Att <- (0:20) / 20
tempfun <- function(p)
  pmap_dbl(draws(resampled_pars), p_NA = 0, p_CA = 0, p_Att = p, life_add)
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
  scale_x_continuous(labels = label_percent(), expand = expansion(0)) +
  geom_ribbon(aes(ymin = min, ymax = max, fill = "95%"), alpha = 0.2) +
  scale_fill_manual(values = "black") +
  labs(
    fill = "Confidence\nBand",
    x = "Probability of Attending\nScheduled Ultrasound Screenings",
    y = "Additional Life Years Gained From\nIntervention per Cohort of 100"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  )

plot_att

ggsave(
  here("figures", "attendance_plot_10_year.jpg"),
  plot_att,
  width = 6,
  height = 5.5
)
ggsave(
  here("figures", "attendance_plot_10_year.png"),
  plot_att,
  width = 6,
  height = 5.5,
  dpi = 600
)
ggsave(
  here("figures", "attendance_plot_10_year.pdf"),
  plot_att,
  width = 6,
  height = 5.5,
  dpi = 600
)

# ----- Stage-specific diagnoses vs attendance --------------------------------
plot_att_stages_10_years <- p_Att_Stages(10)

# Note: stochastic simulation functions (hccsim, us_sim, adh_sim) are
# available in R/simulation.R for exploratory use but are not used in the
# published main analysis.
