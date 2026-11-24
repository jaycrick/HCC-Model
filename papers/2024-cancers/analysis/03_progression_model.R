# papers/2024-cancers/analysis/03_progression_model.R
#
# Progression model figures: stage distribution over time, screening
# frequency effects, and attendance analyses.
# Requires data-derived/resampled_pars_fin.rds (from 01_bootstrap_params.R).
# Bootstrap bands use Sym_B = 8 * Sym_A per draw (prepare_draws(), finding C7).
#
# Figures written to figures/:
#   Untreated_HCC_Distribution_by_Risk_10_Year.jpg
#   Untreated_HCC_over_time_10_Year.jpg
#   HCC_Years_of_Life_Saved_by_Testing_Frequency_10_Year.jpg
#   HCC_Years_of_Life_Saved_by_Testing_Frequency_10_log_Year.jpg
#   HCC_Years_of_Life_Saved_by_p_CA_10_Year.jpg
#   attendance_plot_10_year.{jpg,png,pdf}
#   attendance_plot_stages_10_year.jpg
#   attendance_plot_stages_ribbon_10_year.jpg

# ----- Initialise -------------------------------------------------------------
options(hcc.paper = "2024-cancers") # paths: paper_path() in R/setup.R
source(here::here("R", "setup.R"))
source(here::here("R", "model_core.R"))
source(here::here("R", "plot_theme.R"))
source(here::here("R", "plot_sweeps.R"))
source(here::here("R", "plot_scenarios.R"))
source(here::here("R", "plot_calibration.R"))

years <- 10
months <- years * 12
resampled_pars <- readRDS(paper_path("data-derived", "resampled_pars_fin.rds"))

# ----- Untreated stage distribution -------------------------------------------
# Across a range of annual HCC risk
plot_stage_distribution(
  stage_distribution_by_risk(years),
  x_lab = "Annual HCC Risk Profile",
  title = paste(
    "Model-estimated Distribution of BCLC\nStages From Untreated HCC After",
    years,
    "Years"
  ),
  caption = "With the assumption no-one has\nHCC at the start of the model",
  x_scale = scale_x_continuous(labels = label_percent(), expand = expansion(0))
) |>
  save_figure(
    paste0("Untreated_HCC_Distribution_by_Risk_", years, "_Year"),
    width = 5.5,
    height = 6,
    formats = "jpg",
    dpi = 300
  )

# Over time at the default risk
plot_stage_distribution(
  stage_distribution_over_time(years),
  x_lab = "Years Since Model Start",
  title = paste(
    "Model-estimated Distribution of BCLC\nStages from Untreated HCC Over",
    years,
    "Years"
  ),
  caption = paste(
    "With the assumption no-one has\nHCC at the start of the model and a",
    params_def$p_HCC,
    "annual HCC risk"
  )
) |>
  save_figure(
    paste0("Untreated_HCC_over_time_", years, "_Year"),
    width = 5.5,
    height = 6,
    formats = "jpg",
    dpi = 300
  )

# ----- Screening frequency vs life years --------------------------------------
usfreq <- sweep_life_years(1:24, "m", resampled_pars, t = months)
usfreq_title <- paste(
  "Effect of Ultrasound Screening Frequency\non Life Years over",
  years,
  "Years"
)
base <- paste0("HCC_Years_of_Life_Saved_by_Testing_Frequency_", years)
plot_sweep_band(usfreq, "Ultrasound Frequency", title = usfreq_title) |>
  save_figure(paste0(base, "_Year"), 5.5, 6, "jpg", dpi = 300)
plot_sweep_band(
  usfreq,
  "Ultrasound Frequency (log scale)",
  title = usfreq_title,
  x_scale = scale_x_log10(
    breaks = round(10^(seq(0, log10(24), length.out = 7))),
    expand = expansion(0)
  )
) |>
  save_figure(paste0(base, "_log_Year"), 5.5, 6, "jpg", dpi = 300)

# ----- Complete-adherence proportion vs life years ----------------------------
p_ca <- sweep_life_years(
  seq(0, 1 - params_def$p_NA, length.out = 50),
  "p_CA",
  resampled_pars,
  t = months
)
plot_sweep_band(
  p_ca,
  "Population in Completely Adherent Cohort",
  title = paste("Effect of Complete Adherence\non Life Years over", years, "Years"),
  caption = paste0(
    "People are taken from the partially\nadherent cohort, whilst the ",
    "non-adherent\ncohort is assumed to be stationary"
  ),
  x_scale = scale_x_continuous(labels = label_percent(), expand = expansion(0))
) |>
  save_figure(
    paste0("HCC_Years_of_Life_Saved_by_p_CA_", years, "_Year"),
    5.5,
    9,
    "jpg",
    dpi = 300
  )

# ----- Attendance probability vs life years -----------------------------------
no_adherents <- list(p_NA = 0, p_CA = 0)
attendance <- sweep_life_years(
  (0:20) / 20,
  "p_Att",
  resampled_pars,
  fixed = no_adherents,
  t = months
)
plot_sweep_band(
  attendance,
  "Probability of Attending\nScheduled Ultrasound Screenings",
  x_scale = scale_x_continuous(labels = label_percent(), expand = expansion(0)),
  y_scale = NULL,
  linewidth = 1.5,
  fill_title = "Confidence\nBand"
) |>
  save_figure(paste0("attendance_plot_", years, "_year"), 6, 5.5, dpi = FIGURE_DPI_ALL)

# ----- Stage-specific diagnoses vs attendance ---------------------------------
stages <- sweep_stage_diagnoses(
  (0:20) / 20,
  "p_Att",
  resampled_pars,
  fixed = no_adherents,
  t = months
)
stages_caption <- paste("In 100 people over", years, "years")
plot_stage_diagnoses(
  stages,
  paste0(
    "Effect of Screening Attendance within the\n",
    "Partial Adherence Cohort on Diagnosis of HCC"
  ),
  stages_caption
) |>
  save_figure(paste0("attendance_plot_stages_", years, "_year"), 5, 8, "jpg", dpi = 300)
plot_stage_diagnoses(
  stages,
  paste0(
    "Effect of Screening Attendance\non Stage of Diagnosis within\n",
    "the Partial Adherence Cohort"
  ),
  stages_caption,
  ribbon = TRUE
) |>
  save_figure(
    paste0("attendance_plot_stages_ribbon_", years, "_year"),
    5,
    8,
    "jpg",
    dpi = 300
  )

# Note: stochastic simulation functions (simulate_cohort,
# simulate_by_adherence) in R/simulation.R check the deterministic model and
# are not used in the published main analysis.
