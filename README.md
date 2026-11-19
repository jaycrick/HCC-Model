# HCC Surveillance Markov Model

[![git3moji](https://img.shields.io/badge/git3moji-%E2%9A%A1%EF%B8%8F%F0%9F%90%9B%F0%9F%93%BA%F0%9F%91%AE%F0%9F%94%A4-fffad8.svg?style=flat-square)](https://robinpokorny.github.io/git3moji/)

Code for the model described in:

> Cumming J, Scott N, Howell J, Flores JE, Pavlyshyn D, Hellard ME, Winata LS, Ryan M, Sutherland T, Thompson AJ, Doyle JS, Sacks-Davis R.
> **Improving Hepatocellular Carcinoma Surveillance Outcomes in Patients with Cirrhosis after Hepatitis C Cure: A Modelling Study.**
> *Cancers* 2024, **16**(15), 2745. https://doi.org/10.3390/cancers16152745

This repository contains an R implementation of a deterministic compartmental Markov
model for hepatocellular carcinoma (HCC) surveillance in patients with cirrhosis after
hepatitis C cure. The model tracks cohorts through BCLC disease stages (0, A, B, C, D)
and death, incorporating ultrasound screening sensitivity, patient adherence (non-,
partial, and completely adherent groups), and stage-specific treatment outcomes to
estimate life-years gained under various surveillance-improvement scenarios relative to
a current-practice baseline over a 10-year horizon.


> [!IMPORTANT]
> **Results differ from the publication.** In November 2026 the calibration was
> re-done with the exact likelihood and a better optimiser, and all 10,000 bootstrap
> resamples were refitted (gates G8, G9). The scenario point estimates still match
> the paper (9.8, 17.2, 7.0, 9.0 life-years per 100), but several confidence intervals
> and the partial correlations do not. See
> [Errata and differences from the publication](#errata-and-differences-from-the-publication)
> at the bottom of this file.

---

## Quick start

```bash
git clone https://github.com/jaycrick/HCC-Model.git && cd HCC-Model
git lfs pull                          # download pre-computed bootstrap data
Rscript -e 'renv::restore()'         # install packages at locked versions
Rscript run_all.R                     # reproduce all figures → figures/
```

Open `MODEL.Rproj` in RStudio or Positron to work interactively.

---

## Requirements

- **R** ≥ 4.2 (version used: see `renv.lock`)
- **renv** for reproducible package management
- **Git LFS** to download the pre-computed bootstrap data files

### Install Git LFS (once per machine)

```bash
git lfs install
```

### Restore the R package library

```r
renv::restore()
```

This installs all packages at the exact versions recorded in `renv.lock`.

> **Note on `hypogeom`:** This package was written specifically for this project by
> co-author Damian Pavlyshyn and is only available from GitHub
> (`damian-t-p/hypogeom`). `renv::restore()` will install it automatically. If you
> need to install it manually: `remotes::install_github("damian-t-p/hypogeom")`.

---

## Repository layout

```
.
├── run_all.R               # Master script — reproduces all figures
├── R/                      # Pure function definitions (sourced, not run directly)
│   ├── setup.R             # Libraries, seeds, constants, default parameters
│   ├── model_core.R        # Markov transitions, US detection, adherence model
│   ├── plot_theme.R        # theme_hcc(), save_figure(), interval helpers
│   ├── plot_sweeps.R       # Bootstrap sweeps, band plots, stage distributions
│   ├── plot_scenarios.R    # Tornado, scenario and response-curve plots
│   ├── plot_calibration.R  # Bootstrap calibration figures (analysis/01)
│   ├── simulation.R        # Stochastic check of the cohort model
│   ├── scenarios.R         # Scenario definitions as data; run_scenarios()
│   └── bootstrap.R         # Calibration data, likelihood, bootstrap inputs
├── analysis/               # Numbered analysis scripts (run in order)
│   ├── 01_bootstrap_params.R   # Scenarios + figures from the refitted bootstrap draws (skip by default); reads bootstrap_fit_diagnostics.rds
│   ├── 02_sensitivity.R        # Univariate sensitivity analysis
│   ├── 03_progression_model.R  # Progression & screening-frequency figures
│   └── 04_clinical_scenarios.R # Clinical scenario comparisons
├── data/                   # Read-only parameter inputs
│   ├── Parameters.csv          # Point estimates and uncertainty bounds
│   ├── Parameters_sources.csv  # Source and table reference per parameter
│   └── Parameters_25.csv       # ±25% variation bounds
├── data-derived/           # Pre-computed bootstrap outputs (Git LFS)
│   ├── miltonoutput.rds             # Legacy HPC MLEs, superseded by the local refit (G9); kept for provenance
│   ├── list_of_resampled_milton_output.rds  # Bootstrapped input data
│   ├── resampled_pars.rds           # Refitted bootstrap parameters (n = 10 000; tools/refit_calibration.R)
│   ├── bootstrap_fit_diagnostics.rds # Refit diagnostics (convergence, fallback, at-bound, audit)
│   ├── resampled_pars_fin.rds       # As above + scenario life-year columns
│   └── param_names.rds             # Named list of reported manuscript values
├── output/                 # manuscript_values.csv: every number quoted in the paper
├── tests/testthat/         # Test suite (see Tests)
├── tools/                  # lint_ratchet.R, check_globals.R
└── figures/                # Generated figures (git-ignored; reproduced by pipeline)
```

---

## How to reproduce the figures

Run `Rscript run_all.R` from the project root (see [Quick start](#quick-start)). All figures are written to `figures/`.

### Run order

Scripts must run in numerical order because downstream scripts read from
`data-derived/`:

| Script | Reads | Writes |
|--------|-------|--------|
| `01_bootstrap_params.R` | `resampled_pars.rds`, `bootstrap_fit_diagnostics.rds`, `list_of_resampled_*.rds` | `resampled_pars_fin.rds`, figures |
| `02_sensitivity.R` | `Parameters.csv`, `Parameters_25.csv`, `param_names.rds` | `param_names.rds`, figures |
| `03_progression_model.R` | `resampled_pars_fin.rds` | figures |
| `04_clinical_scenarios.R` | `resampled_pars_fin.rds`, `param_names.rds` | `param_names.rds`, figures |

`run_all.R` enforces this order. It skips `01_bootstrap_params.R` by default because the
bootstrap outputs are already shipped in `data-derived/`. Pass `--bootstrap` to
regenerate them (requires substantial compute; originally run on an HPC cluster):

```bash
Rscript run_all.R --bootstrap
```

The expensive calibration refit (point MLE plus 10 000 bootstrap resamples) is a
separate one-off step, `tools/refit_calibration.R`; it writes
`data-derived/resampled_pars.rds`, `bootstrap_fit_diagnostics.rds` and the
calibrated rows of `data/Parameters.csv` and `Parameters_25.csv`:

```bash
Rscript tools/refit_calibration.R [--n 20] [--cores 10] [--config a|b|c] [--dry-run] [--out DIR]
```

`--dry-run` writes only to `--out` (default a temporary directory). Expected
runtime: config d (default): about 1.5-3 h on a 10-core laptop (estimate); configs a-c are 3-10x slower.

---

## Interactive use

To explore or modify the model outside the pipeline, source the two core files and call functions directly:

```r
source("R/setup.R")      # loads packages; reads Parameters.csv into `params_def`
source("R/model_core.R") # defines progress_untreated(), screen_cohort(), screen_by_adherence(), life_add()
# (R/bootstrap.R and R/scenarios.R are sourced the same way when needed)

# Life-years gained under default parameters (current practice)
life_add()
#> [1] 68.03421

# No surveillance counterfactual
life_add(p_NA = 1, p_CA = 0)
#> [1] 14.86345

# What if ultrasound interval shortened to 3 months?
life_add(m = 3)

# What if complete adherence increases from default to 50%?
life_add(p_CA = 0.5)
```

`life_add()` accepts any combination of the model's 24 parameters as named arguments;
all default to the point estimates in `params_def` (read from `data/Parameters.csv`).
Unknown names are an error (e.g. `life_add(p_att = 1)`), and `n`, `t`, `m` and `params` are
ordinary named arguments. To run on a different parameter set, pass
`params = modify_params(params_def, list(p_HCC = 0.03))` or `load_params(csv)`.

Key functions:

Every function has a roxygen-style header (purpose, `@param`, `@return`) with
a `Ref:` line that points to the section, equation, table or figure of Cumming
et al. 2024 or its Supplement (`Ref: none` where the code has no counterpart in
the publication). Variables follow the paper's notation (`p_HCC`, `p_0`..`p_D`,
`Sens_*`, `Sym_*`, `L_*`, `p_NA`, `p_CA`, `p_Att`).

| Function | File | Purpose |
|----------|------|---------|
| `load_params(csv, column)` | setup.R | Read a Parameters CSV into a named list (by `PARAMETER` name) |
| `modify_params(params, overrides)` | model_core.R | Override named parameters; errors on unknown names |
| `annual_to_monthly(p)` | model_core.R | Convert an annual probability to monthly |
| `transition_matrix(p_HCC_month, p_stage)` | model_core.R | The 7x7 monthly transition matrix |
| `progress_untreated(..., n, t, params)` | model_core.R | Markov progression of an untreated cohort over `t` months |
| `screen_cohort(Pr_Att, ..., n, t, m, params)` | model_core.R | Ultrasound detection for one adherence group; returns `list(dx, undiagnosed)` |
| `screen_by_adherence(..., n, t, m, params)` | model_core.R | Detection across all three adherence groups |
| `life_add(..., n, t, m, params)` | model_core.R | Total life-years gained; wraps `screen_by_adherence` |
| `life_tibble(csv)` | model_core.R | Parameter sweep tibble for the sensitivity tornado plots |
| `prepare_draws(df)`, `draws(df)`, `years(df)` | model_core.R | Bootstrap draws ready for the model (`Sym_B = 8 * Sym_A`); life-years per draw |
| `scenario_defs`, `risk_multipliers`, `run_scenarios(param_rows)` | scenarios.R | Scenarios as data; life-years per scenario x risk level x draw |
| `scenario_columns(long)` | scenarios.R | Pivot `run_scenarios()` output to the legacy wide columns |
| `record_values(values, path)` | setup.R | Merge named values into `output/manuscript_values.csv` |
| `calibration_data` | bootstrap.R | Observed calibration data (Supp. Sec. 2.1.1) |
| `likelihood_matrix(p_0, p_A, sens_0, sens_A)` | bootstrap.R | The 15x15 transition matrix of Supp. Table S2 / matrix S1 |
| `end_state_probs(tm, sens_0, t)` | bootstrap.R | Probability of ending in each state after `t` months |
| `neg_log_lik(par, ...)` | bootstrap.R | Negative log-likelihood of the 8 calibrated parameters |
| `fit_point_estimate(start)` | bootstrap.R | Maximum-likelihood point estimate (thin wrapper of `fit_mle()`: Nelder-Mead then bounded L-BFGS-B polish) |
| `bootstrap_inputs_long()`, `bootstrap_estimates()` | bootstrap.R | Read the shipped bootstrap inputs and MLEs |
| `simulate_cohort(Pr_Att, ...)`, `simulate_by_adherence(...)` | simulation.R | Stochastic (binomial) check of the cohort model; not in the publication |
| `theme_hcc()`, `save_figure()`, `sweep_*()`, `plot_*()` | plot_*.R | Figure helpers: data builders, plot functions and saving |

---

## Data provenance

| File | Source | Description |
|------|--------|-------------|
| `Parameters.csv` | Cumming 2024, Tables 1-3 | Point estimates and 95% bounds for all model parameters (point estimate = published MLE) |
| `Parameters_sources.csv` | Cumming 2024, Tables 2-3 | `PARAMETER`, `SOURCE`, `TABLE_REF` for each row of `Parameters.csv`; unknowns are marked `TODO: confirm with author` |
| `Parameters_25.csv` | Derived from Parameters.csv | ±25% variation bounds for modifiable parameters |
| `miltonoutput.rds` | HPC cluster MLE run | Legacy HPC MLEs, superseded by the local refit (G9); kept for provenance |
| `list_of_resampled_milton_output.rds` | Parametric bootstrap | Bootstrapped summary statistics used as MLE inputs |
| `resampled_pars.rds` | `tools/refit_calibration.R` (G9) | Refitted bootstrap MLEs (8 transition + sensitivity parameters, 10 000 rows) |
| `bootstrap_fit_diagnostics.rds` | `tools/refit_calibration.R` (G9) | Per-resample NLL, convergence code, best start, at-bound flags, old HPC NLL |
| `resampled_pars_fin.rds` | Derived from `resampled_pars.rds` (by `01_bootstrap_params.R`) | As above with scenario life-year columns appended |
| `param_names.rds` | Computed during pipeline | Named list of all values reported in the manuscript text |

---

## Key model outputs (regression anchors)

These values should be reproduced exactly from a clean run:

- `life_add()` (baseline, default parameters) ≈ **68.03421** life years per cohort of 100
- `life_add(p_NA = 1, p_CA = 0)` (no surveillance) ≈ **14.86345** life years

Both are asserted via `stopifnot` in `04_clinical_scenarios.R`.

---

## Contributing

### First-time setup

After cloning, activate the versioned git hooks and install the [air](https://github.com/posit-dev/air) R formatter:

```bash
git config core.hooksPath .githooks

# Install air (if not already present)
curl -fsSL https://github.com/posit-dev/air/releases/latest/download/air-installer.sh | sh
```

### What the hooks enforce

| Hook | What it does |
|------|-------------|
| `pre-commit` | Runs `air format` on staged `.R` files and re-stages them automatically, then runs the lint ratchet on them and blocks the commit if any file has more lints than its baseline; finally `tools/check_globals.R` (codetools) blocks the commit if any project function uses an undefined name |
| `commit-msg` | Enforces [git3moji](https://robinpokorny.github.io/git3moji/) prefixes: `:zap:` `:bug:` `:tv:` `:cop:` `:abc:` |

Format any file manually with `air format <file.R>`, or format the whole project with `air format .`.

### Lint and format

Formatting is owned by [air](https://github.com/posit-dev/air) (`air.toml`, line width 80): `air format .` formats, `air format --check .` checks.
Enable the versioned hooks once per clone with `git config core.hooksPath .githooks`.
`tools/check_globals.R` uses codetools to find any project function that uses an undefined name (a typo or a missing global); the hook and CI both run it.

Lint rules live in `.lintr` (line length 80 to match `air.toml`; `object_usage_linter` is off because it cannot see across the sourced `R/` files, and `tools/check_globals.R` does that job; `object_name_linter` and `indentation_linter` are disabled because the code follows paper notation such as `p_HCC` and air owns formatting).
The legacy code still has many lints, so instead of a hard zero-lint gate the repository uses a ratchet: `tools/lint_baseline.csv` records each file's allowed lint count, and a file may never exceed it.
A file not in the baseline must have zero lints.

```bash
Rscript tools/lint_ratchet.R [files...]   # check (default: all project .R files)
Rscript tools/lint_ratchet.R --update     # lower the baseline after fixing lints
Rscript tools/check_globals.R             # undefined names
```

`--update` refuses to raise any count, so the baseline can only go down.
The target is an empty baseline (all counts 0).

### Tests

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

The suite takes a few seconds.
Never run `run_all.R` as a test; it takes about 1-2 minutes (without `--bootstrap`).
The tests fall into four groups:

- **Golden** (`test-golden.R`, `tests/testthat/fixtures/golden.rds`): pins the code as it is, bugs included. The fixture is generated once by `fixtures/make_golden.R` and regenerated only under an approved, recorded change (see `CHANGES.md`). Doubles are compared at relative tolerance **1e-12**: portable across platforms (libm/BLAS noise is about 1e-15), while real model changes are at least 1e-3.
- **Invariants** (`test-invariants.R`, `test-core.R`, `test-simulation.R`, `test-likelihood.R`): properties that must hold whatever the numbers are, e.g. rows of the transition matrix sum to 1, people are conserved (`sum(dx) + sum(undiagnosed) == n`), the stochastic model agrees with the deterministic one in expectation, and the likelihood matches its closed form.
- **Paper anchors** (`test-paper.R`, `test-scenarios.R`): published values (abstract, Sec. 3, Table 3) after rounding to the printed precision.
- **Data and ledger** (`test-data.R`, `test-ledger.R`, `test-plotting.R`): CSV structure (including that `Parameters_sources.csv` lists the same parameters in the same order), ledger format, plot helpers.

### Continuous integration

`.github/workflows/check.yml` runs on every push and pull request (Ubuntu, R version and packages from `renv.lock`).
It runs `air format --check .` (air pinned to the version in the workflow), the lint ratchet, `tools/check_globals.R` and the test suite.

---

## Outputs

- `output/manuscript_values.csv` (columns `name`, `value`) is the ledger of every number quoted in the manuscript, written by `record_values()`. It is sorted by name and written at a precision that round-trips exactly, so it is text-diffable.
- `data-derived/` is read-only input (the shipped bootstrap data). Nothing writes to it except the deliberate `Rscript run_all.R --bootstrap` regeneration.
- `figures/` is generated by the pipeline and git-ignored.

---

## Citation

Please cite the associated manuscript:

```bibtex
@Article{cancers16152745,
  author   = {Cumming, Jacob and Scott, Nick and Howell, Jessica and
              Flores, Joan Ericka and Pavlyshyn, Damian and
              Hellard, Margaret E. and Winata, Leon Shin-han and
              Ryan, Marno and Sutherland, Tom and Thompson, Alexander J. and
              Doyle, Joseph S. and Sacks-Davis, Rachel},
  title    = {Improving Hepatocellular Carcinoma Surveillance Outcomes in
              Patients with Cirrhosis after Hepatitis C Cure: A Modelling Study},
  journal  = {Cancers},
  volume   = {16},
  year     = {2024},
  number   = {15},
  article-number = {2745},
  url      = {https://www.mdpi.com/2072-6694/16/15/2745},
  doi      = {10.3390/cancers16152745},
  pmid     = {39123472},
  issn     = {2072-6694},
}
```

See also `CITATION.cff` for a machine-readable citation.

---

## License

MIT — see `LICENSE`.

---

## Errata and differences from the publication

Cumming et al. 2024, *Cancers* 16(15):2745, and its Supplement. Full details,
before/after tables and the complete ledger diff are in `CHANGES.md` (gates G8, G9;
findings P1-P7). The published values are kept, frozen, in
`data-derived/param_names.rds` and in the comments of `tests/testthat/test-paper.R`.

### A. Published numbers the current code no longer reproduces (G8, G9)

Two methodological fixes change the calibrated parameters:

- **G8 (finding P5):** the likelihood uses the exact absorption probabilities of the
  absorbing Markov chain instead of the T^100 approximation.
- **G9 (finding P7):** the maximum-likelihood fit is Nelder-Mead followed by a
  bounded L-BFGS-B on the logit scale, and all 10,000 bootstrap resamples were
  refitted locally (`tools/refit_calibration.R`). Every refit is at least as good as
  the original HPC fit. The original fits had stopped early in the tails of the
  bootstrap distribution: in 418 resamples (4.2%) the best fit has ultrasound
  sensitivity in BCLC 0 (Sens_0) at its upper limit of 1.

As a result, the point estimates barely move but the bootstrap intervals widen:

| Quantity | Published | Current code | Location in paper |
|---|---|---|---|
| Sens_0, 95% CI | 0.08 to 0.54 | 0.08 to 1.00 | Table 3 |
| p_0, 95% CI | 0.031 to 0.240 | 0.030 to 0.482 | Table 3 |
| p_B, 95% CI | 0.20 to 0.74 | 0.21 to 0.89 | Table 3 |
| p_C, 95% CI | 0.22 to 0.57 | 0.21 to 0.51 | Table 3 |
| Realistic adherence, 95% CI | 7.9 to 11.6 | 7.9 to 11.7 | Abstract, Sec. 3.3 |
| Optimistic adherence, 95% CI | 13.9 to 20.3 | 13.8 to 20.5 | Abstract, Sec. 3.3 |
| Diagnosis sensitivity, 95% CI | 3.6 to 13.8 | **2.5** to 14.1 | Abstract, Sec. 3.3 |
| Treatment, 95% CI | 7.5 to 10.3 | 7.4 to 10.3 | Abstract, Sec. 3.3 |
| Optimistic > diagnosis (%) | 98.4 | 98.1 | Table 4 |
| Realistic > diagnosis (%) | 80.1 | 77.4 | Table 4 |
| Treatment > diagnosis (%) | 74.7 | 71.4 | Table 4 |
| Realistic > treatment (%) | 73.6 | 73.9 | Table 4 |
| Mean months in stage, lower CI (0, B, C, D) | 4.2, 1.4, 1.7, 6.1 | 2.1, 1.1, 2.0, 6.2 | Sec. 3.1, Fig. 2A |
| Mean months in stage, upper CI (0, B, C, D) | 32.6, 4.9, 4.5, 7.8 | 33.0, 4.8, 4.8, 7.7 | Sec. 3.1, Fig. 2A |
| Bootstraps with stage A longer than stage 0 (%) | 73 | 70 | Sec. 3.1 |
| SE of baseline life-years | 5.7 | 6.0 | Sec. 3.4 |
| Partial correlation, life-years ~ logit(p_0) | -0.87 | +0.29 | Sec. 3.4 |
| Partial correlation, life-years ~ logit(Sens_0) | 0.92 | -0.02 | Sec. 3.4 |
| +/-25% sensitivity, L_A upper (years) | 5.7 | 5.6 | Sec. 3.2 |

Unchanged at the paper's printed precision: the scenario point estimates
(9.8, 17.2, 7.0 and 9.0 life-years per cohort of 100), the other Table 3 intervals and
the other +/-25% results. `life_add()` moves from 68.023 to 68.034.

**Partial correlations.** The reversal is an artefact of the logit scale: the boundary
fits (Sens_0 near 1) become extreme outliers on that scale and dominate the partial
correlations. Rank-based (Spearman) partial correlations are robust to this and agree
with the published conclusion: life-years ~ p_0 **-0.85**, ~ Sens_0 **+0.90**,
~ Sens_A **+0.51** (ledger names `pcor_rank_0`, `pcor_rank_US_0`, `pcor_rank_US_A`).
The logit analysis is kept as the published method.

### B. Errors and inconsistencies in the publication text

The items below are differences between the paper's text and the analysis that
produced its numbers. Except where marked closed, the code keeps the analysis
behaviour. Each item is also a comment at the relevant line of code.

| ID | Where | What |
|----|-------|------|
| P1 | Supplement p. 6 | The printed Pr(d_alpha) omits the (1 - Sym) factors; the code uses the correct form (0.3289 at the MLE against 52/158 = 0.329 observed). |
| P2 | Paper Sec. 2.5 | The text gives low / high annual HCC risks of 1.1% / 2.2% (2.2% is also the baseline); the code, and Supplement Sec. 5 ("half and double"), use x0.5 and x2 of p_HCC (1.1% and 4.4%). The Supplement Figs. S6-S7 captions say 1% and 4%. |
| P3 | Paper Sec. 2.2.2 | "52 of n = 156"; the Supplement (s-bar = 158) and the code use 158. |
| P4 | Supplement Fig. S3 caption | Says N(mu = 68.0); 68.0 is the point estimate `life_add()`. The fitted normal in the code uses the bootstrap mean (about 67.6 published; 67.5 after G9). |
| P5 | Supplement Table S2 | Lists 14 states; the code adds a 15th absorbing state `past_A` (progressed past BCLC A undiagnosed). Also, Pr(d_0) uses the t^(100) approximation (error about 1.5e-6) rather than the exact closed form, because it generated the published MLE. **Closed by G8**: the code now uses exact absorption probabilities. See `CHANGES.md`. |
| P6 | Paper Sec. 2.4, scenario 3 | The text lists L_A +1 year and L_B +6 months only; the code, and the published 9.0 life-years, also add L_C + 0.5. |
| P7 | Table 3 / `Parameters.csv` | The published MLE is Nelder-Mead's stopping point at its default tolerance, not an exact optimum (NLL 2e-4 above a tight BFGS optimum); negligible next to the bootstrap intervals. **Closed by G9**: `fit_mle()` polishes with a bounded L-BFGS-B, and the whole bootstrap was refitted locally. See `CHANGES.md`. |
