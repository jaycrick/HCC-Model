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
│   ├── progression.R       # Figure-generating analysis functions
│   ├── simulation.R        # Stochastic simulation functions (exploratory)
│   └── bootstrap.R         # Bootstrap likelihood and scenario functions
├── analysis/               # Numbered analysis scripts (run in order)
│   ├── 01_bootstrap_params.R   # Parameter estimation (HPC; skip with shipped data)
│   ├── 02_sensitivity.R        # Univariate sensitivity analysis
│   ├── 03_progression_model.R  # Progression & screening-frequency figures
│   └── 04_clinical_scenarios.R # Clinical scenario comparisons
├── data/                   # Read-only parameter inputs
│   ├── Parameters.csv          # Point estimates and uncertainty bounds
│   └── Parameters_25.csv       # ±25% variation bounds
├── data-derived/           # Pre-computed bootstrap outputs (Git LFS)
│   ├── miltonoutput.rds             # Raw MLE outputs from HPC cluster
│   ├── list_of_resampled_milton_output.rds  # Bootstrapped input data
│   ├── resampled_pars.rds           # Bootstrapped transition parameters (n = 10 000)
│   ├── resampled_pars_fin.rds       # As above + scenario life-year columns
│   └── param_names.rds             # Named list of reported manuscript values
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
| `01_bootstrap_params.R` | `miltonoutput.rds`, `list_of_resampled_*.rds` | `resampled_pars_fin.rds`, `param_names.rds` |
| `02_sensitivity.R` | `Parameters.csv`, `Parameters_25.csv`, `param_names.rds` | `param_names.rds`, figures |
| `03_progression_model.R` | `resampled_pars_fin.rds` | figures |
| `04_clinical_scenarios.R` | `resampled_pars_fin.rds`, `param_names.rds` | `param_names.rds`, figures |

`run_all.R` enforces this order. It skips `01_bootstrap_params.R` by default because the
bootstrap outputs are already shipped in `data-derived/`. Pass `--bootstrap` to
regenerate them (requires substantial compute; originally run on an HPC cluster):

```bash
Rscript run_all.R --bootstrap
```

---

## Interactive use

To explore or modify the model outside the pipeline, source the two core files and call functions directly:

```r
source("R/setup.R")      # loads packages, reads Parameters.csv, sets *_def globals
source("R/model_core.R") # defines fromf4(), us_theor_pa(), adh_theor(), life_add()

# Life-years gained under default parameters (current practice)
life_add()
#> [1] 68.02326

# No surveillance counterfactual
life_add(p_NA = 1, p_CA = 0)
#> [1] 14.85677

# What if ultrasound interval shortened to 3 months?
life_add(m = 3)

# What if complete adherence increases from default to 50%?
life_add(p_CA = 0.5)
```

`life_add()` accepts any combination of the model's ~25 parameters as named arguments;
all default to the `*_def` point estimates read from `data/Parameters.csv`.

Key functions:

| Function | Purpose |
|----------|---------|
| `fromf4(n, t, ...)` | Markov progression for `n` people over `t` months |
| `us_theor_pa(n, t, m, Pr_Att, ...)` | Ultrasound detection for a single adherence group |
| `adh_theor(n, t, m, p_NA, p_CA, p_Att, ...)` | Detection across all three adherence groups |
| `life_add(...)` | Total life-years gained; wraps `adh_theor` |
| `life_tibble(csv)` | Builds parameter sweep tibble for sensitivity tornado plots |

---

## Data provenance

| File | Source | Description |
|------|--------|-------------|
| `Parameters.csv` | Literature / expert elicitation | Point estimates and 95% bounds for all model parameters |
| `Parameters_25.csv` | Derived from Parameters.csv | ±25% variation bounds for modifiable parameters |
| `miltonoutput.rds` | HPC cluster MLE run | 10 000 MLE optimisation results from bootstrapped input data |
| `list_of_resampled_milton_output.rds` | Parametric bootstrap | Bootstrapped summary statistics used as MLE inputs |
| `resampled_pars.rds` | Derived from `miltonoutput.rds` | Extracted bootstrapped transition + sensitivity parameters |
| `resampled_pars_fin.rds` | Derived from `resampled_pars.rds` | As above with scenario life-year columns appended |
| `param_names.rds` | Computed during pipeline | Named list of all values reported in the manuscript text |

---

## Key model outputs (regression anchors)

These values should be reproduced exactly from a clean run:

- `life_add()` (baseline, default parameters) ≈ **68.02326** life years per cohort of 100
- `life_add(p_NA = 1, p_CA = 0)` (no surveillance) ≈ **14.85677** life years

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
| `pre-commit` | Runs `air format` on staged `.R` files and re-stages them automatically, then runs the lint ratchet on them and blocks the commit if any file has more lints than its baseline |
| `commit-msg` | Enforces [git3moji](https://robinpokorny.github.io/git3moji/) prefixes: `:zap:` `:bug:` `:tv:` `:cop:` `:abc:` |

Format any file manually with `air format <file.R>`, or format the whole project with `air format .`.

### Lint ratchet

Lint rules live in `.lintr` (line length 80 to match `air.toml`; `object_name_linter` and `indentation_linter` are disabled because the code follows paper notation such as `p_HCC` and air owns formatting).
The legacy code still has many lints, so instead of a hard zero-lint gate the repository uses a ratchet: `tools/lint_baseline.csv` records each file's allowed lint count, and a file may never exceed it.
A file not in the baseline must have zero lints.

```bash
Rscript tools/lint_ratchet.R [files...]   # check (default: all project .R files)
Rscript tools/lint_ratchet.R --update     # lower the baseline after fixing lints
```

`--update` refuses to raise any count, so the baseline can only go down.
The target is an empty baseline (all counts 0).

### Tests

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

The suite takes a few seconds.
Never run `run_all.R` as a test; it takes about 30 minutes.

### Continuous integration

`.github/workflows/check.yml` runs on every push and pull request (Ubuntu, R version and packages from `renv.lock`).
It runs `air format --check .`, the lint ratchet, and the test suite.

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
