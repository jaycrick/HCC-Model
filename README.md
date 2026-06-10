# HCC Surveillance Markov Model

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

```bash
git clone <repo-url>
cd HCC-Model
git lfs pull            # download data-derived/ files
Rscript -e 'renv::restore()'
Rscript run_all.R
```

All figures are written to `figures/`.

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
