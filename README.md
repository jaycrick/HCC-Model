# HCC Surveillance Markov Model

Code for the HCC surveillance model described in: **[TODO: add full citation and DOI]**

This repository contains an R implementation of a Markov model for hepatocellular
carcinoma (HCC) surveillance in cirrhotic patients. The model tracks cohorts through
BCLC disease stages (0, A, B, C, D) and death, incorporating ultrasound screening
sensitivity, patient adherence (non-, partial, and completely adherent groups), and
stage-specific treatment outcomes. It estimates life-years gained under various
surveillance-improvement scenarios relative to a current-practice baseline.

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

This installs all packages at the exact versions recorded in `renv.lock`, including
`hypogeom` which is sourced from GitHub (`damian-t-p/hypogeom`).

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

- `life_add()` (baseline, default parameters) ≈ **66.4987** life years per cohort of 100
- `life_add(p_NA = 1, p_CA = 0)` (no surveillance) ≈ **14.19459** life years

Both are asserted via `stopifnot` in `04_clinical_scenarios.R`.

---

## Citation

Please cite the associated manuscript (DOI: TODO) and this software repository:

```
@software{hcc_model,
  title  = {HCC Surveillance Markov Model},
  author = {[TODO]},
  year   = {2024},
  url    = {https://github.com/TODO/HCC-Model},
  doi    = {TODO}
}
```

See also `CITATION.cff` for a machine-readable citation.

---

## License

MIT — see `LICENSE`.
