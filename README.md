# HCC Surveillance Markov Model

[![git3moji](https://img.shields.io/badge/git3moji-%E2%9A%A1%EF%B8%8F%F0%9F%90%9B%F0%9F%93%BA%F0%9F%91%AE%F0%9F%94%A4-fffad8.svg?style=flat-square)](https://robinpokorny.github.io/git3moji/)

An R implementation of a compartmental Markov model for hepatocellular
carcinoma (HCC) surveillance in patients with cirrhosis after hepatitis C cure.
The model tracks cohorts through BCLC disease stages (0, A, B, C, D) and death,
and includes ultrasound screening sensitivity, patient adherence (non-, partial
and completely adherent groups) and stage-specific treatment outcomes. It
estimates life-years gained under surveillance-improvement scenarios relative to
current practice.

The model code in [`R/`](R) is shared by every paper that uses it. Each paper
has its own folder under [`papers/`](papers) with its inputs, analysis scripts
and outputs.

## Papers

| Folder | Paper | Status |
|--------|-------|--------|
| [`papers/2024-cancers`](papers/2024-cancers) | Cumming et al. **Improving Hepatocellular Carcinoma Surveillance Outcomes in Patients with Cirrhosis after Hepatitis C Cure: A Modelling Study.** *Cancers* 2024, 16(15), 2745. https://doi.org/10.3390/cancers16152745 | Published |
| [`papers/stochastic`](papers/stochastic) | Stochastic model with additional compartments (working title) | In preparation |

## Versions and archives

Each released version is a git tag and a
[GitHub release](https://github.com/jaycrick/HCC-Model/releases). Tags never
move. To get the code exactly as it was for a paper, use its tag, not `main`:
`main` keeps changing as later papers extend the model.

To cite the code, cite the paper it was used for (its DOI below). When you
refer to a specific version of the code, also give the tag, because the
paper's DOI identifies the paper, not a version of the code.

| Paper | Paper DOI | Tag | Results match publication? |
|-------|-----------|-----|----------------------------|
| Cumming et al. 2024, *Cancers* | [10.3390/cancers16152745](https://doi.org/10.3390/cancers16152745) | [`v1.0.0-published`](https://github.com/jaycrick/HCC-Model/tree/v1.0.0-published) | **Yes**: the original scripts and data, as used for the paper |
| Cumming et al. 2024, *Cancers* | [10.3390/cancers16152745](https://doi.org/10.3390/cancers16152745) | [`v1.1.0`](https://github.com/jaycrick/HCC-Model/tree/v1.1.0) | **Yes**: refactored and tested code; reproduces the published numbers |
| Cumming et al. 2024, *Cancers* | [10.3390/cancers16152745](https://doi.org/10.3390/cancers16152745) | [`v1.2.0`](https://github.com/jaycrick/HCC-Model/tree/v1.2.0) | **No**: corrected calibration (G8, G9); point estimates match, some intervals do not ([errata](papers/2024-cancers/README.md#errata-and-differences-from-the-publication)) |
| Cumming et al. 2024, *Cancers* | [10.3390/cancers16152745](https://doi.org/10.3390/cancers16152745) | [`v1.3.0`](https://github.com/jaycrick/HCC-Model/tree/v1.3.0) | **No**: same results as `v1.2.0`; per-paper layout, irreplaceable data in plain git (in the download) |
| *(development)* | n/a | `main` | Same results as `v1.3.0` |

A release's source download (zip or tar.gz) does **not** contain the `.rds`
files that can be regenerated from the code: in their place it has Git LFS
pointer files. [Regenerating the .rds files](papers/2024-cancers/README.md#regenerating-the-rds-files)
gives the exact commands for each version. Files that cannot be regenerated are
plain git files, so the downloads contain them. A `git clone` followed by
`git lfs pull` gets everything.

Tags up to and including `v1.3.0` use renv (`renv::restore()`); `main` uses
rig (`rig proj sync`).

---

## Quick start

Assumes rig and Git LFS are already installed (see [Requirements](#requirements)).

```bash
git clone https://github.com/jaycrick/HCC-Model.git && cd HCC-Model
git lfs pull                                   # download the derived .rds files
rig proj sync                                  # install R 4.6 and the locked packages into .rvenv/
rig run papers/2024-cancers/run_all.R          # reproduce paper 1's figures
```

`rig run`, and any R started in the repository root (including RStudio or
Positron via `MODEL.Rproj`), uses the project library through `.Renviron`.
Restart R after `rig proj sync`.

---

## Requirements

- **rig** ≥ 0.11.1, which installs and manages R and the packages. You do not
  need to install R separately.
- **Git LFS** to download the pre-computed bootstrap data files

### Install rig and R (once per machine)

See [rig's install page](https://rig.r-lib.org/install.html). The user-level
installers below need no administrator rights; any other method on that page
(installer, winget, Homebrew, ...) also works, including with admin rights.

```bash
# macOS / Linux
curl -LsSf https://r-lib.github.io/rig/install.sh | sh
```

```powershell
# Windows (PowerShell); then open a new terminal
irm https://r-lib.github.io/rig/install.ps1 | iex
```

Then make R and packages install into your home folder, without admin rights:

```bash
rig system user-mode
```

On Windows PowerShell only, also run `rig system fix-r-alias` (otherwise `r`
means `Invoke-History`), then open a new terminal. Details: the
[rig user-mode tutorial](https://rig.r-lib.org/tutorial-user.html).

### Install Git LFS (once per machine)

```bash
git lfs install
```

### Managing packages

- `rig proj add <pkg>` adds a package (`--dev` for test and lint tooling
  only); `rig proj remove <pkg>` removes one. Both update `rproj.toml` and
  `rproj.lock`.
- Commit both files. Never commit `.rvenv/`: it is machine-specific and
  rebuilt by `rig proj sync`.
- `rig proj lock --upgrade` updates package versions; the golden and ledger
  tests must still pass afterwards.

> **Note on `hypogeom`:** This package was written specifically for this project by
> co-author Damian Pavlyshyn and is only available from GitHub
> (`damian-t-p/hypogeom`). `rig proj sync` installs it from GitHub at the commit
> pinned in `rproj.lock`.

---

## Repository layout

```
.
├── R/                      # Shared model core (sourced, not run directly)
│   ├── setup.R             # Libraries, seeds, constants, paper_path(), default parameters
│   ├── model_core.R        # Markov transitions, US detection, adherence model
│   ├── simulation.R        # Stochastic (binomial) version of the cohort model
│   ├── scenarios.R         # Scenario definitions as data; run_scenarios()
│   ├── bootstrap.R         # Calibration data, likelihood, bootstrap fitting
│   ├── plot_theme.R        # theme_hcc(), save_figure(), interval helpers
│   ├── plot_sweeps.R       # Bootstrap sweeps, band plots, stage distributions
│   ├── plot_scenarios.R    # Tornado, scenario and response-curve plots
│   └── plot_calibration.R  # Bootstrap calibration figures
├── papers/
│   ├── 2024-cancers/       # Paper 1 (published): see its README
│   │   ├── run_all.R           # Reproduces all of the paper's figures
│   │   ├── analysis/           # Numbered analysis scripts
│   │   ├── data/               # Parameter CSVs
│   │   ├── data-derived/       # Bootstrap data (.rds)
│   │   ├── output/             # manuscript_values.csv (git-ignored)
│   │   ├── figures/            # Generated figures (git-ignored)
│   │   └── releases/           # Release notes for the v1.x tags
│   └── stochastic/         # Paper 2 (in preparation)
├── tests/testthat/         # Test suite: core invariants and paper 1's published values
├── tools/                  # refit_calibration.R, lint_ratchet.R, check_globals.R
├── CHANGES.md              # Every change that affects results, with before/after values
└── CITATION.cff            # Machine-readable citation
```

### How the core finds a paper's files

The core never hard-codes a paper's folder. It reads and writes paper files
through `paper_path()` (in `R/setup.R`), which resolves to
`papers/<id>/...` for the paper chosen with `options(hcc.paper = "<id>")`.
Every script in a paper folder sets this option before it sources
`R/setup.R`. If the option is not set, the default is `2024-cancers`, the
published model, which is also what the test suite checks.

### Adding a paper

1. Create `papers/<id>/` with `data/`, `analysis/`, a `run_all.R` and a
   `README.md`, following `papers/2024-cancers/`.
2. Start every script with `options(hcc.paper = "<id>")`, then source the core
   files from `R/`.
3. Extend the core in `R/` so that paper 1 remains a special case. For example,
   a new compartment whose flows are zero gives back the published 7-state
   chain. The golden and paper-anchor tests then show that paper 1 is
   unchanged. If a change to the core does move paper 1's results, record it
   in `CHANGES.md` as a gate, as G1-G9 were.
4. When the paper is submitted, make a release (see [Releasing](#releasing)).

---

## Interactive use

To explore or modify the model outside the pipeline, source the two core files and call functions directly:

```r
source("R/setup.R")      # loads packages; reads the paper's Parameters.csv into `params_def`
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
all default to the point estimates in `params_def` (read from the active paper's `data/Parameters.csv`; by default `papers/2024-cancers/data/Parameters.csv`).
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
| `record_values(values, path)` | setup.R | Merge named values into the paper's `output/manuscript_values.csv` |
| `calibration_data` | bootstrap.R | Observed calibration data (Supp. Sec. 2.1.1) |
| `likelihood_matrix(p_0, p_A, sens_0, sens_A)` | bootstrap.R | The 15x15 transition matrix of Supp. Table S2 / matrix S1 |
| `end_state_probs(tm, sens_0, t)` | bootstrap.R | Probability of ending in each state after `t` months |
| `neg_log_lik(par, ...)` | bootstrap.R | Negative log-likelihood of the 8 calibrated parameters |
| `fit_point_estimate(start)` | bootstrap.R | Maximum-likelihood point estimate (thin wrapper of `fit_mle()`: Nelder-Mead then bounded L-BFGS-B polish) |
| `bootstrap_inputs_long()`, `bootstrap_estimates()` | bootstrap.R | Read the shipped bootstrap inputs and MLEs |
| `simulate_cohort(Pr_Att, ...)`, `simulate_by_adherence(...)` | simulation.R | Stochastic (binomial) check of the cohort model; not in the publication |
| `theme_hcc()`, `save_figure()`, `sweep_*()`, `plot_*()` | plot_*.R | Figure helpers: data builders, plot functions and saving |

---

## Contributing

### First-time setup

After cloning, activate the versioned git hooks and install the [air](https://github.com/posit-dev/air) R formatter:

```bash
git config core.hooksPath .githooks
# Make `git blame` skip bulk formatting commits (listed in .git-blame-ignore-revs).
git config blame.ignoreRevsFile .git-blame-ignore-revs

# Install air 0.13.0 (the version pinned in CI) if not already present.
# macOS:
brew install air
# Other OSes: see https://posit-dev.github.io/air/ (install page) and install
# version 0.13.0 so local formatting matches CI.
```

### What the hooks enforce

| Hook | What it does |
|------|-------------|
| `pre-commit` | Runs `air format` on staged `.R` files and re-stages them automatically, then runs the lint ratchet on them and blocks the commit if any file has more lints than its baseline; finally `tools/check_globals.R` (codetools) blocks the commit if any project function uses an undefined name |
| `commit-msg` | Enforces [git3moji](https://robinpokorny.github.io/git3moji/) prefixes: `:zap:` `:bug:` `:tv:` `:cop:` `:abc:` |

Format any file manually with `air format <file.R>`, or format the whole project with `air format .`.

### Lint and format

Formatting is owned by [air](https://github.com/posit-dev/air) (`air.toml`): line width 88, LF line endings, and `persistent-line-breaks = false`, so any call that fits on one line is collapsed onto it however it was written. `tribble()` calls are aligned as tables (`default-table = true`); write one row per line and air pads the columns. `air format .` formats, `air format --check .` checks.
Enable the versioned hooks once per clone with `git config core.hooksPath .githooks`.
`tools/check_globals.R` uses codetools to find any project function that uses an undefined name (a typo or a missing global); the hook and CI both run it.

Lint rules live in `.lintr` (line length 88 to match `air.toml`; `commas_linter` is off because air's aligned `tribble()` tables put spaces before commas, and air already enforces comma spacing everywhere else; `object_usage_linter` is off because it cannot see across the sourced `R/` files, and `tools/check_globals.R` does that job; `object_name_linter` and `indentation_linter` are disabled because the code follows paper notation such as `p_HCC` and air owns formatting).
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
Never run a paper's `run_all.R` as a test; it takes about 1-2 minutes (without `--bootstrap`).
The tests fall into four groups:

- **Golden** (`test-golden.R`, `tests/testthat/fixtures/golden.rds`): pins the code as it is, bugs included. The fixture is generated once by `fixtures/make_golden.R` and regenerated only under an approved, recorded change (see `CHANGES.md`). Doubles are compared at relative tolerance **1e-12**: portable across platforms (libm/BLAS noise is about 1e-15), while real model changes are at least 1e-3.
- **Invariants** (`test-invariants.R`, `test-core.R`, `test-simulation.R`, `test-likelihood.R`): properties that must hold whatever the numbers are, e.g. rows of the transition matrix sum to 1, people are conserved (`sum(dx) + sum(undiagnosed) == n`), the stochastic model agrees with the deterministic one in expectation, and the likelihood matches its closed form.
- **Paper anchors** (`test-paper.R`, `test-scenarios.R`): published values (abstract, Sec. 3, Table 3) after rounding to the printed precision.
- **Data and ledger** (`test-data.R`, `test-ledger.R`, `test-plotting.R`): CSV structure (including that `Parameters_sources.csv` lists the same parameters in the same order), ledger format, plot helpers.

### Continuous integration

`.github/workflows/check.yml` runs on every push and pull request (Ubuntu and Windows, R and packages from `rproj.lock` via rig).
It runs `air format --check .` (air pinned to the version in the workflow), the lint ratchet, `tools/check_globals.R` and the test suite.

### Releasing

Make a release whenever a paper is submitted, accepted or corrected:

1. Add an entry to `CHANGES.md`. In `CITATION.cff`, set `version` and
   `date-released`, and bump `version` in `rproj.toml`.
2. Write release notes in `papers/<id>/releases/<tag>.md`, including how to
   regenerate any `.rds` files that are in Git LFS.
3. Tag the commit and publish the release:
   ```bash
   git tag -a vX.Y.Z -m "<one line>"
   git push origin vX.Y.Z
   gh release create vX.Y.Z --title "<title>" --notes-file papers/<id>/releases/vX.Y.Z.md
   ```
   In the paper's code-availability statement, give the repository URL and
   the tag, e.g. "https://github.com/jaycrick/HCC-Model, release v2.0.0".
4. When the paper is published, add its DOI to the version table above and to
   the release notes on GitHub (release notes can be edited; tags cannot). In
   `CITATION.cff`, make the new paper the `preferred-citation` (the one
   GitHub's "Cite this repository" button shows) and move the previous one
   under `references:`.
5. Never move or delete a tag. If a released version needs a fix that `main`
   should not get, create a maintenance branch from its tag
   (`git switch -c paper1-maint v1.2.0`), fix it there and release
   `v1.2.1` from that branch. Only do this when such a fix is needed.

Optional backup: after a release, use "Save code now" at
[Software Heritage](https://archive.softwareheritage.org/save/) on the
repository URL. It archives the full history and every tag independently of
GitHub, free and without an account (Git LFS content excluded).

The rule for `.gitattributes`: Git LFS holds only files that can be
regenerated. Any file that cannot be regenerated must be a plain git file so
that release downloads contain it.

---

## Citation

Please cite the paper whose results you use. For paper 1:

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

