# Paper 1: Cumming et al. 2024, *Cancers*

> Cumming J, Scott N, Howell J, Flores JE, Pavlyshyn D, Hellard ME, Winata LS, Ryan M, Sutherland T, Thompson AJ, Doyle JS, Sacks-Davis R.
> **Improving Hepatocellular Carcinoma Surveillance Outcomes in Patients with Cirrhosis after Hepatitis C Cure: A Modelling Study.**
> *Cancers* 2024, **16**(15), 2745. https://doi.org/10.3390/cancers16152745

This folder holds everything specific to this paper: its parameter inputs, the
bootstrap data, the analysis scripts and the outputs. The model code they use
is the shared core in [`R/`](../../R) at the repository root.

> [!IMPORTANT]
> **Results on `main` differ from the publication.** In November 2026 the
> calibration was redone with the exact likelihood and a better optimiser, and
> all 10,000 bootstrap resamples were refitted (gates G8, G9). The scenario
> point estimates still match the paper (9.8, 17.2, 7.0, 9.0 life-years per
> 100), but several confidence intervals and the partial correlations do not.
> See [Errata and differences from the publication](#errata-and-differences-from-the-publication).
> To run the code exactly as published, use release
> [`v1.0.0-published`](https://github.com/jaycrick/HCC-Model/releases/tag/v1.0.0-published)
> (original scripts) or
> [`v1.1.0`](https://github.com/jaycrick/HCC-Model/releases/tag/v1.1.0)
> (refactored and tested; reproduces the published numbers). The version table
> is in the [root README](../../README.md#versions-and-archives).

---

## Layout

```
papers/2024-cancers/
├── run_all.R                  # Runs the analyses in order → figures/
├── analysis/                  # Numbered analysis scripts
│   ├── 01_bootstrap_params.R   # Scenarios + figures from the bootstrap draws (skipped by default)
│   ├── 02_sensitivity.R        # Univariate sensitivity analysis
│   ├── 03_progression_model.R  # Progression and screening-frequency figures
│   └── 04_clinical_scenarios.R # Clinical scenario comparisons
├── data/                      # Read-only parameter inputs (CSV)
├── data-derived/              # Bootstrap data (.rds; see "Regenerating the .rds files")
├── output/                    # manuscript_values.csv (git-ignored; written by the scripts)
├── figures/                   # Generated figures (git-ignored)
└── releases/                  # Release notes for the v1.x tags (GitHub releases)
```

Every script sets `options(hcc.paper = "2024-cancers")` before it sources
`R/setup.R`, so the shared core reads and writes this folder through
`paper_path()`.

---

## How to reproduce the figures

From the repository root (see the root README's Quick start for setup):

```bash
Rscript papers/2024-cancers/run_all.R
```

Figures are written to `papers/2024-cancers/figures/`. The run takes about two
minutes.

### Run order

The scripts must run in numerical order because later scripts read from
`data-derived/`:

| Script | Reads | Writes |
|--------|-------|--------|
| `01_bootstrap_params.R` | `resampled_pars.rds`, `bootstrap_fit_diagnostics.rds`, `list_of_resampled_*.rds` | `resampled_pars_fin.rds`, ledger, figures |
| `02_sensitivity.R` | `Parameters.csv`, `Parameters_25.csv` | ledger, figures |
| `03_progression_model.R` | `resampled_pars_fin.rds` | figures |
| `04_clinical_scenarios.R` | `resampled_pars_fin.rds` | ledger, figures |

"Ledger" is `output/manuscript_values.csv` (see [Outputs](#outputs)).

`run_all.R` enforces this order. It skips `01_bootstrap_params.R` by default
because its output is shipped through Git LFS. Pass `--bootstrap` to run it:

```bash
Rscript papers/2024-cancers/run_all.R --bootstrap
```

The expensive calibration refit (point MLE plus 10,000 bootstrap resamples) is a
separate one-off step, `tools/refit_calibration.R` at the repository root. It
writes `data-derived/resampled_pars.rds`, `data-derived/bootstrap_fit_diagnostics.rds`
and the calibrated rows of `data/Parameters.csv` and `data/Parameters_25.csv`:

```bash
Rscript tools/refit_calibration.R [--n 20] [--cores 10] [--config a|b|c|d] [--dry-run] [--out DIR]
```

`--dry-run` writes only to `--out` (default: a temporary directory).

---

## Regenerating the .rds files

A release's source download (zip or tar.gz from the GitHub release page) does
not contain the files stored in Git LFS: those are replaced by small text
pointer files.
Only files that can be regenerated are stored in LFS; the three that cannot are
plain git files and are in every archive from `v1.3.0` on. A `git clone`
followed by `git lfs pull` gets every file, so none of this is needed when
working from a clone.

| File | In archive? | How to get it |
|------|-------------|---------------|
| `data-derived/list_of_resampled_milton_output.rds` | Yes | Cannot be regenerated (HPC bootstrap resamples, 2024)\* |
| `data-derived/miltonoutput.rds` | Yes | Cannot be regenerated (HPC bootstrap MLEs, 2024; starting values for the refit) |
| `data-derived/param_names.rds` | Yes | Cannot be regenerated (frozen published values) |
| `data-derived/resampled_pars.rds` | No | Step 1 |
| `data-derived/bootstrap_fit_diagnostics.rds` | No | Step 1 |
| `data-derived/resampled_pars_fin.rds` | No | Step 2 |
| `tests/testthat/fixtures/golden.rds` | No | Step 3 |
| `tests/testthat/fixtures/ledger_baseline.rds` | No | Step 3 (after step 2) |

\* The resampling code is in `bootstrapped_param_calcs.R` of
`v1.0.0-published` (commented out, `set.seed(314)`). Rerunning it reproduces
the first five of the 10,000 resamples exactly and then diverges, because the
exact probability vector used for the HPC run was not recorded.

From the repository root, with the packages from `rproj.lock`
(`rig proj sync`):

```bash
# 1. Calibration refit: point MLE + all 10,000 bootstrap resamples.
#    Writes data-derived/resampled_pars.rds and bootstrap_fit_diagnostics.rds,
#    and rewrites the calibrated rows of data/Parameters.csv and
#    data/Parameters_25.csv (these come out byte-identical to the committed files).
#    No random numbers; the result does not depend on --cores.
Rscript tools/refit_calibration.R --cores 10

# 2. Scenario life-years for every draw (set.seed(315)) -> resampled_pars_fin.rds,
#    then all figures and the full ledger output/manuscript_values.csv
Rscript papers/2024-cancers/run_all.R --bootstrap

# 3. Test fixtures, then the tests
Rscript tests/testthat/fixtures/make_golden.R
Rscript tests/testthat/fixtures/make_ledger_baseline.R
Rscript -e 'testthat::test_dir("tests/testthat")'
```

These commands were checked by regenerating everything from the three inputs
with the v1.2.0 code, started from the same fixed point (bit-identical start,
see `CHANGES.md` S13). Every file and both CSVs were reproduced exactly;
the diagnostics file differs only in the run timestamp. The refit takes about
2.7 h on 10 cores. Details in [`releases/v1.2.0.md`](releases/v1.2.0.md).

The older releases store every `.rds` in LFS. Their release notes
([`releases/`](releases)) give the commands for each one, including which
archive to download the irreplaceable files from.

---

## Data provenance

| File | Source | Can be regenerated? | Description |
|------|--------|---------------------|-------------|
| `data/Parameters.csv` | Cumming 2024, Tables 1-3; calibrated rows from `tools/refit_calibration.R` (G9) | n/a (plain git) | Point estimates and 95% bounds for all model parameters |
| `data/Parameters_sources.csv` | Cumming 2024, Tables 2-3 | n/a (plain git) | `PARAMETER`, `SOURCE`, `TABLE_REF` for each row of `Parameters.csv`; unknowns are marked `TODO: confirm with author` |
| `data/Parameters_25.csv` | Derived from `Parameters.csv` | n/a (plain git) | ±25% variation bounds for modifiable parameters |
| `data-derived/list_of_resampled_milton_output.rds` | Parametric bootstrap, HPC cluster (2024) | **No** (plain git) | Bootstrapped summary statistics, the inputs to the bootstrap MLEs |
| `data-derived/miltonoutput.rds` | HPC cluster MLE run (2024) | **No** (plain git) | Original HPC bootstrap MLEs; superseded by the refit (G9) but used as its starting values |
| `data-derived/param_names.rds` | Original 2024 pipeline | **No** (plain git) | Frozen reference: the values as published. Never rewritten |
| `data-derived/resampled_pars.rds` | `tools/refit_calibration.R` (G9) | Yes (Git LFS) | Refitted bootstrap MLEs (8 transition and sensitivity parameters, 10,000 rows) |
| `data-derived/bootstrap_fit_diagnostics.rds` | `tools/refit_calibration.R` (G9) | Yes (Git LFS) | Per-resample NLL, convergence code, best start, at-bound flags, old HPC NLL, audit |
| `data-derived/resampled_pars_fin.rds` | `analysis/01_bootstrap_params.R` | Yes (Git LFS) | `resampled_pars.rds` plus scenario life-year columns |

The test fixtures `tests/testthat/fixtures/golden.rds` and
`tests/testthat/fixtures/ledger_baseline.rds` (repository root) are also in Git
LFS and can also be regenerated.

---

## Key model outputs (regression anchors)

These values should be reproduced exactly from a clean run:

- `life_add()` (baseline, default parameters) ≈ **68.03421** life-years per cohort of 100
- `life_add(p_NA = 1, p_CA = 0)` (no surveillance) ≈ **14.86345** life-years

Both are checked with `stopifnot` in `04_clinical_scenarios.R`. The test suite
at the repository root (`tests/testthat/`) checks this paper's published values
(`test-paper.R`, `test-scenarios.R`), the golden values and the ledger on every
push.

---

## Outputs

- `output/manuscript_values.csv` (columns `name`, `value`) is the ledger of every number quoted in the manuscript, written by `record_values()`. It is sorted by name and written at a precision that round-trips exactly, so it is text-diffable.
- `data-derived/` is read-only input. Only the deliberate regeneration steps described above write to it.
- `figures/` is generated by the pipeline and git-ignored.

---

## Errata and differences from the publication

Cumming et al. 2024, *Cancers* 16(15):2745, and its Supplement. Full details,
before/after tables and the complete ledger diff are in [`CHANGES.md`](../../CHANGES.md) (gates G8, G9;
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
| P5 | Supplement Table S2 | Lists 14 states; the code adds a 15th absorbing state `past_A` (progressed past BCLC A undiagnosed). Also, Pr(d_0) uses the t^(100) approximation (error about 1.5e-6) rather than the exact closed form, because it generated the published MLE. **Closed by G8**: the code now uses exact absorption probabilities. See [`CHANGES.md`](../../CHANGES.md). |
| P6 | Paper Sec. 2.4, scenario 3 | The text lists L_A +1 year and L_B +6 months only; the code, and the published 9.0 life-years, also add L_C + 0.5. |
| P7 | Table 3 / `Parameters.csv` | The published MLE is Nelder-Mead's stopping point at its default tolerance, not an exact optimum (NLL 2e-4 above a tight BFGS optimum); negligible next to the bootstrap intervals. **Closed by G9**: `fit_mle()` polishes with a bounded L-BFGS-B, and the whole bootstrap was refitted locally. See [`CHANGES.md`](../../CHANGES.md). |
