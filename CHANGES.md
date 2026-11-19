# Changes

## Structural changes (no numerical effect; S11 below 1e-12, see entry), S1-S11

- S1: characterisation, invariant and paper-anchor test suite (golden values
  in `tests/testthat/fixtures/golden.rds`); lint ratchet, pre-commit check and
  GitHub Actions CI.
- S2: core model on one named parameter list read by name from
  `data/Parameters.csv` (`load_params()`), not by position.
- S3: intervention scenarios defined once as data (`R/scenarios.R`).
- S4-S5: plotting helpers split into data builders, plot functions and
  `save_figure()`; `tools/check_globals.R` (undefined names); one
  `prepare_draws()` for every sweep.
- S6: manuscript numbers recorded in one ledger,
  `output/manuscript_values.csv` (`record_values()`); `data-derived/` is
  read-only input.
- S7: likelihood and calibration restructured (`likelihood_matrix()`,
  `end_state_probs()`, `neg_log_lik()`, `fit_point_estimate()`).
- S8: stochastic simulation mirrors the cohort model (C10, test-only effect).
- S9: model functions documented with paper and supplement references.
- S10: lint baseline is 0 (`tools/lint_ratchet.R`); no function over 60 lines
  (`neg_log_lik()` gained the helper `loglik_medians()`); `R/plotting.R` split
  into `R/plot_theme.R`, `R/plot_sweeps.R`, `R/plot_scenarios.R` and
  `R/plot_calibration.R`; repeated literal vectors replaced by the constants
  `STAGE_IDS`, `STAGE_P_NAMES`, `STAGE_SYM_NAMES`, `STAGE_SENS_NAMES`,
  `STAGE_L_NAMES`, `ADHERENCE_PARAMS`, `FIGURE_FORMATS`, `FIGURE_DPI_ALL`,
  `NM_START`; `renv.lock` no longer records ggsci and mc2d. No golden value
  changed.
- S11: vectorised cohort engine. `screen_cohort_mat()` runs N parameter sets at
  once on an N x 7 state matrix, with element-wise updates in place of
  `vec %*% tmat`. The three adherence groups are stacked as 3N rows
  (`screen_by_adherence_mat()`), and `life_add_vec()` gives life-years per
  draw. `screen_cohort()`, `screen_by_adherence()` and `life_add()` are thin
  N = 1 wrappers with unchanged signatures. `sweep_model()` makes one call per
  grid value instead of one per draw, and `run_scenarios()` one call per
  scenario x risk.
  - Numerical effect: NOT bit-identical. On this machine (arm64, reference
    BLAS) `%*%` rounds differently from the element-wise sum (probably FMA
    contraction in the BLAS; not verified). Measured: the `life_add()`
    baseline golden is bit-identical; the other goldens differ by up to 3.6e-16
    relative, the full-size sweep objects by up to 6.5e-16 and
    `run_scenarios()` by up to 5.8e-16. In `output/manuscript_values.csv`, 11
    of 85 values change, by at most 4.3e-15 relative. The largest are
    differences of nearby quantities (cancellation), with absolute differences
    about 1e-14. All are far below the 1e-12 gate. Two JPGs change only in
    anti-aliasing pixels (attendance by p_CA, increasing adherence); the PDFs
    already differed between identical runs (embedded timestamps). Before this
    change the ledger already differed by 1 ulp from an earlier run
    (`real_adh_lower_high_risk`).
  - Timing (arm64 macOS, R 4.6):

    | step | before | after |
    |---|---|---|
    | analysis/03 | 1778 s (under Rprof) | 64 s |
    | `Rscript run_all.R` (no `--bootstrap`) | ~30-37 min | 71 s |
    | `run_scenarios()`, 200 draws | 4.8 s | 0.18 s |
    | `run_scenarios()`, 10,000 draws (01 `--bootstrap`) | ~4 min (extrapolated) | 6.7 s |
    | `fit_point_estimate()` | 4.0 s | 4.0 s (unchanged) |

    `fit_point_estimate()` spends 90% in hypogeom `cfgeom`/fft via
    `DescTools::dOrder` (external code).

# Changes that alter model output

## G1 - untreated progression double conversion (finding C1)

`fromf4()` (now `progress_untreated()`) converted `p_0..p_D` from annual to
monthly, but those parameters are already monthly (Cumming et al. 2024,
Cancers 16:2745, Table 1; `data/Parameters.csv`; the conversion was removed
from the screening engine in 2023 but not here). The untreated-progression
figures (`risk_no_int`, `years_no_int`) were therefore about 12x too slow.
Published scenario results (`life_add()` and friends) are unaffected.

`progress_untreated(n = 100, t = 120)`, before and after:

| State | Before | After |
|---|---|---|
| Cirrhosis | 80.05501585576 | 80.0550158558 |
| BCLC 0 | 13.38155807791 | 1.9015484395 |
| BCLC A | 5.40591015046 | 2.8226122219 |
| BCLC B | 0.50530853727 | 0.4573605702 |
| BCLC C | 0.30876059923 | 0.4483577380 |
| BCLC D | 0.25227005443 | 0.9921202126 |
| Mortality | 0.09117672493 | 13.3229849620 |

Test: `progress_untreated(n, t)` must equal the undiagnosed states plus deaths
of `screen_cohort(Pr_Att = 0, Sym_* = 0)` (tests/testthat/test-core.R). Only the
`progress_untreated_100_120` entry of `tests/testthat/fixtures/golden.rds` was
regenerated (key renamed from `fromf4_100_120` earlier in the same step with its
value untouched).

## G2 - p_B axis label (finding C5)

`TeX_labels["p_B"]` was built from `c("p", "C")` and rendered as p_C on every
tornado plot that labels parameters with it. The special case is removed, so
every label is built from its own parameter name.

| | Before | After |
|---|---|---|
| `deparse(TeX_labels[["p_B"]])` | `italic(p)[C]` | `italic(p)[B]` |

Test: for every entry of `TeX_labels` the deparsed plotmath equals the one
built independently from the parameter name (tests/testthat/test-plotting.R);
three expectations failed before the fix. No number changes.

## G3 - "Upper p_HCC" tornado (finding C6)

`sens_anal_f4_0_upper.jpg` plotted `years_mid_comp` (and sorted by
`years_low_comp`) instead of `years_high_comp`; its `Type` case_when matched
nothing and is deleted. Now `plot_bound_tornado(param_tbl, "high")` plots and
sorts by `years_high_comp`, as the "low" plot does with `years_low_comp`.

| Largest bars | Before (years_mid_comp) | After (years_high_comp) |
|---|---|---|
| Sens_0 Upper | 33.617 | 56.350 |
| p_0 Lower | 20.223 | 34.683 |
| Sens_0 Lower | -19.134 | -32.159 |
| p_0 Upper | -17.183 | -29.241 |

Test: the plotted x of the bar layer equals `years_high_comp` and the factor
order is the abs-sum order of `years_high_comp`; failed before.
The `sens_anal_25_f4_0_upper/lower.jpg` figures were already correct and are
byte-identical to the baseline.

## G4 - Sym_B in the bootstrap sweeps (finding C7)

The sweeps (`usfreq_graph`, `pafunc`, the attendance analysis in 03 and
`p_Att_Stages`) passed each draw's Sym_A but the default Sym_B, while the
scenarios use Sym_B = `SYM_B_RATIO` * Sym_A (8x). One helper,
`prepare_draws(df)` (R/model_core.R), now sets it per draw; `run_scenarios()`
and every sweep (`sweep_model()` in R/plot_sweeps.R) use it. The point-estimate
line is unchanged (default parameters; tested identical).

Bands, all 10 000 draws (scratch script, 2.5 / 97.5 %; IQR for stages):

| Sweep | x | Point | Before | After |
|---|---|---|---|---|
| usfreq m | 1 | 129.841 | 114.794, 142.003 | 114.783, 141.990 |
| usfreq m | 6 | 68.023 | 55.882, 78.642 | 55.881, 78.609 |
| usfreq m | 12 | 48.907 | 39.929, 56.577 | 39.929, 56.593 |
| usfreq m | 24 | 35.210 | 28.939, 40.321 | 28.785, 40.591 |
| p_CA | 0 | 53.476 | 43.747, 61.657 | 43.723, 61.656 |
| p_CA | 0.336 | 67.726 | 55.628, 78.305 | 55.632, 78.271 |
| p_CA | 0.685 | 82.571 | 68.114, 95.599 | 68.044, 95.603 |
| attendance p_Att | 0 | 14.857 | 11.380, 17.768 | 10.675, 18.825 |
| attendance p_Att | 0.5 | 76.049 | 61.746, 88.767 | 61.780, 88.727 |
| attendance p_Att | 1 | 113.683 | 93.312, 132.291 | 93.263, 132.288 |
| stage BCLC B, p_Att 0 | 0 | 2.661 | 2.080, 3.191 | 2.280, 3.032 |
| stage BCLC C, p_Att 0 | 0 | 5.792 | 5.105, 6.364 | 5.138, 6.316 |
| stage BCLC B, p_Att 0.5 | 0.5 | 2.270 | 1.844, 2.621 | 1.932, 2.586 |
| stage BCLC D, p_Att 1 | 1 | 0.420 | 0.334, 0.497 | 0.334, 0.494 |

(BCLC 0 and A bands are unchanged at the shown points: Sym_B affects only
later stages. All 15 stage rows are in the scratch output.)

Tests: `prepare_draws` sets Sym_B; a sweep on 5 draws equals `life_add` with an
explicit Sym_B (m, p_CA, p_Att and the stage sweep); `run_scenarios` shares it.

## G5 - labels (findings C8, C9) and file name

(a) Ribbon legends said "50%" for the 2.5-97.5 % band; they now come from
`band_label(probs)`: "95%". `p_Att_Stages` uses the 25-75 % band (quantiles
unchanged); its ribbon figure had no band label, so the caption now says
"Ribbon: IQR across bootstrap draws". Tests: `band_label`, legend fill, caption.

(b) The facet level "Pr" never occurred (`word("p_NA", "_")` is "p"), so the
adherence parameters sat under "Monthly Risk". One lookup `param_lookup`
(parameter, group, label) in R/setup.R gives each parameter exactly one group;
`variable_labels`, `facet_labels` and the facet order are derived from it, and
the Adherence facet (p_NA, p_CA, p_Att) now comes first. Tests: every parameter
has one group; p_NA, p_CA, p_Att are Adherence; labels derive from the lookup.
Affects sens_anal_univariate, sens_anal_f4_0_both, sens_anal_mod and
sens_anal_25_univariate.

(c) `Untreated_HCp_Cistribution_by_Risk_10_Year.jpg` is now
`Untreated_HCC_Distribution_by_Risk_10_Year.jpg`.

## G6 - `mlen` for even n (closed, not a bug)

`calibration_data$stage$mlen` is n - 1 for even n (12 -> 11, 138 -> 137,
210 -> 209), as Supp. p.3 footnote 2 prescribes. No code change, no number
changes.

## G7 - stale hard-coded point estimates in analysis/01

`analysis/01_bootstrap_params.R` hard-coded `hypo_probs`, which was stale.
`data/Parameters.csv` POINT equals the paper's Tables 1-2, and equals the
Nelder-Mead MLE from `optim()` in 01 (max |diff| 4.6e-9). 01 now takes the
point estimates from the CSV (`params_def`) and compares the optimiser's result
to them instead of overwriting them.

| | Before (`hypo_probs`) | After (CSV) |
|---|---|---|
| NLL | 12.587514779190872 | 12.549739472600052 |
| p_0 | 0.07982013 | 0.07982013 |
| p_A | 0.05540597 | 0.05540597 |
| p_B | 0.33922926 | 0.34353183 |
| p_C | 0.35681215 | 0.35197652 |
| p_D | 0.16023503 | 0.16046211 |
| Sens_0 | 0.21878495 | 0.20916773 |
| Sens_A | 0.54017653 | 0.53267066 |
| Sym_A | 0.01015650 | 0.01015362 |

Downstream uses of the point estimates in 01 that change: the red dots / lines
in `HCC_Rates_Range_Plot`, `HCC_Rates_Hist`, `HCC_Rates_Violin`; the points in
`HCC_Time_Range_Plot` and the panel of `multipanel_cals`; the cross in
`US_sensitivity_Cal` (this already used the CSV); the ledger names
`mean_0` .. `mean_D` (the `_Lower` / `_Upper` ones come from the bootstrap and
do not change); and the `point` column of the results table (point-estimate
scenarios). Ledger values `mean_*` (years in stage, 1/p/12), before -> after:

| | Before | After |
|---|---|---|
| mean_0 | 1.0440140016 | 1.0440140016 |
| mean_A | 1.5040497140 | 1.5040497140 |
| mean_B | 0.2456549100 | 0.2425782011 |
| mean_C | 0.2335495956 | 0.2367582171 |
| mean_D | 0.5200693839 | 0.5193334011 |

Tests: the golden entry `neg_log_lik_hypo_probs` (12.587514779190872) became
`neg_log_lik_csv` (value above); it is the only golden entry regenerated.
`test-likelihood.R` checks that Nelder-Mead from the documented start
reproduces the CSV (1e-6 per parameter, 1e-8 in NLL) and that NLL(CSV) <
NLL(old `hypo_probs`).

### P7 - the published MLE is Nelder-Mead's stopping point

The CSV / paper values are where Nelder-Mead stops at its default tolerance,
not an exact optimum. A tight BFGS on the logit scale reaches NLL 12.5495403,
2e-4 lower (max |delta logit| 4.6e-3; p_B 0.34353 vs 0.34307). That is
negligible next to the bootstrap intervals, so the published values stand.

### Other likelihood notes (S7)

- `mlen` is n - 1 for even n (12 -> 11, 138 -> 137, 210 -> 209), per Supp. p.3
  footnote 2; not a bug (G6 closed).
- Supp. p.6 prints Pr(d_alpha) without the (1 - Sym) factors; the code's form
  is correct (0.3289 vs 52/158 = 0.329 observed).
- Pr(d_0) uses T^100 for the limit. Error against the exact closed form, and
  against T^1000: 1.445e-6 (below the 1e-5 tested).
- The 15x15 matrix is built by `likelihood_matrix()` with named states.

## Final sweep (S10): ledger and figure diff against the baseline

`Rscript run_all.R` (36 min 36 s, exit 0), plus `analysis/01` run without its
two `saveRDS()` calls so `data-derived/` stays untouched. The 33 ledger values
written by 01 are compared too. `output/manuscript_values.csv` against the
frozen baseline (relative tolerance 1e-12): 85 names compared, 3 differ, all
from G7.

| Name | Before | After | Gate |
|---|---|---|---|
| mean_B | 0.245654909996 | 0.242578201075 | G7 |
| mean_C | 0.233549595588 | 0.236758217092 | G7 |
| mean_D | 0.520069383913 | 0.519333401096 | G7 |

The baseline holds 6 further names (`Pr_Com_Adh_Lower/Upper`,
`Pr_Non_Adh_Lower/Upper`, `Pr_PA_Att_Lower/Upper`) that no code in this repo
writes. They are stale entries from the pre-publication naming, which was
renamed to `p_CA`/`p_NA`/`p_Att` in the original Dropbox repository (commit
e72e6aa). The old read-modify-write of `param_names.rds` never deleted them.
The same quantities are written as `p_CA_*`, `p_NA_*`, `p_Att_*`, and those
match the baseline `p_*` entries.

Figures (files written by this run): the same set as the baseline except
G5c (`Untreated_HCp_Cistribution_...` is now `Untreated_HCC_Distribution_...`).
19 files changed hash; 3 of them are metadata only.

| Figures | Gate |
|---|---|
| `scenarios_new.pdf`, `scenarios_new_low_risk.pdf`, `scenarios_new_high_risk.pdf` | none: only /CreationDate and /ModDate differ (12 bytes; same size, same text, identical render at 100 dpi) |
| `Untreated_HCC_over_time_10_Year.jpg`, `Untreated_HCC_Distribution_by_Risk_10_Year.jpg` | G1 (and G5c name) |
| `attendance_plot_10_year.{jpg,png,pdf}`, `attendance_plot_stages_ribbon_10_year.jpg`, `HCC_Years_of_Life_Saved_by_p_CA_10_Year.jpg`, `HCC_Years_of_Life_Saved_by_Testing_Frequency_10_Year.jpg`, `..._10_log_Year.jpg` | G4 (bands), G5a (legend/caption) |
| `sens_anal_univariate.jpg`, `sens_anal_mod.{jpg,png,pdf}`, `sens_anal_25_univariate.jpg`, `sens_anal_f4_0_both.jpg` | G5b (facets), G2 (p_B label) |
| `sens_anal_f4_0_upper.jpg` | G3 |
| `sens_anal_f4_0_lower.jpg` | G2 |
