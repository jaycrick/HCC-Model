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
- S12: per-paper layout and releases (after tag `v1.2.0`).
  - Tags `v1.0.0-published` (06703e6, the code and data as published),
    `v1.1.0` (6902e05, refactored, reproduces the publication) and `v1.2.0`
    (3eab3fd, G8 + G9). The root README has the version table.
  - Paper 1's files moved to `papers/2024-cancers/` (`analysis/`, `data/`,
    `data-derived/`, `output/`, `figures/`, `run_all.R`, and a README holding
    the paper-specific documentation and errata). `R/` is the shared core,
    which reads and writes paper files only through `paper_path()`, selected
    by `options(hcc.paper = "<id>")` (default `2024-cancers`).
    `papers/stochastic/` holds paper 2.
  - `.gitattributes`: the three files that cannot be regenerated
    (`list_of_resampled_milton_output.rds`, `miltonoutput.rds`,
    `param_names.rds`) are now plain git files, so GitHub release downloads
    contain them; the rest stay in Git LFS. Contents unchanged (same SHA-256).
  - Numerical effect: none. All tests pass unchanged, and the 58 ledger values
    written by `run_all.R` are byte-identical to the run before the move.
- S13: fixed starting point for the refit's point fit.
  - Problem: `tools/refit_calibration.R` started its point fit from the POINT
    column of `data/Parameters.csv`, which the same script rewrites. The
    shipped G9 refit started from the published MLE (the CSV before G9); a
    rerun started from the refitted CSV (start shifted by about 1e-9). In 253
    of the 10 000 bootstrap resamples two starts reach optima whose NLL
    differs by about 1e-11, so the tiny shift changed which start won and moved
    those fits by up to 2e-5 relative (NLL unchanged to 1.1e-8). Downstream,
    ledger values moved by up to 2.3e-3 relative (`pcor_logit_US_0`). The
    refit itself is deterministic: 1 and 9 cores, and repeated runs, give
    bit-identical results.
  - Fix: the point fit starts from `NM_START` and the constant
    `PUBLISHED_MLE` (`R/bootstrap.R`), the published MLE as stored in the CSV
    up to `v1.1.0`. It is bit-identical to the start the shipped refit used, so
    the shipped files are unchanged and a rerun now reproduces them.
  - Tests: `PUBLISHED_MLE` names and rounded values; the point fit from the
    fixed starts reproduces the CSV to 1e-7.

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

**Closed by G9** (below): the optimiser is now Nelder-Mead followed by a bounded
logit L-BFGS-B polish, and the CSV holds the polished MLE. The paragraph above
is kept as history.

### Other likelihood notes (S7)

- `mlen` is n - 1 for even n (12 -> 11, 138 -> 137, 210 -> 209), per Supp. p.3
  footnote 2; not a bug (G6 closed).
- Supp. p.6 prints Pr(d_alpha) without the (1 - Sym) factors; the code's form
  is correct (0.3289 vs 52/158 = 0.329 observed).
- Pr(d_0) uses T^100 for the limit. Error against the exact closed form, and
  against T^1000: 1.445e-6 (below the 1e-5 tested). **Closed by G8** (below):
  `end_state_probs()` now returns the exact absorption probabilities; the
  T^100 power survives only as `t = 100` in a test that documents this error.
- The 15x15 matrix is built by `likelihood_matrix()` with named states.

## G8 - exact absorption probabilities (finding P5)

`end_state_probs(tm, sens_0)` returned `t0 %*% T^100`, an approximation of the
limit (Supp. Sec. 2.1.3 p. 4, footnote 3). It now solves the absorbing chain
exactly: with transient states 3:14 and absorbing states 1, 2, 15,
B = (I - Q)^-1 R and the end-state vector is `t0[absorbing] + t0[transient] %*% B`.
A numeric `t` still gives the power form (used only by a test). Pr(d_0) equals
the closed form to 1e-12 and the T^1000 power to 1e-9; the probabilities sum to
1 to 1e-12.

NLL at the old point (`data/Parameters.csv` before G9): 12.5497394726 (T^100)
-> 12.5497425317 (exact), relative change 2.4e-7. Pr(d_0) changes by a relative
4.4e-6 and Pr(d_A) by 1.4e-5.

## G9 - optimiser and local bootstrap refit (finding P7)

**Method.** `fit_mle()` (`R/bootstrap.R`) runs Nelder-Mead (probability scale,
`maxit = 10000`, the old behaviour) and then polishes with L-BFGS-B on the
scaled-logit scale z = qlogis(par / ub) (`factr = 1e3`, `maxit = 1000`). A
non-finite NLL becomes 1e10 inside the polish only. The polish is kept only if
it lowers the NLL. The box (`FIT_LB`, `FIT_UB`) has these reasons:

- H1: the five p_* have lower bound 1e-3. A tiny p makes the hypogeometric FFT
  support huge (up to 8 s and many GB per evaluation at p = 1e-4; an unbounded
  BFGS from `NM_START` exhausted 18 GB). Sens_0, Sens_A, Sym_A have 1e-6.
- H2: Sym_A < 1/8 (upper bound (1/8)(1 - 1e-6)), because Sym_B = 8 Sym_A must
  be at most 1 (otherwise the NLL is NaN).
- H3: the NLL is Inf near p_D = 1, hence the 1e10 penalty.
- H4: fits can sit at p = 1, so the upper bound is 1 - 1e-6.

`at_bound` flags a parameter on or outside either limit (`fit_at_bound()`).
The point fit is multi-start: `fit_mle(list(NM_START, old CSV point))`. It
reproduces the new CSV to 5e-13 (NLL 12.54954; 12.54974 was the T^100 NLL at
the old point).

**Bootstrap refit** (`tools/refit_calibration.R`, config d): each of the 10,000
resamples is polished (no Nelder-Mead) from its old HPC estimate, clipped into
the box. If that fit does not converge, sits on the box or fails, the resample
is refitted with NM + polish from the new point and from the HPC estimate, and
the better result is kept (`fallback`). `best_start` is 1 (HPC polish), 2 or 3
(the fallback starts; see `fit_resample()`). A 200-resample audit reruns every
start of config a on evenly spaced resamples and compares. Runtime: 114 min on
10 cores (0.68 s per resample). Outputs: `data-derived/resampled_pars.rds`,
`bootstrap_fit_diagnostics.rds`, and the calibrated rows of `data/Parameters.csv`
and `Parameters_25.csv` (LOWER/UPPER = bootstrap 2.5% / 97.5%, Sym_B = 8 Sym_A).

**Point estimates and intervals**, old -> new:

| Parameter | Old POINT | New POINT |
|---|---|---|
| p_0 | 0.07982013 | 0.079927088 |
| p_A | 0.05540597 | 0.055445835 |
| p_B | 0.34353183 | 0.34306511 |
| p_C | 0.35197652 | 0.35207011 |
| p_D | 0.16046211 | 0.16056321 |
| Sens_0 | 0.20916773 | 0.209582 |
| Sens_A | 0.53267066 | 0.53149421 |
| Sym_A | 0.01015362 | 0.010156768 |
| Sym_B | 0.08122896 | 0.081254146 |

| Parameter | Old 95% CI | New 95% CI |
|---|---|---|
| p_0 | 0.03066517 to 0.24007789 | 0.030301295 to 0.4823508 |
| p_A | 0.04257023 to 0.08169524 | 0.042556688 to 0.081737821 |
| p_B | 0.2034034 to 0.7371478 | 0.20698394 to 0.88773261 |
| p_C | 0.2237918 to 0.5741617 | 0.20680464 to 0.51179899 |
| p_D | 0.1283536 to 0.1628741 | 0.13069987 to 0.16067135 |
| Sens_0 | 0.07997944 to 0.53655688 | 0.077870324 to 0.999998 |
| Sens_A | 0.3606599 to 0.8198396 | 0.35518875 to 0.81687098 |
| Sym_A | 0.00670664 to 0.0156842 | 0.0067225479 to 0.015943164 |
| Sym_B | 0.05365312 to 0.1254736 | 0.053780383 to 0.12754531 |

**Diagnostics.**

- Fallback used in 506 of 10,000 resamples. `best_start` was 1 in 9665, 2 in
  178, 3 in 157.
- Convergence code 52 (L-BFGS-B abnormal line-search termination) in 17
  resamples (0.17%); each is better than its HPC fit (by 4e-5 to 1.26 NLL).
- At the box: Sens_0 418, p_B 21, Sens_A 11, p_C 3, others 0.
- Old HPC NLL minus new NLL, under the new likelihood: minimum 1.8e-7, median
  0.027. The new fit is never worse.
- Audit (200 resamples, config a): max gap 7.9e-11.

**Boundary fits are genuine.** In resamples 42, 47, 71 and 77 the NLL profile
decreases monotonically up to Sens_0 = 1. These resamples have a low m_0
(21-29 against 38 observed), so p_0 is high (0.42-0.55). The old HPC fits there
had Sens_0 0.39-0.64 and an NLL 0.18-1.27 worse: the HPC bootstrap
under-converged in the tails. This is why the Sens_0 interval widens to
0.078-1.000 and the p_0 upper limit rises from 0.24 to 0.48.

**HPC reproduction caveat.** Nelder-Mead from `NM_START` does not reproduce the
HPC fits exactly (resample 1: max |diff| 0.0039; HPC NLL 12.87703 vs local NM
12.87689). The HPC start or settings are unknown, so the old HPC estimates are
used only as starting points and for the comparison above.

**Partial correlations (flag).** The published logit partial correlations of
Baseline with p_0, Sens_0 and Sens_A are unchanged as a method, but their
values moved: pcor_logit_0 -0.873 -> +0.288, pcor_logit_US_0 +0.920 -> -0.018,
pcor_logit_US_A 0.221 -> 0.187. The cause is that 460 resamples now have
|logit| > 8 (65 before), mostly Sens_0 at the upper bound with logit about 13,
which dominate the logit-scale correlation. `analysis/01` therefore also
records rank-based (Spearman) partial correlations on the raw parameters
(ranks do not depend on the logit): pcor_rank_0 -0.845, pcor_rank_US_0 +0.902,
pcor_rank_US_A +0.508, close to the published -0.87 and 0.92. The fit
diagnostics `n_at_bound_Sens_0` (418), `n_fit_nonconverged` (17) and
`n_fit_fallback` (506) are in the ledger too.

**Published numbers changed** (the anchors in `test-paper.R` now pin the refit;
point-estimate anchors that still round to the published value are unchanged):

| Published | G9 refit | Location in paper |
|---|---|---|
| 11.6 | 11.7 | Abstract / Sec 3.3, realistic 95% CI upper |
| 13.9 | 13.8 | Abstract / Sec 3.3, optimistic 95% CI lower |
| 20.3 | 20.5 | Abstract / Sec 3.3, optimistic 95% CI upper |
| 3.6 | 2.5 | Abstract / Sec 3.3, diagnosis 95% CI lower |
| 13.8 | 14.1 | Abstract / Sec 3.3, diagnosis 95% CI upper |
| 7.5 | 7.4 | Abstract / Sec 3.3, treatment 95% CI lower |
| 98.4 | 98.1 | Table 4, optimistic > diagnosis (%) |
| 80.1 | 77.4 | Table 4, realistic > diagnosis (%) |
| 74.7 | 71.4 | Table 4, treatment > diagnosis (%) |
| 73.6 | 73.9 | Table 4, realistic > treatment (%) |
| 4.2, 1.4, 1.7, 6.1 | 2.1, 1.1, 2.0, 6.2 | Sec 3.1 / Fig 2A, lower CI of mean months, stages 0, B, C, D |
| 32.6, 4.9, 4.5, 7.8 | 33.0, 4.8, 4.8, 7.7 | Sec 3.1 / Fig 2A, upper CI of mean months, stages 0, B, C, D |
| 73 | 70 | Sec 3.1, % of bootstraps with stage A longer than stage 0 |
| 5.7 | 6.0 | Sec 3.4, SE of Baseline |
| -0.87 | 0.29 | Sec 3.4, partial correlation Baseline ~ logit(p_0) |
| 0.92 | -0.02 | Sec 3.4, partial correlation Baseline ~ logit(Sens_0) |
| 5.7 | 5.6 | Sec 3.2, +/-25% L_A upper (years) |

**Ledger diff.** `output/manuscript_values.csv` before and after, 83 of 85
compared names changed (relative change above 1e-12), sorted by relative
change. Scenario point estimates move by at most 3e-4 relative; most of the
large changes are interval limits.

| Name | Old | New | Rel. change |
|---|---|---|---|
| pcor_logit_0 | -0.87291738 | 0.28766045 | 1.33 |
| pcor_logit_US_0 | 0.91955120 | -0.01781958 | 1.02 |
| mean_0_Lower | 0.34710958 | 0.17276499 | 0.502 |
| diag_scen_lower_high_risk | 6.46512578 | 4.60027683 | 0.288 |
| diag_scen_lower | 3.57197107 | 2.54556846 | 0.287 |
| diag_scen_lower_low_risk | 1.87886156 | 1.33926838 | 0.287 |
| mean_B_Lower | 0.11304834 | 0.09387211 | 0.17 |
| pcor_logit_US_A | 0.22093606 | 0.18669003 | 0.155 |
| mean_C_Lower | 0.14513914 | 0.16282434 | 0.122 |
| mean_C_Upper | 0.37236999 | 0.40295678 | 0.0821 |
| treat_vs_diag | 0.74660000 | 0.71350000 | 0.0443 |
| sd_noint | 5.71575241 | 5.95422621 | 0.0417 |
| muA_gt_mu0 | 0.72990000 | 0.70460000 | 0.0347 |
| real_vs_diag | 0.80140000 | 0.77450000 | 0.0336 |
| diag_scen_upper_high_risk | 25.23104702 | 25.78137224 | 0.0218 |
| diag_scen_upper | 13.82891637 | 14.12452272 | 0.0214 |
| diag_scen_upper_low_risk | 7.24354013 | 7.39733291 | 0.0212 |
| mean_D_Upper | 0.64924802 | 0.63759307 | 0.018 |
| mean_B_Upper | 0.40969494 | 0.40260773 | 0.0173 |
| mean_D_Lower | 0.51164257 | 0.51865707 | 0.0137 |
| mean_0_Upper | 2.71752357 | 2.75015753 | 0.012 |
| real_adh_upper_low_risk | 6.06775164 | 6.12347677 | 0.00918 |
| opt_adh_lower_low_risk | 7.28577850 | 7.22236758 | 0.0087 |
| real_adh_upper | 11.55696285 | 11.65706081 | 0.00866 |
| opt_adh_lower_high_risk | 25.45082331 | 25.23385620 | 0.00852 |
| real_adh_lower_low_risk | 4.14581196 | 4.11123807 | 0.00834 |
| opt_adh_upper_low_risk | 10.67418331 | 10.76087027 | 0.00812 |
| treat_scens_lower_low_risk | 3.88028671 | 3.84943204 | 0.00795 |
| real_adh_lower | 7.92163520 | 7.85929528 | 0.00787 |
| treat_scens_upper_low_risk | 5.35910795 | 5.40079646 | 0.00778 |
| real_adh_upper_high_risk | 20.99218777 | 21.15517411 | 0.00776 |
| treat_scens_lower | 7.46717389 | 7.40949097 | 0.00772 |
| opt_adh_lower | 13.91739613 | 13.81164934 | 0.0076 |
| treat_scens_lower_high_risk | 13.84559121 | 13.74106451 | 0.00755 |
| real_adh_lower_high_risk | 14.47795837 | 14.37043309 | 0.00743 |
| opt_adh_upper | 20.34234559 | 20.49309009 | 0.00741 |
| treat_scens_upper | 10.26094268 | 10.33120968 | 0.00685 |
| opt_adh_upper_high_risk | 36.95714819 | 37.20856285 | 0.0068 |
| treat_scens_upper_high_risk | 18.82840479 | 18.94920921 | 0.00642 |
| real_vs_treat | 0.73570000 | 0.73900000 | 0.00449 |
| opt_vs_diag | 0.98390000 | 0.98120000 | 0.00274 |
| L_B_Upper | 0.95073591 | 0.95282937 | 0.0022 |
| L_B_Lower | -0.95073591 | -0.95282937 | 0.0022 |
| Sens_B_Lower | -0.26196522 | -0.26247325 | 0.00194 |
| Sens_B_Upper | 0.19706595 | 0.19743744 | 0.00189 |
| mean_noint | 67.56751588 | 67.45410153 | 0.00168 |
| mean_B | 0.24257820 | 0.24290821 | 0.00136 |
| mean_0 | 1.04401400 | 1.04261691 | 0.00134 |
| Sens_C_Upper | 0.00525347 | 0.00525995 | 0.00123 |
| Sens_C_Lower | -0.01183993 | -0.01185451 | 0.00123 |
| L_A_Lower | -5.65336840 | -5.64673525 | 0.00117 |
| L_A_Upper | 5.65336840 | 5.64673525 | 0.00117 |
| Sens_A_Lower | -2.07463683 | -2.07274445 | 0.000912 |
| diag_scen_high_risk | 12.68456378 | 12.67350088 | 0.000872 |
| diag_scen | 6.99067098 | 6.98459887 | 0.000869 |
| diag_scen_low_risk | 3.67232174 | 3.66913885 | 0.000867 |
| Sens_0_Lower | -7.19653625 | -7.20218646 | 0.000785 |
| Sens_0_Upper | 6.57151249 | 6.57651154 | 0.000761 |
| mean_A | 1.50404971 | 1.50296832 | 0.000719 |
| L_0_Upper | 10.20174165 | 10.20894754 | 0.000706 |
| L_0_Lower | -10.20174165 | -10.20894754 | 0.000706 |
| mean_D | 0.51933340 | 0.51900640 | 0.00063 |
| Sens_A_Upper | 1.70860369 | 1.70759455 | 0.000591 |
| mean_A_Lower | 1.02005129 | 1.01951988 | 0.000521 |
| p_CA_Upper | 3.63687539 | 3.63866845 | 0.000493 |
| p_CA_Lower | -3.63687539 | -3.63866845 | 0.000493 |
| treat_scens_high_risk | 16.58339269 | 16.57533742 | 0.000486 |
| treat_scens | 8.99979444 | 8.99545591 | 0.000482 |
| treat_scens_low_risk | 4.69064603 | 4.68839341 | 0.00048 |
| p_Att_Upper | 3.53955114 | 3.54086483 | 0.000371 |
| L_C_Lower | -0.19997014 | -0.20004124 | 0.000356 |
| L_C_Upper | 0.19997014 | 0.20004124 | 0.000356 |
| mean_A_Upper | 1.95754953 | 1.95817244 | 0.000318 |
| real_adh_low_risk | 5.12504061 | 5.12656581 | 0.000298 |
| real_adh | 9.77259767 | 9.77547110 | 0.000294 |
| real_adh_high_risk | 17.79080223 | 17.79591092 | 0.000287 |
| opt_adh_low_risk | 9.01248131 | 9.01503881 | 0.000284 |
| opt_adh | 17.18642613 | 17.19123928 | 0.00028 |
| opt_adh_high_risk | 31.29154449 | 31.30008340 | 0.000273 |
| mean_C | 0.23675822 | 0.23669528 | 0.000266 |
| p_Att_Lower | -3.94904525 | -3.95004503 | 0.000253 |
| p_NA_Upper | -4.43596458 | -4.43563197 | 7.5e-05 |
| p_NA_Lower | 4.43596458 | 4.43563197 | 7.5e-05 |

**Tests and fixtures regenerated under G8 + G9:**
`tests/testthat/fixtures/golden.rds` (via `make_golden.R`; only the entries that
depend on the calibrated likelihood or parameters changed, and the `param_names`
names are identical); the plan anchors in `test-golden.R`; the published-value
anchors in `test-paper.R` (table above); and a NEW fixture
`tests/testthat/fixtures/ledger_baseline.rds` (generator
`fixtures/make_ledger_baseline.R`), which `test-ledger.R` now compares the 02 and
04 ledger values against (tolerance still 1e-12). It holds the 85 ledger values
whose names are in the frozen reference. `data-derived/param_names.rds` is NOT
modified: it remains the frozen record of the published values, and a test now
checks that it is unchanged (91 values; published `diag_scen_lower` 3.57197...).
Both fixtures are `*.rds`, so they are stored through Git LFS like `golden.rds`.
`test-likelihood.R` has new tests for the exact absorption probabilities,
`fit_mle`, `fit_at_bound`, polish-only starts and `refit_bootstrap`
determinism (cores 1 vs 2). The old G7 test ("Nelder-Mead from the documented
start reproduces the CSV") was REPLACED, not kept: it is now "fit_mle from
c(NM_START, CSV) reproduces the CSV". No tolerance was loosened.

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
