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
and every sweep (`sweep_model()` in R/plotting.R) use it. The point-estimate
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
