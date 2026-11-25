# Paper 2: stochastic model with additional compartments (in preparation)

Working folder for the second paper. It uses the stochastic version of the model
(`R/simulation.R`) and extends the model structure with additional compartments.
Rename the folder (for example `papers/2027-<journal>/`) once the target is known,
and update the tables in the root README.

Status: the adherence model is fitted (`analysis/01_fit_adherence.R`); the
compartment extension and the scenario analyses are still to come.

Layout, following [`papers/2024-cancers`](../2024-cancers):

```
papers/stochastic/
├── README.md
├── run_all.R          # options(hcc.paper = "stochastic"), then the analyses in order
├── analysis/
│   └── 01_fit_adherence.R   # adherence fit, CIs, bootstrap, GOF, sensitivity
├── data/              # this paper's Parameters CSVs and the aggregate adherence tables
│   └── README.md      # how the tables were cleaned
├── data-derived/      # adherence_bootstrap.rds (regenerable, Git LFS)
├── output/            # git-ignored (ledger, model comparison, sensitivity)
└── figures/           # git-ignored
```

Run from the repository root: `Rscript papers/stochastic/run_all.R`
(`--bootstrap` reruns the bootstrap instead of reusing the stored draws).

## Adherence model ("due until attend")

Paper 1 cut the observed attendance rate into three groups, and each group
attended every scheduled 6-monthly scan with a fixed probability. Paper 2 uses
a latent-class model instead (`adherence = "due"` in `R/model_core.R` and
`R/simulation.R`):

- each person belongs to a latent class Low, Med or High with shares
  `Pi_Low`, `Pi_Med`, `Pi_High` (sum 1);
- a person starts **due**; each month while due they attend with probability
  `Att_<class>`; after an attended scan they are not due for m = 6 months, then
  due again.

The likelihood depends on each person only through A (scans attended while due)
and F (months due without a scan):
l = sum n log sum_c pi_c q_c^A (1 - q_c)^F. It is fitted by EM (R/adherence.R),
with classes ordered Att_Low < Att_Med < Att_High. **The chosen model has 3
classes with Att_Low fixed at 0** (see "Chosen model" below). Intervals in
`data/Parameters.csv` are nonparametric bootstrap percentile 95% CIs (people
resampled, B = 1000, each resample refitted with Att_Low = 0). Profile-likelihood
CIs for the four free parameters (Pi_Low, Pi_High, Att_Med, Att_High; computed
under the ordering Att_Low <= Att_Med <= Att_High so classes cannot swap labels)
are recorded in the ledger. Only 2 of 1000 resamples have a free parameter at a
bound (it was 330 of 1000 under the free model).

**Base cohort.** The base analysis follows people to 2020-03-01 at most (as in
paper 1), N = 108. Follow-up to the 2022-03-18 census is the sensitivity
variant `withcovid`. The tables carry the exposure T (months of follow-up after
the shift rule) next to A and F; the likelihood ignores T, which is used for
simulation (goodness of fit, parametric bootstrap) with each person's exact
follow-up.

### Model selection (number of classes)

K = 1 to 5 classes, with and without the lowest class's q fixed at 0.
df = 2K - 1 (K - 1 shares + K q's), or 2K - 2 with the zero class.
For every criterion n = number of **people** (108), never person-months or
cells, because people are the independent units. BIC is the primary criterion
for the number of classes (Keribin 2000); AICc and ICL are shown alongside, and
AIC alone overfits mixtures. Singular BIC (Drton & Plummer 2017) is
theoretically better suited to mixtures but is not used: its learning
coefficients for this geometric-type mixture are not readily available.

Base (pre-COVID) table; `0` = lowest class q fixed at 0; `deg` = degenerate fit
(a share < 1e-4, two q's within 1e-4, or a q on the boundary):

| K | zero | loglik | df | AIC | AICc | BIC | ICL | deg |
|---|---|---|---|---|---|---|---|---|
| 1 |  | -968.7 | 1 | 1939.5 | 1939.5 | 1942.2 | 1942.2 |  |
| 2 |  | -805.8 | 3 | 1617.7 | 1617.9 | 1625.7 | 1637.9 |  |
| 2 | 0 | -846.3 | 2 | 1696.5 | 1696.7 | 1701.9 | 1703.2 |  |
| 3 |  | -796.6 | 5 | 1603.3 | 1603.9 | 1616.7 | 1664.3 |  |
| 3 | 0 | -797.2 | 4 | 1602.5 | 1602.8 | 1613.2 | 1655.0 |  |
| 4 |  | -795.9 | 7 | 1605.8 | 1606.9 | 1624.6 | 1700.9 |  |
| 4 | 0 | -796.3 | 6 | 1604.6 | 1605.4 | 1620.6 | 1694.1 |  |
| 5 |  | -795.9 | 9 | 1609.8 | 1611.7 | 1634.0 | 1743.4 | deg |
| 5 | 0 | -795.9 | 8 | 1607.8 | 1609.3 | 1629.3 | 1717.3 |  |

With-COVID table:

| K | zero | loglik | df | AIC | AICc | BIC | ICL | deg |
|---|---|---|---|---|---|---|---|---|
| 1 |  | -1342.4 | 1 | 2686.8 | 2686.8 | 2689.5 | 2689.5 |  |
| 2 |  | -1108.9 | 3 | 2223.9 | 2224.1 | 2231.9 | 2239.4 |  |
| 2 | 0 | -1175.4 | 2 | 2354.9 | 2355.0 | 2360.2 | 2361.7 |  |
| 3 |  | -1089.5 | 5 | 2188.9 | 2189.5 | 2202.3 | 2236.5 |  |
| 3 | 0 | -1095.4 | 4 | 2198.8 | 2199.2 | 2209.5 | 2241.2 |  |
| 4 |  | -1081.6 | 7 | 2177.2 | 2178.3 | 2196.0 | 2253.4 |  |
| 4 | 0 | -1081.6 | 6 | 2175.2 | 2176.1 | 2191.3 | 2245.2 |  |
| 5 |  | -1081.1 | 9 | 2180.3 | 2182.1 | 2204.4 | 2296.4 | deg |
| 5 | 0 | -1081.1 | 8 | 2178.3 | 2179.7 | 2199.7 | 2291.7 |  |

On the base table BIC is lowest for K = 3 with Att_Low = 0 (BIC 1613.2),
3.5 below the free 3-class model (1616.7); K = 3 is clearly supported over K = 1
and K = 2 (BIC differences 329 and 12.6), and K = 4 and K = 5 are worse on BIC,
AICc and ICL (the free K = 5 fit is degenerate).
AIC, AICc and ICL also prefer Att_Low = 0 (AICc 1602.8 vs 1603.9; ICL 1655.0 vs 1664.3).
So **under the pre-COVID base Att_Low = 0 is preferred and is not rejected**:
LR = 1.178, boundary p = 0.139 (0.5 * P(chi2_1 > LR), Self & Liang 1987; valid
here because the other two classes stay identifiable under the null), parametric
bootstrap LRT p = 0.095 (B = 200, simulated from the Att_Low = 0 fit with each
person's exact T). The two p-values disagree mildly; neither rejects at 5%.
Only one person has A >= 1 and a posterior probability above 0.5 of the low class
under the free fit (A = 1, F = 72); occasional scans by such people drive
Att_Low > 0. On the with-COVID table the free 3-class model is preferred to
Att_Low = 0 (BIC 2202.3 vs 2209.5), but K = 4 with a zero class has the lowest BIC
(2191.3), so with COVID-era follow-up K = 3 is not the BIC choice.

### Chosen model: 3 classes, Att_Low fixed at 0

`Parameters.csv` carries the 3-class model with Att_Low = 0 (POINT = LOWER =
UPPER = 0). The free 3-class model is a sensitivity analysis. Reasons:

- BIC, AICc and ICL all prefer it (table above);
- the LRTs do not reject Att_Low = 0 (boundary p = 0.14, parametric bootstrap
  p = 0.095);
- under the free model only 1 person with A >= 1 (A = 1, F = 72) has a
  posterior probability above 0.5 of the low class, so Att_Low > 0 rests on one
  occasional scan;
- the low class (Pi_Low = 0.27) matches paper 1's non-adherent group
  (p_NA = 0.31, attendance rate < 0.1, "never attend");
- the free fit's Att_Low is 0.0025 and its interval reaches the bound, and
  resamples hit a bound in 330 of 1000 draws versus 2 of 1000 here.

Caveat: scans at other hospitals are not observed, so "never attends" means
"never attends at this hospital"; those people may be scanned elsewhere.

| Parameter | Estimate | Bootstrap 95% CI | Profile 95% CI |
|---|---|---|---|
| Pi_Low | 0.266 | 0.148 to 0.362 | 0.165 to 0.366 |
| Pi_Med | 0.275 | 0.181 to 0.396 | (derived) |
| Pi_High | 0.460 | 0.361 to 0.568 | 0.357 to 0.563 |
| Att_Low | 0 (fixed) | 0 to 0 | |
| Att_Med | 0.0516 | 0.0244 to 0.0796 | 0.0311 to 0.0764 |
| Att_High | 0.299 | 0.260 to 0.347 | 0.263 to 0.339 |

Sensitivity (all Att_Low = 0 except the last row; shares and monthly
attendance probabilities):

| Variant | Pi_Low | Pi_Med | Pi_High | Att_Med | Att_High | loglik |
|---|---|---|---|---|---|---|
| base | 0.266 | 0.275 | 0.460 | 0.0516 | 0.299 | -797.2 |
| withcovid | 0.278 | 0.265 | 0.457 | 0.0464 | 0.268 | -1095.4 |
| diagscan | 0.266 | 0.275 | 0.459 | 0.0518 | 0.301 | -798.4 |
| noshift | 0.235 | 0.271 | 0.495 | 0.0354 | 0.230 | -786.5 |
| floor | 0.263 | 0.280 | 0.457 | 0.0511 | 0.311 | -785.5 |
| free Att_Low | 0.318 | 0.236 | 0.446 | 0.0668 | 0.304 | -796.6 |

(free Att_Low also has Att_Low = 0.0025.) Goodness of fit (plug-in posterior
predictive, 500 replicates): p = 0.95 (A bins), p = 0.84 (rate bins).

### Parameter names

Paper 2 `Parameters.csv` has **no** `p_NA`, `p_CA` or `p_Att`, so the two
adherence models cannot be mixed. Scripts pass `adherence = "due"` explicitly.

| Name | Used by | Meaning |
|------|---------|---------|
| `p_NA` | `"scheduled"` (paper 1) | Share with attendance rate < 0.1 (never attend) |
| `p_CA` | `"scheduled"` | Share with attendance rate >= 0.8 (attend every scan) |
| `p_Att` | `"scheduled"` | Probability a middle-group person attends a scheduled scan |
| `Pi_Low`, `Pi_Med`, `Pi_High` | `"due"` (paper 2) | Latent-class shares (sum 1) |
| `Att_Low`, `Att_Med`, `Att_High` | `"due"` | Monthly probability of attending while due |

The values are therefore not comparable between papers: `Att_*` are monthly
rates while due, not probabilities per scheduled scan.

### Data provenance

The model is fitted to unpublished tertiary-hospital attendance data. The raw
person-level data are held by the authors and are **not public**. Only the
aggregate tables `data/adherence_counts*.csv` (columns A, F, T, n; no dates or
identifiers) are committed, and `01_fit_adherence.R` reads only those, so the
analysis runs from a fresh clone. How the tables are built is documented in
[`data/README.md`](data/README.md). The authors regenerate them with a private
preparation script (kept outside the repository) that reads the raw folder
named by the environment variable `HCC_ADH_DATA_DIR`, using the file manifest in
that folder, and reuses the cleaning functions in `R/adherence_data.R`.

Ground rules (see "Adding a paper" in the root README):

- Change the model in the shared core (`R/`), not in copies here.
- Paper 1 must stay a special case of the extended model: the test suite checks
  its published values on every push.
