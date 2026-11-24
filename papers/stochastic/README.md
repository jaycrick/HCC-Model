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
with classes ordered Att_Low < Att_Med < Att_High. Intervals are
profile-likelihood 95% CIs (Pi_Med, being derived, uses bootstrap percentile
limits), checked against a nonparametric bootstrap of people (B = 1000).

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
aggregate tables `data/adherence_counts*.csv` (columns A, F, n; no dates or
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
