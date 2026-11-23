# Paper 2: stochastic model with additional compartments (in preparation)

Working folder for the second paper. It uses the stochastic version of the model
(`R/simulation.R`) and extends the model structure with additional compartments.
Rename the folder (for example `papers/2027-<journal>/`) once the target is known,
and update the tables in the root README.

Planned layout, following [`papers/2024-cancers`](../2024-cancers):

```
papers/stochastic/
├── README.md
├── run_all.R          # options(hcc.paper = "stochastic"), then the analyses in order
├── analysis/
├── data/              # this paper's Parameters CSVs
├── data-derived/
├── output/            # git-ignored
└── figures/           # git-ignored
```

Ground rules (see "Adding a paper" in the root README):

- Change the model in the shared core (`R/`), not in copies here.
- Paper 1 must stay a special case of the extended model: the test suite checks
  its published values on every push.
