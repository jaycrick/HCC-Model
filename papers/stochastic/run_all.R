# Run all analyses for paper 2, stochastic model (file papers/stochastic/run_all.R)
# Run from the project root (where MODEL.Rproj lives):
#
#   Rscript papers/stochastic/run_all.R [--bootstrap]
#
# Uses only the committed aggregate table data/adherence_counts.csv (and its
# sensitivity variants); the raw data are not needed. The adherence bootstrap
# in data-derived/ is reused if present; --bootstrap reruns it.

options(hcc.paper = "stochastic")
library(here)

paper_dir <- here("papers", "stochastic")

message("==> Running 01_fit_adherence.R ...")
source(file.path(paper_dir, "analysis", "01_fit_adherence.R"))

message("==> Done. Figures written to ", file.path(paper_dir, "figures"))
