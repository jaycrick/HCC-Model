# ONE-OFF generator for tests/testthat/fixtures/likelihood_golden.rds.
#
# RULE: run ONCE against hypogeom damian-t-p@2b329c4 (before the bump to the
# exact-prefix fork). NEVER regenerate to make a failing test pass. Regenerate
# only under an approved gate and record the gate in the commit message.
#
# Usage (from repo root):  Rscript tests/testthat/fixtures/make_likelihood_golden.R

source(here::here("tests", "testthat", "helper-setup.R"))

default_data <- list(
  meds = calibration_data$stage$med,
  Dx_0 = calibration_data$dx[["Dx_0"]],
  Dx_A = calibration_data$dx[["Dx_A"]],
  Dx_sym_early = calibration_data$dx[["Dx_sym_early"]]
)
point <- function(label, par, data = default_data) {
  list(label = label, par = par, data = data)
}

points <- list(point("csv point", csv_par))

# (b) 40 deterministic relative perturbations (up to +-2%, alternating sign)
for (i in 1:40) {
  d <- ((i - 20.5) / 20) * 0.02
  par <- csv_par * (1 + d * (-1)^(seq_len(8) + i))
  par <- pmin(pmax(par, FIT_LB), FIT_UB)
  points[[length(points) + 1]] <- point(paste("perturbation", i), par)
}

# (c) first 20 bootstrap resamples at their HPC estimates
inputs <- readRDS(paper_path("data-derived", "list_of_resampled_milton_output.rds"))[
  1:20
]
hpc <- bootstrap_estimates()
for (b in 1:20) {
  x <- inputs[[b]]
  data <- list(meds = x[1:5], Dx_0 = x[6], Dx_A = x[7], Dx_sym_early = x[8])
  par <- unlist(hpc[b, NLL_PAR_NAMES])
  points[[length(points) + 1]] <- point(paste("resample", b), par, data)
}

# (d) box corner, (e) near-1
corner <- csv_par
corner[STAGE_P_NAMES] <- FIT_LB[STAGE_P_NAMES]
points[[length(points) + 1]] <- point("box corner", corner)
near1 <- csv_par
near1[["p_B"]] <- 1 - 1e-6
near1[["p_C"]] <- 0.999
points[[length(points) + 1]] <- point("near 1", near1)

for (k in seq_along(points)) {
  pt <- points[[k]]
  points[[k]]$nll <- do.call(neg_log_lik, c(list(pt$par), pt$data))
  points[[k]]$ll_med <- loglik_medians(pt$par, pt$data$meds)
}

f <- fit_mle(list(NM_START))
pkg <- packageDescription("hypogeom")
# The installed DESCRIPTION has no RemoteSha; take the pinned SHA from the lock.
lock <- readLines(here::here("rproj.lock"))
sha_line <- grep("RemoteRepo = \"damian-t-p/hypogeom\"", lock, value = TRUE)[1]
lock_sha <- sub('.*RemoteSha = "([0-9a-f]+)".*', "\\1", sha_line)
golden <- list(
  points = points,
  fit = list(par = f$par, value = f$value),
  provenance = list(
    hypogeom_version = pkg$Version,
    hypogeom_sha = if (is.null(pkg$RemoteSha)) lock_sha else pkg$RemoteSha,
    r_version = R.version.string
  )
)

saveRDS(golden, here::here("tests", "testthat", "fixtures", "likelihood_golden.rds"))
message("wrote likelihood_golden.rds (", length(points), " points)")
