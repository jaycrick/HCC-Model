#!/usr/bin/env Rscript
# One-off calibration refit (gates G8 + G9; plan EFF-2).
#
# Refits the point MLE and all bootstrap resamples with fit_mle() (exact
# absorption probabilities, Nelder-Mead + bounded logit L-BFGS-B polish),
# then regenerates the calibrated rows of the parameter CSVs.
#
# Usage (from the project root):
#   Rscript tools/refit_calibration.R [--n <int>] [--cores <int>]
#                                     [--config a|b|c|d]
#                                     [--audit <int>] [--dry-run] [--out <dir>]
#
#   --n        refit only the first n resamples (default: all 10 000)
#   --cores    forked workers (default: detectCores() - 1); results do not
#              depend on it
#   --config   start strategy per resample (default d):
#                a  NM + polish from NM_START, the new point MLE and the old
#                   HPC estimate;
#                b  NM + polish from the new point MLE and the old HPC estimate;
#                c  NM + polish from NM_START; polish only from the new point
#                   MLE and from the old HPC estimate;
#                d  polish only from the old HPC estimate; if that fit does not
#                   converge, sits on the box or the polish fails, rerun the
#                   resample with config b and keep the better result
#                   (diag$fallback; best_start codes the start, see
#                   fit_resample() in R/bootstrap.R).
#              The old HPC estimate is clipped into the box (FIT_LB, FIT_UB).
#   --audit    run config a on this many evenly spaced resamples after the main
#              refit and compare NLLs (default 200; 0 disables); a gap above
#              1e-6 prints a WARNING but does not stop the script
#   --dry-run  write NOTHING to data/ or data-derived/; outputs go to --out
#   --out      output directory for a dry run (default: tempdir())
#
# Real run writes:
#   data-derived/resampled_pars.rds              (n x 8, NLL_PAR_NAMES)
#   data-derived/bootstrap_fit_diagnostics.rds   (diag + old HPC NLL under the
#                                                 new likelihood, config, time)
#   data/Parameters.csv                          calibrated rows: POINT = MLE,
#                                                LOWER/UPPER = bootstrap 2.5 /
#                                                97.5 % (summarise_ci()),
#                                                Sym_B = 8 x Sym_A
#   data/Parameters_25.csv                       calibrated rows: POINT and
#                                                POINT * (1 -/+ 0.25)
# The refit runs in chunks of 500 resamples; each chunk is checkpointed to
# refit_checkpoint/chunk_<k>.rds (in data-derived/, or --out under --dry-run)
# and skipped on a rerun. The directory is deleted after the final write.
# Values are written with 8 significant digits; every other line, the BOM and
# the line endings are kept byte for byte.

suppressPackageStartupMessages({
  library(hypogeom)
  library(DescTools)
})
source(here::here("R", "setup.R"))
source(here::here("R", "model_core.R"))
source(here::here("R", "bootstrap.R"))
source(here::here("R", "plot_theme.R"))

# ----- Command line -----------------------------------------------------------

parse_args <- function(args) {
  get <- function(flag, default) {
    i <- match(flag, args)
    if (is.na(i)) default else args[[i + 1]]
  }
  list(
    n = as.integer(get("--n", NA)),
    cores = as.integer(get(
      "--cores",
      max(1L, parallel::detectCores() - 1L)
    )),
    config = get("--config", "d"),
    audit = as.integer(get("--audit", 200L)),
    dry_run = "--dry-run" %in% args,
    out = get("--out", tempdir())
  )
}

# ----- Start strategies -------------------------------------------------------

clip_box <- function(par) pmin(pmax(par, FIT_LB), FIT_UB)

make_starts_fn <- function(config, point_par, hpc) {
  hpc_b <- function(b) clip_box(unlist(hpc[b, NLL_PAR_NAMES]))
  switch(
    config,
    a = function(b) list(NM_START, point_par, hpc_b(b)),
    b = function(b) list(point_par, hpc_b(b)),
    c = function(b) {
      list(
        NM_START,
        list(par = point_par, nm = FALSE),
        list(par = hpc_b(b), nm = FALSE)
      )
    },
    d = function(b) list(list(par = hpc_b(b), nm = FALSE)),
    stop("--config must be a, b, c or d")
  )
}

# Fallback starts: config d falls back to config b (NM + polish from the point
# MLE and from the clipped HPC estimate).
make_fallback_fn <- function(config, point_par, hpc) {
  if (config == "d") make_starts_fn("b", point_par, hpc) else NULL
}

# NLL of the old HPC estimates under the NEW likelihood, per resample.
old_hpc_nll <- function(inputs, hpc, cores) {
  vals <- parallel::mclapply(
    seq_along(inputs),
    function(b) {
      x <- inputs[[b]]
      neg_log_lik(
        unlist(hpc[b, NLL_PAR_NAMES]),
        meds = x[1:5],
        Dx_0 = x[6],
        Dx_A = x[7],
        Dx_sym_early = x[8]
      )
    },
    mc.cores = cores
  )
  unlist(vals)
}

# ----- CSV rows ---------------------------------------------------------------

CAL_ROWS <- c(NLL_PAR_NAMES, "Sym_B")

fmt8 <- function(x) as.character(signif(x, 8))

# Replace the calibrated lines of a parameter CSV; all other bytes unchanged.
replace_rows <- function(path, new) {
  txt <- readLines(path, warn = FALSE) # BOM is dropped; see write_lines_keep()
  nm <- sub(",.*$", "", txt)
  for (r in rownames(new)) {
    i <- match(r, nm)
    stopifnot(!is.na(i))
    txt[i] <- paste(c(r, fmt8(new[r, ])), collapse = ",")
  }
  txt
}

write_lines_keep <- function(txt, path, src) {
  # Same terminator as the shipped files: LF after every line. readLines()
  # drops the UTF-8 BOM, so restore it if the source file had one.
  bom <- identical(readBin(src, "raw", 3), as.raw(c(0xef, 0xbb, 0xbf)))
  out <- charToRaw(paste0(paste(txt, collapse = "\n"), "\n"))
  writeBin(c(if (bom) as.raw(c(0xef, 0xbb, 0xbf)), out), path)
}

main_rows <- function(point_par, pars) {
  ci <- vapply(pars[NLL_PAR_NAMES], summarise_ci, numeric(2))
  m <- rbind(POINT = point_par, ci)
  m <- t(m)
  colnames(m) <- c("POINT", "LOWER", "UPPER")
  m <- rbind(m, Sym_B = SYM_B_RATIO * m["Sym_A", ])
  m[CAL_ROWS, , drop = FALSE]
}

rows_25 <- function(point_par) {
  p <- c(point_par, Sym_B = SYM_B_RATIO * point_par[["Sym_A"]])[CAL_ROWS]
  cbind(POINT = p, LOWER = 0.75 * p, UPPER = 1.25 * p)
}

before_after <- function(old, new) {
  cmp <- cbind(
    old_POINT = old[CAL_ROWS, "POINT"],
    new_POINT = new[CAL_ROWS, "POINT"],
    old_LOWER = old[CAL_ROWS, "LOWER"],
    new_LOWER = new[CAL_ROWS, "LOWER"],
    old_UPPER = old[CAL_ROWS, "UPPER"],
    new_UPPER = new[CAL_ROWS, "UPPER"]
  )
  print(signif(cmp, 6))
}

# ----- Main -------------------------------------------------------------------

CHUNK <- 500L

# Refit all `ids` in chunks, with checkpoints, progress and ETA.
refit_chunked <- function(inputs, starts_fn, fallback_fn, a, ck_dir) {
  ids <- seq_along(inputs)
  chunks <- split(ids, ceiling(ids / CHUNK))
  dir.create(ck_dir, showWarnings = FALSE, recursive = TRUE)
  t0 <- Sys.time()
  st <- new.env()
  st$done <- 0
  out <- lapply(seq_along(chunks), function(k) {
    f <- file.path(ck_dir, paste0("chunk_", k, ".rds"))
    r <- if (file.exists(f)) readRDS(f)
    if (!is.null(r) && identical(r$key, c(a$config, length(chunks[[k]])))) {
      message("chunk ", k, "/", length(chunks), ": resumed from checkpoint")
      return(r$res)
    }
    res <- refit_bootstrap(
      inputs,
      starts_fn,
      cores = a$cores,
      fallback_fn = fallback_fn,
      ids = chunks[[k]]
    )
    saveRDS(list(key = c(a$config, length(chunks[[k]])), res = res), f)
    st$done <- st$done + length(chunks[[k]])
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    message(
      "chunk ",
      k,
      "/",
      length(chunks),
      " done; elapsed ",
      round(el / 60, 1),
      " min; ETA ",
      round(el / st$done * (length(ids) - st$done) / 60, 1),
      " min"
    )
    res
  })
  list(
    pars = dplyr::bind_rows(lapply(out, `[[`, "pars")),
    diag = dplyr::bind_rows(lapply(out, `[[`, "diag"))
  )
}

# Config a on an evenly spaced subset; compare with the main result.
run_audit <- function(inputs, hpc, point_par, res, a) {
  n <- length(inputs)
  ids <- unique(round(seq(1, n, length.out = min(a$audit, n))))
  au <- refit_bootstrap(
    inputs,
    make_starts_fn("a", point_par, hpc),
    cores = a$cores,
    ids = ids
  )$diag
  out <- tibble::tibble(
    b = ids,
    audit_value = au$value,
    main_value = res$diag$value[ids],
    gap = main_value - audit_value
  )
  nbad <- sum(out$gap > 1e-6)
  message(
    "Audit of ",
    length(ids),
    " resamples: max gap ",
    signif(max(out$gap), 3),
    ", ",
    nbad,
    " with gap > 1e-6"
  )
  if (nbad > 0) {
    message(
      "WARNING: *** ",
      nbad,
      " audited resample(s) have a better optimum ",
      "under config a than the main refit (gap > 1e-6): b = ",
      paste(out$b[out$gap > 1e-6], collapse = ", "),
      " ***"
    )
  }
  out
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  a <- parse_args(args)
  t0 <- Sys.time()
  inputs <- readRDS(paper_path(
    "data-derived",
    "list_of_resampled_milton_output.rds"
  ))
  hpc <- bootstrap_estimates()
  n <- if (is.na(a$n)) length(inputs) else a$n
  inputs <- inputs[seq_len(n)]

  message("Point fit ...")
  # Fixed starts (not the CSV, which this script rewrites): see PUBLISHED_MLE.
  point <- fit_mle(list(NM_START, PUBLISHED_MLE))
  print(point[c("par", "value", "convergence", "at_bound")])
  message(
    "Refit of ",
    n,
    " resamples, config ",
    a$config,
    ", ",
    a$cores,
    " cores ..."
  )
  ck_dir <- file.path(
    if (a$dry_run) a$out else paper_path("data-derived"),
    "refit_checkpoint"
  )
  tm <- system.time(
    res <- refit_chunked(
      inputs,
      make_starts_fn(a$config, point$par, hpc),
      make_fallback_fn(a$config, point$par, hpc),
      a,
      ck_dir
    )
  )
  message(
    "Elapsed ",
    round(tm[["elapsed"]], 1),
    " s (",
    round(tm[["elapsed"]] / n, 2),
    " s per resample)"
  )
  diag <- res$diag |>
    dplyr::mutate(
      old_hpc_value = old_hpc_nll(inputs, hpc, a$cores),
      config = a$config,
      timestamp = format(t0, "%Y-%m-%d %H:%M:%S")
    )
  audit <- if (a$audit > 0) run_audit(inputs, hpc, point$par, res, a)
  write_outputs(a, point, res$pars, diag, audit)
  unlink(ck_dir, recursive = TRUE)
  invisible(list(point = point, res = res, diag = diag, audit = audit))
}

write_outputs <- function(a, point, pars, diag, audit = NULL) {
  dd <- if (a$dry_run) a$out else paper_path("data-derived")
  dt <- if (a$dry_run) a$out else paper_path("data")
  dir.create(a$out, showWarnings = FALSE, recursive = TRUE)
  saveRDS(pars, file.path(dd, "resampled_pars.rds"))
  saveRDS(
    list(diag = diag, audit = audit),
    file.path(dd, "bootstrap_fit_diagnostics.rds")
  )
  saveRDS(point, file.path(a$out, "point_fit.rds"))

  main_path <- paper_path("data", "Parameters.csv")
  new <- main_rows(point$par, pars)
  before_after(read.csv(main_path, row.names = 1), new)
  path_25 <- paper_path("data", "Parameters_25.csv")
  write_lines_keep(
    replace_rows(main_path, new),
    file.path(dt, "Parameters.csv"),
    main_path
  )
  write_lines_keep(
    replace_rows(path_25, rows_25(point$par)),
    file.path(dt, "Parameters_25.csv"),
    path_25
  )
  message("Outputs written to ", dd, " and ", dt)
}

if (sys.nframe() == 0) {
  main()
}
