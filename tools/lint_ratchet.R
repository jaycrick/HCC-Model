#!/usr/bin/env Rscript
# Lint ratchet: per-file lint counts may never exceed tools/lint_baseline.csv.
#
# Usage:
#   Rscript tools/lint_ratchet.R [files...]   check (default: all project .R)
#   Rscript tools/lint_ratchet.R --update     lower the baseline (never raise)
#
# Exit status is non-zero if any file has more lints than its baseline, or a
# file absent from the baseline has any lint.
library(lintr)

baseline_path <- "tools/lint_baseline.csv"

args <- commandArgs(trailingOnly = TRUE)
update <- "--update" %in% args
files <- setdiff(args, "--update")
explicit <- length(files) > 0

if (!explicit) {
  files <- system2(
    "git",
    c("ls-files", "-co", "--exclude-standard", "--", shQuote("*.R")),
    stdout = TRUE
  )
}
files <- unique(sub("^\\./", "", files))
files <- files[grepl("\\.R$", files)]
files <- files[file.exists(files)]

read_baseline <- function() {
  if (!file.exists(baseline_path)) {
    return(setNames(integer(0), character(0)))
  }
  b <- utils::read.csv(baseline_path, stringsAsFactors = FALSE)
  setNames(as.integer(b$count), b$file)
}
initial <- !file.exists(baseline_path)
baseline <- read_baseline()

lints <- lapply(files, lintr::lint)
names(lints) <- files
counts <- vapply(lints, length, integer(1))

if (update) {
  new <- baseline
  # Without explicit files, drop entries for files that no longer exist.
  if (!explicit) {
    new <- new[names(new) %in% files]
  }
  raised <- character(0)
  for (f in files) {
    if (f %in% names(baseline)) {
      if (counts[[f]] > baseline[[f]]) {
        raised <- c(raised, sprintf("%s: %d > %d", f, counts[[f]], baseline[[f]]))
      } else {
        new[[f]] <- counts[[f]]
      }
    } else if (counts[[f]] > 0 && !initial) {
      raised <- c(raised, sprintf("%s: new file with %d lint(s)", f, counts[[f]]))
    } else {
      # Zero-count entries may be added; the very first run seeds any count.
      new[[f]] <- counts[[f]]
    }
  }
  if (length(raised) > 0) {
    message("REFUSING to raise the baseline. Fix these lints instead:")
    message(paste0("  ", raised, collapse = "\n"))
    quit(status = 1)
  }
  new <- new[order(names(new))]
  utils::write.csv(
    data.frame(file = names(new), count = as.integer(new)),
    baseline_path,
    row.names = FALSE,
    quote = FALSE
  )
  cat(sprintf("Baseline written: %d files, %d lints total.\n", length(new), sum(new)))
  quit(status = 0)
}

failed <- FALSE
for (f in files) {
  allowed <- if (f %in% names(baseline)) baseline[[f]] else 0L
  if (counts[[f]] > allowed) {
    failed <- TRUE
    cat(sprintf("\n%s: %d lint(s), baseline allows %d\n", f, counts[[f]], allowed))
    print(lints[[f]])
  }
}

if (failed) {
  cat(
    "\nLint ratchet FAILED: lint counts may not increase.",
    "Fix the lints above (see .lintr).\n"
  )
  quit(status = 1)
}
cat(sprintf(
  "Lint ratchet OK: %d file(s), %d lint(s) (baseline %d).\n",
  length(files),
  sum(counts),
  sum(baseline[names(baseline) %in% files])
))
