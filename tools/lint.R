#!/usr/bin/env Rscript
# Lint gate: every project .R file must have zero lints (rules in .lintr).
#
# Usage:
#   Rscript tools/lint.R [files...]   check (default: all project .R files)
#
# Exit status is non-zero if any file has a lint. To exempt a line, use a
# `# nolint` comment; to exempt a whole file, add it to `exclusions:` in .lintr.
library(lintr)

files <- commandArgs(trailingOnly = TRUE)
if (length(files) == 0) {
  # Tracked and untracked-but-not-ignored files, so renv/ and .rvenvlib/ stay out.
  files <- system2(
    "git",
    c("ls-files", "-co", "--exclude-standard", "--", shQuote("*.R")),
    stdout = TRUE
  )
}
files <- unique(sub("^\\./", "", files))
files <- files[grepl("\\.R$", files) & file.exists(files)]

lints <- lapply(files, lintr::lint)
n <- vapply(lints, length, integer(1))

for (i in which(n > 0)) {
  print(lints[[i]])
}

if (sum(n) > 0) {
  cat(sprintf(
    "\nLint FAILED: %d lint(s) in %d file(s) (see .lintr).\n",
    sum(n),
    sum(n > 0)
  ))
  quit(status = 1)
}
cat(sprintf("Lint OK: %d file(s), 0 lints.\n", length(files)))
