# Manuscript-values ledger: record_values() mechanics, and the values that can
# be recomputed cheaply (02 and 04) against the ledger baseline fixture
# (fixtures/ledger_baseline.rds, regenerated under G9). The published values
# stay frozen in data-derived/param_names.rds, which is never rewritten.

test_that("record_values round-trips exactly", {
  path <- file.path(withr::local_tempdir(), "ledger.csv")
  v <- c(a = 1 / 3, b = pi * 1e-9, c = 123456789.123456789, d = 0)
  record_values(v, path = path)
  txt <- read.csv(path, colClasses = "character")
  expect_identical(
    txt$name,
    sort(names(v), method = "radix"),
    label = "ledger names"
  )
  # The text is exact: it is a correctly rounded decimal of each double that
  # (at 17 significant digits) identifies it uniquely. R's own parser is not
  # always correctly rounded (it can be 1 ulp off), so reading back is checked
  # to 1e-15 relative, and the text is checked against the 17-digit form.
  sorted <- v[txt$name]
  expect_equal(
    as.numeric(txt$value),
    unname(sorted),
    tolerance = 1e-15,
    label = "values read back from the ledger CSV"
  )
  expect_identical(
    txt$value[txt$name == "c"],
    "123456789.12345679",
    label = "ledger text for c is the 17-digit decimal"
  )
})

test_that("record_values sorts by name and is deterministic", {
  path <- file.path(withr::local_tempdir(), "ledger.csv")
  record_values(c(zeta = 1, alpha = 2), path = path)
  record_values(c(mid = 3), path = path)
  expect_identical(
    read.csv(path)$name,
    c("alpha", "mid", "zeta"),
    label = "ledger row order"
  )
  first <- readLines(path)
  record_values(c(mid = 3), path = path)
  expect_identical(
    readLines(path),
    first,
    label = "ledger text after re-recording identical values"
  )
})

test_that("record_values overwrites a repeated name and keeps the rest", {
  path <- file.path(withr::local_tempdir(), "ledger.csv")
  record_values(c(a = 1, b = 2), path = path)
  record_values(c(a = 10), path = path)
  back <- read.csv(path)
  expect_identical(
    setNames(as.numeric(back$value), back$name),
    c(a = 10, b = 2),
    label = "ledger after overwriting a"
  )
})

test_that("record_values rejects unnamed or duplicated values", {
  path <- file.path(withr::local_tempdir(), "ledger.csv")
  expect_error(record_values(c(1, 2), path = path), label = "unnamed values")
  expect_error(
    record_values(c(a = 1, a = 2), path = path),
    label = "duplicated names"
  )
})

# ----- Baseline comparison ----------------------------------------------------
baseline <- readRDS(
  here::here("tests", "testthat", "fixtures", "ledger_baseline.rds")
)

test_that("the frozen published reference param_names.rds is unchanged", {
  frozen <- readRDS(paper_path("data-derived", "param_names.rds"))
  expect_length(frozen, 91)
  # published value as stored (commit 6902e05); G9 changed the ledger, not this
  expect_equal(
    frozen[["diag_scen_lower"]],
    3.57197106509,
    tolerance = 1e-9,
    label = "frozen param_names.rds diag_scen_lower (published, pre-G9)"
  )
  expect_true(
    all(names(baseline) %in% names(frozen)),
    label = "baseline names are a subset of the frozen names"
  )
})

# Names that may legitimately differ from the baseline (none observed so far).
# G1 (untreated progression) and G4 (Sym_B in sweeps) would be listed here
# with a comment citing the gate; nothing is loosened silently.
gate_changed <- character(0)

expect_matches_baseline <- function(got, what) {
  expect_true(
    all(names(got) %in% names(baseline)),
    label = paste(what, "names all present in baseline")
  )
  for (nm in setdiff(names(got), gate_changed)) {
    expect_equal(
      unname(got[[nm]]),
      unname(baseline[[nm]]),
      tolerance = 1e-12,
      label = paste0(what, " ledger value ", nm)
    )
  }
}

test_that("02 sensitivity ledger values equal the baseline", {
  path <- file.path(withr::local_tempdir(), "ledger.csv")
  tbl <- modifiable_sensitivity(life_tibble(csv_25))
  record_values(
    setNames(tbl$years_mid_comp, paste(tbl$Parameter, tbl$Bound, sep = "_")),
    path = path
  )
  back <- read.csv(path)
  expect_matches_baseline(setNames(back$value, back$name), "02")
})

test_that("04 scenario ledger values equal the baseline", {
  path <- file.path(withr::local_tempdir(), "ledger.csv")
  point <- run_scenarios(tibble::tibble(.rows = 1), derive_sym_b = FALSE)
  for (risk in names(risk_multipliers)) {
    sfx <- ifelse(risk == "base", "", paste0("_", risk))
    d <- scenario_diff_data(point, resampled_pars_fin, risk)
    key <- d$points$ledger_key
    record_values(
      setNames(
        c(d$points[[2]], d$points$lowCI, d$points$highCI),
        c(
          paste0(key, sfx),
          paste0(key, "_lower", sfx),
          paste0(key, "_upper", sfx)
        )
      ),
      path = path
    )
  }
  back <- read.csv(path)
  expect_gt(nrow(back), 0, label = "04 ledger rows")
  expect_matches_baseline(setNames(back$value, back$name), "04")
})
