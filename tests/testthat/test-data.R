# Data tests for the two parameter CSVs.

main <- read_params(csv_main)
p25 <- read_params(csv_25)

prob_params <- c(
  "p_HCC",
  "p_0",
  "p_A",
  "p_B",
  "p_C",
  "p_D",
  "Sens_0",
  "Sens_A",
  "Sens_B",
  "Sens_C",
  "Sens_D",
  "p_NA",
  "p_CA",
  "p_Att",
  "Sym_0",
  "Sym_A",
  "Sym_B",
  "Sym_C",
  "Sym_D"
)

test_that("both CSVs list PARAMETER in identical order", {
  expect_identical(
    rownames(main),
    rownames(p25),
    label = "PARAMETER order: data/Parameters.csv vs data/Parameters_25.csv"
  )
  expect_identical(
    colnames(main),
    colnames(p25),
    label = "column names: Parameters.csv vs Parameters_25.csv"
  )
})

test_that("LOWER <= POINT <= UPPER in both CSVs", {
  for (nm in c("main", "p25")) {
    d <- get(nm)
    expect_true(
      all(d$LOWER <= d$POINT),
      label = paste0(
        "LOWER <= POINT in ",
        nm,
        " (violations: ",
        paste(rownames(d)[d$LOWER > d$POINT], collapse = ","),
        ")"
      )
    )
    expect_true(
      all(d$POINT <= d$UPPER),
      label = paste0(
        "POINT <= UPPER in ",
        nm,
        " (violations: ",
        paste(rownames(d)[d$POINT > d$UPPER], collapse = ","),
        ")"
      )
    )
  }
})

test_that("probability parameters lie in [0,1]", {
  for (nm in c("main", "p25")) {
    d <- get(nm)
    expect_true(
      all(prob_params %in% rownames(d)),
      label = paste0("all probability parameters present in ", nm)
    )
    vals <- as.matrix(d[prob_params, c("POINT", "LOWER", "UPPER")])
    expect_true(
      all(vals >= 0 & vals <= 1),
      label = paste0("probability parameters in [0,1] in ", nm)
    )
  }
})

test_that("Parameters_25 POINT == Parameters POINT", {
  expect_equal(
    p25$POINT,
    main$POINT,
    tolerance = 0,
    label = "POINT column: Parameters_25.csv == Parameters.csv"
  )
})
