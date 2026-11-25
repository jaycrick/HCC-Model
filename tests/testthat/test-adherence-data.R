# Unit tests for R/adherence_data.R on synthetic sequences (no raw data).

reg <- regularise_us_months

test_that("well-spaced scans are unchanged", {
  r <- reg(c(1, 7, 14), 40)
  expect_equal(r$us_months, c(1, 7, 14))
  expect_equal(r$end_month, 40)
})

test_that("a gap < 3 months merges and shifts later scans and the end earlier", {
  r <- reg(c(4, 6, 12), 40)
  expect_equal(r$us_months, c(4, 10))
  expect_equal(r$end_month, 38)
})

test_that("a gap in [3, 6) months pushes to +6 and shifts later dates later", {
  r <- reg(c(4, 8, 20), 40)
  expect_equal(r$us_months, c(4, 10, 22))
  expect_equal(r$end_month, 42)
})

test_that("gaps of exactly 3 push (not merge); 2 merges; 6 is unchanged", {
  expect_equal(reg(c(1, 4), 30)$us_months, c(1, 7))
  expect_equal(reg(c(1, 4), 30)$end_month, 33)
  expect_equal(reg(c(1, 3), 30)$us_months, 1)
  expect_equal(reg(c(1, 3), 30)$end_month, 28)
  expect_equal(reg(c(1, 7), 30)$us_months, c(1, 7))
  expect_equal(reg(c(1, 7), 30)$end_month, 30)
})

test_that("Jan, Feb, Aug -> Jan, Jul with the end one month earlier", {
  r <- reg(c(1, 2, 8), 24)
  expect_equal(r$us_months, c(1, 7))
  expect_equal(r$end_month, 23)
  # by hand: F = 1 (month 0 to 1), d = 7; scan at 7 adds 0, d = 13;
  # end 23 adds 10 -> F = 11; two scans attended
  expect_equal(due_exposure(r$us_months, r$end_month), c(A = 2, F = 11))
})

test_that("Jan, Jun, Dec -> Jan, Jul, Jan (+1 year) with the end one month later", {
  r <- reg(c(1, 6, 12), 24)
  expect_equal(r$us_months, c(1, 7, 13))
  expect_equal(r$end_month, 25)
  expect_equal(r$us_months, reg(c(1, 7, 13), 24)$us_months)
  # no drift: 5-month gaps are pushed by exactly one month each time
  expect_equal(reg(c(0, 5, 10, 15), 30)$us_months, c(0, 6, 12, 18))
})

test_that("shifts cascade: a push can turn a later gap into a merge or a push", {
  # scan 2 pushed 4 -> 6 (+2), so scan 3 moves 7 -> 9; gap to scan 2 is 3 -> push
  r <- reg(c(0, 4, 7), 30)
  expect_equal(r$us_months, c(0, 6, 12))
  expect_equal(r$end_month, 30 + 2 + 3)
  # push then merge: scan 3 at 6 becomes 8, gap 2 -> merged, end shifts +2 -2
  r2 <- reg(c(0, 4, 6), 30)
  expect_equal(r2$us_months, c(0, 6))
  expect_equal(r2$end_month, 30)
  # merge then push: scan 2 merged (1), scan 3 at 4 becomes 3 -> pushed
  r3 <- reg(c(0, 1, 4), 30)
  expect_equal(r3$us_months, c(0, 6))
  expect_equal(r3$end_month, 30 - 1 + 3)
})

test_that("two scans in the same month merge", {
  r <- reg(c(3, 3, 10), 20)
  expect_equal(r$us_months, c(3, 10))
  expect_equal(r$end_month, 20)
})

test_that("the first scan is never shifted", {
  expect_equal(reg(c(1, 2), 50)$us_months, 1)
  expect_equal(reg(2, 50)$us_months, 2)
})

test_that("regularise_us_months handles no scans", {
  r <- reg(numeric(0), 50)
  expect_length(r$us_months, 0)
  expect_equal(r$end_month, 50)
})

test_that("due_exposure matches hand calculations", {
  expect_equal(due_exposure(numeric(0), 40), c(A = 0, F = 40))
  expect_equal(due_exposure(0, 20), c(A = 1, F = 14))
  expect_equal(due_exposure(0, 6), c(A = 1, F = 0))
  # scans at 3 and 12: F = 3 + (12 - 9) = 6, then d = 18, end 30 -> +12
  expect_equal(due_exposure(c(3, 12), 30), c(A = 2, F = 18))
  # censored mid-refractory period: F never negative
  expect_equal(due_exposure(c(3, 12), 15), c(A = 2, F = 6))
})

test_that("due_exposure rejects scans inside the refractory period", {
  expect_error(due_exposure(c(3, 8), 30))
})

test_that("adherence_counts aggregates and sorts", {
  tab <- adherence_counts(c(1, 0, 1, 1, 0), c(5, 9, 5, 2, 9), c(12, 9, 12, 8, 9))
  expect_equal(names(tab), c("A", "F", "T", "n"))
  expect_equal(tab$A, c(0, 1, 1))
  expect_equal(tab$F, c(9, 2, 5))
  expect_equal(tab$T, c(9, 8, 12))
  expect_equal(tab$n, c(2, 1, 2))
  expect_equal(sum(tab$n), 5)
  # same (A, F) with different T stays in separate cells
  t2 <- adherence_counts(c(1, 1), c(5, 5), c(12, 14))
  expect_equal(nrow(t2), 2)
})

test_that("check_disclosure accepts a clean table and rejects extras", {
  tab <- adherence_counts(c(1, 0), c(5, 9), c(12, 9))
  expect_message(check_disclosure(tab), "min n = 1")
  expect_error(suppressMessages(check_disclosure(cbind(tab, id = 1:2))))
  expect_error(suppressMessages(check_disclosure(transform(tab, n = n + 0.5))))
  expect_error(suppressMessages(check_disclosure(transform(tab, n = -n))))
})

test_that("drop_in_gap ignores scans inside the not-due period", {
  expect_equal(drop_in_gap(c(0, 3, 6, 11, 12, 18), 6), c(0, 6, 12, 18))
  expect_equal(drop_in_gap(numeric(0), 6), numeric(0))
  expect_equal(drop_in_gap(c(2, 8), 6), c(2, 8))
  # exact spacing is kept; the dropped scan does not reset the clock
  expect_equal(drop_in_gap(c(1, 5, 7), 6), c(1, 7))
})
