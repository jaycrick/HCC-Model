# Adherence data helpers (file R/adherence_data.R)
# Pure functions that turn a person's scan dates into the two sufficient
# statistics of the "due-until-attend" adherence model: A = scans attended
# while due, F = months due without a scan. No file I/O, no randomness.
# Used by the private data-preparation script; unit-tested on
# synthetic sequences only.

# Average month length in days (365.25 / 12).
MONTH_DAYS <- 30.4375

#' Regularise scan months so consecutive scans are at least `gap` months apart.
#'
#' Works on whole months: convert scans and the end of follow-up to integer
#' months first (round(days / MONTH_DAYS), or floor), then apply this rule.
#' Scans are processed in order, one at a time, each compared with the previous
#' KEPT scan after earlier shifts. With g the gap in months, g < merge_below
#' drops the scan (merged into the previous) and shifts all later scans and the
#' end earlier by g; merge_below <= g < gap moves the scan to previous + gap and
#' shifts all later scans and the end later by gap - g; g >= gap leaves it. The
#' first scan is never shifted (people are due at the start). Two scans in the
#' same month have g = 0 and merge.
#'
#' Ref: plan 2026-11-24, Sec. 7 amendments (whole-month shift rule).
#'
#' @param us_months   Sorted integer months since start, duplicates allowed.
#' @param end_month   Integer end of follow-up, in months since start.
#' @param gap         Required gap in months (the model's m).
#' @param merge_below Gaps below this many months are merged.
#' @return list(us_months, end_month) after regularisation.
regularise_us_months <- function(us_months, end_month, gap = 6, merge_below = 3) {
  x <- us_months
  n <- length(x)
  keep <- rep(TRUE, n)
  if (n >= 2) {
    prev <- x[1]
    for (j in 2:n) {
      g <- x[j] - prev
      later <- seq_len(n) > j
      if (g < merge_below) {
        keep[j] <- FALSE
        x[later] <- x[later] - g
        end_month <- end_month - g
      } else if (g < gap) {
        shift <- gap - g
        x[j] <- prev + gap
        x[later] <- x[later] + shift
        end_month <- end_month + shift
        prev <- x[j]
      } else {
        prev <- x[j]
      }
    }
  }
  list(us_months = x[keep], end_month = end_month)
}

#' Drop scans that fall inside the not-due period after the previous kept scan.
#'
#' The "no shift" alternative to regularise_us_months(): a scan less than `gap` units
#' after the previous KEPT scan is ignored (it neither counts nor resets the
#' clock); the end of follow-up is not moved. Used on integer months (gap = 6).
#'
#' Ref: plan 2026-11-24, finding F10 (sensitivity analysis); no paper reference.
#'
#' @param us   Sorted integer scan months since start.
#' @param gap  Minimum spacing between kept scans, in months.
#' @return The kept scan times.
drop_in_gap <- function(us, gap) {
  keep <- rep(TRUE, length(us))
  prev <- -Inf
  for (j in seq_along(us)) {
    if (us[j] - prev < gap - 1e-9) {
      keep[j] <- FALSE
    } else {
      prev <- us[j]
    }
  }
  us[keep]
}

#' Attended scans (A) and due-without-scan months (F) for one person.
#'
#' People are due from month 0. Walk: d = 0; for each scan month s, F += s - d,
#' A += 1, d = s + m (m months not due after attending). After the last scan,
#' F += max(0, end_month - d).
#'
#' Ref: plan 2026-11-24, Sec. 4 S1 step 3.
#'
#' @param us_months  Integer scan months since start, sorted, spaced >= m.
#' @param end_month  End of follow-up in months since start.
#' @param m          Months not due after an attended scan.
#' @return Named numeric vector c(A, F).
due_exposure <- function(us_months, end_month, m = 6) {
  d <- 0
  n_att <- 0
  n_due <- 0
  for (s in us_months) {
    stopifnot(s >= d)
    n_due <- n_due + (s - d)
    n_att <- n_att + 1
    d <- s + m
  }
  n_due <- n_due + max(0, end_month - d)
  c(A = n_att, F = n_due)
}

#' Aggregate per-person (A, F) into a count table.
#'
#' Ref: none (aggregation helper).
#'
#' @param A       Numeric vector: scans attended while due, per person.
#' @param F_due   Numeric vector of equal length: due months without a scan.
#' @return data.frame(A, F, n) sorted by A then F; n = people with that pair.
adherence_counts <- function(A, F_due) {
  stopifnot(length(A) == length(F_due))
  tab <- as.data.frame(table(A = A, Fd = F_due), stringsAsFactors = FALSE)
  tab <- tab[tab$Freq > 0, ]
  out <- data.frame(
    A = as.integer(tab$A),
    F = as.integer(tab$Fd),
    n = as.integer(tab$Freq)
  )
  out <- out[order(out$A, out$F), ]
  rownames(out) <- NULL
  out
}

#' Check that a count table is safe to release.
#'
#' Stops unless the columns are exactly A, F, n and all values are
#' non-negative whole numbers. Reports the smallest cell and the cell count.
#'
#' Ref: none (disclosure safeguard).
#'
#' @param tab  Data frame to check.
#' @return `tab`, invisibly.
check_disclosure <- function(tab) {
  stopifnot(identical(names(tab), c("A", "F", "n")))
  ok <- vapply(
    tab,
    function(v) is.numeric(v) && all(!is.na(v) & v >= 0 & v == round(v)),
    logical(1)
  )
  stopifnot(all(ok))
  message("Disclosure check: ", nrow(tab), " cells, min n = ", min(tab$n))
  invisible(tab)
}
