# Cleaned adherence data (aggregate tables)

The adherence model (`R/adherence.R`) is fitted to **aggregate count tables only**.
The person-level source data are unpublished tertiary-hospital attendance data.
They are held by the authors and are not public.
The scripts that build the tables from them are kept outside this repository.
The cleaning rules are implemented in `R/adherence_data.R` and unit-tested in
`tests/testthat/test-adherence-data.R`.

## Tables

| File | Contents |
|------|----------|
| `adherence_counts.csv` | Base analysis |
| `adherence_counts_diagscan.csv` | Sensitivity: count the diagnostic scan |
| `adherence_counts_noshift.csv` | Sensitivity: no shift, ignore in-gap scans |
| `adherence_counts_floor.csv` | Sensitivity: `floor()` month conversion |
| `gap_counts.csv` | Gaps between consecutive scans before the shift rule (descriptive only) |

Every `adherence_counts*` table has exactly the columns `A`, `F`, `n`, one row per distinct pair:

- `A`: number of scans attended while due;
- `F`: number of months due without a scan;
- `n`: number of people with that pair.

`sum(n)` is the cohort size, N = 108.
The table is a sufficient statistic for the likelihood
(l = sum n log sum_c pi_c q_c^A (1 - q_c)^F), so it reproduces the fit exactly.
It holds no dates or identifiers.
Single-person cells exist (minimum cell count 1) and are covered by the approval to release.

`gap_counts.csv` (columns `gap_months`, `n`) counts the gaps between consecutive
in-window scans, in whole months, **before** the shift rule (344 gaps, 108 people).
The start-to-first-scan interval is not a gap.
Gaps of 0 to 12 months are always listed; longer gaps are listed only where present.
The counts after the rule follow from it:
gaps of 0 to 2 months are merged away, gaps of 3 to 5 months become 6, and gaps of
6 months or more are unchanged.

## Cohort (aggregate CONSORT counts)

The cohort follows the original attendance analysis (Cumming 2024, Supplement Sec. 3.1).
Rows (= people) after each step:

| Step | Rows | Excluded |
|------|------|----------|
| Joins of the four source files | 218 | |
| Cirrhosis recorded as "No" | 215 | 3 |
| SVR recorded as "No" | 196 | 19 |
| Lost to follow-up | 152 | 44 |
| HCC diagnosed before 2016-07-02 | 139 | 13 |
| Fewer than 3 expected scans (see below) | 118 | 21 |
| No scan history available (missing attendance rate) | 108 | 10 |

Rows with a missing value in the screening variable are kept at the cirrhosis, SVR and lost-to-follow-up steps.
Expected scans = floor((end - start) / 183 days); people with 2 or fewer are excluded.

## Dates

- Start = max(SVR date, completion of antiviral therapy + 84 days).
  If neither exists, start = max(earliest scan - 183 days, 2016-07-01).
- End = earliest of death (any source) and HCC diagnosis (any source).
  If neither, end = the census date, 2022-03-18 (or the registry consensus date if later).
- Scan window: scans strictly after the start and strictly before the end.
  Scans outside it (mostly before the start) are discarded.
- Scan dates are sorted and de-duplicated within a person.

## Making A and F

1. **Month conversion first**: month = round(days / 30.4375), for every scan and for the end of follow-up.
2. **Shift rule** (`regularise_us_months`) on the integer months, applied one scan at a time in date order, each scan compared with the previous *kept* scan after earlier shifts.
   With g the gap in whole months:
   - g < 3: the scan is merged into the previous one (dropped), and every later scan and the end shift earlier by g;
   - 3 <= g < 6: the scan moves to previous + 6, and every later scan and the end shift later by 6 - g (cascading);
   - g >= 6: unchanged.
   The first scan is never shifted. Two scans in the same month have g = 0 and merge.
   Worked examples (scans in calendar months, start in month 0):
   - Jan, Feb, Aug becomes Jan, Jul, and the end of follow-up moves 1 month earlier;
   - Jan, Jun, Dec becomes Jan, Jul, Jan (next year), and the end moves 1 month later.
   Converting to months before shifting avoids drift from fractional months.
3. **Due at start** (`due_exposure`): a person is due at month 0.
   Walk the scan months s in order: F += s - d, A += 1, d = s + 6 (six months not due after an attended scan).
   After the last scan, F += max(0, end - d).

## Variants

- `diagscan`: when follow-up ends at an HCC diagnosis, a scan on the diagnosis date is also counted (window `u <= end`).
  A tie between death and diagnosis is treated as an HCC end.
- `noshift`: no shift rule. Scans are converted to integer months, then a scan less than 6 months after the previous kept scan is ignored: it neither counts nor resets the clock (`drop_in_gap`), and the end is not moved.
- `floor`: as the base, but month = floor(days / 30.4375) before the same integer rule.
