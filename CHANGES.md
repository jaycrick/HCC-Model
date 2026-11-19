# Changes that alter model output

## G1 - untreated progression double conversion (finding C1)

`fromf4()` (now `progress_untreated()`) converted `p_0..p_D` from annual to
monthly, but those parameters are already monthly (Cumming et al. 2024,
Cancers 16:2745, Table 1; `data/Parameters.csv`; the conversion was removed
from the screening engine in 2023 but not here). The untreated-progression
figures (`risk_no_int`, `years_no_int`) were therefore about 12x too slow.
Published scenario results (`life_add()` and friends) are unaffected.

`progress_untreated(n = 100, t = 120)`, before and after:

| State | Before | After |
|---|---|---|
| Cirrhosis | 80.05501585576 | 80.0550158558 |
| BCLC 0 | 13.38155807791 | 1.9015484395 |
| BCLC A | 5.40591015046 | 2.8226122219 |
| BCLC B | 0.50530853727 | 0.4573605702 |
| BCLC C | 0.30876059923 | 0.4483577380 |
| BCLC D | 0.25227005443 | 0.9921202126 |
| Mortality | 0.09117672493 | 13.3229849620 |

Test: `progress_untreated(n, t)` must equal the undiagnosed states plus deaths
of `screen_cohort(Pr_Att = 0, Sym_* = 0)` (tests/testthat/test-core.R). Only the
`progress_untreated_100_120` entry of `tests/testthat/fixtures/golden.rds` was
regenerated (key renamed from `fromf4_100_120` earlier in the same step with its
value untouched).
