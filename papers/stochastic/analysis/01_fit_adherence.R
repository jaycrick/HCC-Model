# papers/stochastic/analysis/01_fit_adherence.R
#
# Fit the 3-class "due-until-attend" adherence model (R/adherence.R) to the
# committed aggregate table data/adherence_counts.csv (columns A, F, T, n; the
# BASE is the pre-COVID cohort, follow-up capped at 2020-03-01); no raw data are
# read. Steps: model comparison (K = 1..5 classes, with and without
# Att_Low = 0; AIC, AICc, BIC, ICL), boundary and parametric-bootstrap LRT for
# the zero low class, profile-likelihood CIs, nonparametric bootstrap, plug-in
# posterior-predictive goodness of fit, sensitivity to the data-preparation
# choices, and the paper-2 Parameters CSVs.
#
# Writes: data/Parameters.csv, data/Parameters_sources.csv (from paper 1's
# files with the adherence rows replaced), data-derived/adherence_bootstrap.rds
# (regenerable; Git LFS via *.rds; pass --bootstrap to force a rerun),
# output/*.csv (model comparison, sensitivity; git-ignored), the ledger
# output/manuscript_values.csv and figures/adherence_gof.png.

# ----- Initialise -------------------------------------------------------------
# setup.R reads data/Parameters.csv when sourced, and this script creates that
# file, so source it under paper 1 and then switch to this paper.
options(hcc.paper = "2024-cancers")
source(here::here("R", "setup.R"))
options(hcc.paper = "stochastic")
source(here::here("R", "adherence.R"))
source(here::here("R", "plot_theme.R"))
library(ggplot2)

M_GAP <- 6 # months not due after an attended scan (the model's m)
B_BOOT <- 1000 # bootstrap resamples
run_bootstrap <- "--bootstrap" %in% commandArgs(trailingOnly = TRUE)
SEED <- 20261124

tab <- read.csv(paper_path("data", "adherence_counts.csv"))
stopifnot(identical(names(tab), c("A", "F", "T", "n")))
N <- sum(tab$n)
# The likelihood ignores T: fits and the bootstrap use the table aggregated over
# T (same likelihood, and the same RNG stream as a table without T).
tab_af <- aggregate(n ~ A + F, data = tab, FUN = sum)
tab_af <- tab_af[order(tab_af$A, tab_af$F), ]
message("Cohort N = ", N, ", cells = ", nrow(tab))
ledger <- list()

# ----- Model comparison -------------------------------------------------------
# K = 1..5 classes, with and without the lowest class's q fixed at 0 (K >= 2).
# df = 2K - 1 (K - 1 shares + K q's), or 2K - 2 with zero_low. For EVERY
# criterion n = number of PEOPLE (sum(tab$n) = N), never person-months or cells:
# people are the independent units. BIC is the primary criterion for the number
# of classes (Keribin 2000); AICc and ICL are reported alongside; AIC alone
# overfits mixtures. Singular BIC (Drton & Plummer 2017) is theoretically
# better suited to mixtures but is not used: the learning coefficients for this
# geometric-type mixture are not readily available.
criteria_table <- function(tb) {
  out <- do.call(
    rbind,
    lapply(1:5, function(K) {
      do.call(
        rbind,
        lapply(c(FALSE, TRUE), function(z) {
          if (z && K < 2) {
            return(NULL)
          }
          adherence_criteria(fit_adherence_k(tb, K, zero_low = z), tb)
        })
      )
    })
  )
  out$dBIC <- out$BIC - min(out$BIC)
  out
}
comp <- criteria_table(tab_af)
print(comp, digits = 6, row.names = FALSE)
comp_cov <- criteria_table(read.csv(paper_path(
  "data",
  "adherence_counts_withcovid.csv"
)))
message("With-COVID table:")
print(comp_cov, digits = 6, row.names = FALSE)
# The chosen model is the 3-class model with Att_Low fixed at 0 (preferred by
# BIC, AICc and ICL; see the README). The free 3-class fit is kept for the LRT,
# the comparison and the sensitivity table.
FIX <- list(Att_Low = 0)
fit3 <- fit_adherence(tab_af)
fit3_0 <- fit_adherence(tab_af, fix = FIX)
chosen <- fit3_0

# Att_Low = 0 against Att_Low > 0 (K = 3). Asymptotic: the null is on the
# boundary, so the reference is a 50:50 mixture of chi-square(0) and
# chi-square(1) (Self & Liang 1987, JASA 82:605-610); valid here only because
# the other two classes stay identifiable under the null.
lr <- 2 * (fit3$loglik - fit3_0$loglik)
lrt_p <- 0.5 * pchisq(lr, 1, lower.tail = FALSE)
message(sprintf("LRT Att_Low = 0: LR = %.3f, boundary p = %.3g", lr, lrt_p))
message(
  "Free 3-class fit at a bound: ",
  paste(names(fit3$par)[fit3$at_bound], collapse = ", ")
)
message("Distinct optima across starts: ", length(unique(round(fit3$starts$loglik, 4))))

# Parametric bootstrap LRT: simulate people from the null fit (Att_Low = 0)
# with each person's exact follow-up T (months after the shift rule), refit
# both models and recompute the LR. B = 200.
B_LRT <- 200
T_people <- rep(tab$T, tab$n)
set.seed(SEED + 2)
lr_boot <- vapply(
  seq_len(B_LRT),
  function(b) {
    sim <- simulate_adherence_people(fit3_0$par, T_people, M_GAP)
    tb <- aggregate(n ~ A + F, data = transform(sim, n = 1), FUN = sum)
    f1 <- fit_adherence(tb, n_starts = 6)
    f0 <- fit_adherence(tb, fix = list(Att_Low = 0), n_starts = 6)
    max(0, 2 * (f1$loglik - f0$loglik))
  },
  numeric(1)
)
lrt_p_boot <- (1 + sum(lr_boot >= lr)) / (B_LRT + 1)
message(sprintf(
  "Parametric bootstrap LRT: p = %.3f (B = %d); bootstrap LR median %.3f, 95%% %.3f",
  lrt_p_boot,
  B_LRT,
  median(lr_boot),
  quantile(lr_boot, 0.95)
))

# Who drives Att_Low > 0: people with A >= 1 and a high posterior probability
# of the low class under the FREE fit (fit3); this is what justifies fixing
# Att_Low = 0 in the chosen model.
post <- adherence_posterior(fit3$par, tab)
drv <- tab$A >= 1 & post[, "Low"] > 0.5
n_drv <- sum(tab$n[drv])
F_drv <- rep(tab$F[drv], tab$n[drv])
A_drv <- rep(tab$A[drv], tab$n[drv])
message(sprintf(
  "People with A >= 1 and P(low class) > 0.5: %d; their F: median %s, range %s-%s",
  n_drv,
  median(F_drv),
  if (n_drv) min(F_drv) else NA,
  if (n_drv) max(F_drv) else NA
))
message(
  "Their A (value:count): ",
  paste(names(table(A_drv)), table(A_drv), sep = ":", collapse = " ")
)
dir.create(paper_path("output"), showWarnings = FALSE, recursive = TRUE)
write.csv(
  comp,
  paper_path("output", "adherence_model_comparison.csv"),
  row.names = FALSE
)
write.csv(
  comp_cov,
  paper_path("output", "adherence_model_comparison_withcovid.csv"),
  row.names = FALSE
)

# ----- Profile-likelihood CIs (chosen model: Att_Low fixed at 0) ---------------
# Only the 4 free parameters are profiled; Att_Low is fixed, Pi_Med derived.
est <- chosen$par
PROF_NAMES <- c("Pi_Low", "Pi_High", "Att_Med", "Att_High")
prof <- lapply(setNames(PROF_NAMES, PROF_NAMES), function(p) {
  adherence_profile_ci(chosen, tab_af, p, fix = FIX)
})
profile_ci <- t(vapply(prof, function(r) c(r$lower, r$upper), numeric(2)))
colnames(profile_ci) <- c("LOWER", "UPPER")
print(cbind(POINT = est[PROF_NAMES], profile_ci), digits = 4)

# ----- Nonparametric bootstrap --------------------------------------------------
boot_path <- paper_path("data-derived", "adherence_bootstrap.rds")
boot <- if (!run_bootstrap && file.exists(boot_path)) {
  tryCatch(readRDS(boot_path), error = function(e) NULL)
} else {
  NULL
}
if (
  is.null(boot) ||
    boot$B != B_BOOT ||
    boot$N != N ||
    !identical(unlist(boot$fix), unlist(FIX))
) {
  set.seed(SEED)
  counts <- draw_bootstrap_counts(tab_af, B_BOOT)
  t0 <- Sys.time()
  draws <- adherence_bootstrap(tab_af, B_BOOT, list(n_starts = 10, fix = FIX), counts)
  message("Bootstrap: ", B_BOOT, " resamples in ", format(round(Sys.time() - t0, 1)))
  boot <- list(draws = draws, B = B_BOOT, N = N, seed = SEED, n_starts = 10, fix = FIX)
  dir.create(dirname(boot_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(boot, boot_path)
}
bd <- boot$draws[, ADH_PAR_NAMES]
# Pi_Med is derived, so it has no profile CI here: its CSV limits are the
# bootstrap percentile limits.
boot_ci <- t(apply(bd, 2, quantile, probs = c(0.025, 0.975)))
colnames(boot_ci) <- c("LOWER", "UPPER")
print(cbind(POINT = est, boot_ci), digits = 4)
bd_free <- bd[, setdiff(ADH_PAR_NAMES, "Att_Low")]
message(
  "Bootstrap draws with a parameter at a bound (<1e-6 or >1-1e-6), ",
  "excluding the fixed Att_Low: ",
  sum(apply(bd_free, 1, function(r) any(r <= 1e-6 | r >= 1 - 1e-6)))
)

# ----- Goodness of fit: plug-in posterior-predictive check ----------------------
# Each person's exact follow-up T (months after the shift rule) comes from the
# table, so simulated people have the observed exposure. The check conditions on
# each person's T, not on a fully independent test.
# Each person is simulated once per replicate from the fitted model (plug-in
# MLE, no parameter uncertainty); the replicated cohort has the observed size.
n_rep <- 500
people <- tab[rep(seq_len(nrow(tab)), tab$n), c("A", "F", "T")]
T_hat <- people$T
rate_cut <- function(A, T_m) {
  r <- pmin(A * M_GAP / T_m, 1) # attendance-rate proxy (scans per due-scan slot)
  cut(r, c(-Inf, 0.1, 0.8, Inf), right = FALSE, labels = c("<0.1", "0.1-0.8", ">=0.8"))
}
a_cut <- function(A) {
  cut(A, c(-Inf, 0, 2, 5, 8, Inf), labels = c("0", "1-2", "3-5", "6-8", ">8"))
}
set.seed(SEED + 1)
sims <- lapply(seq_len(n_rep), function(i) {
  s <- simulate_adherence_people(est, T_hat, M_GAP)
  list(a = table(a_cut(s$A)), rate = table(rate_cut(s$A, T_hat)))
})
sim_mat <- function(key) {
  t(vapply(sims, function(s) as.numeric(s[[key]]), numeric(length(sims[[1]][[key]]))))
}
obs <- list(
  a = as.numeric(table(a_cut(people$A))),
  rate = as.numeric(table(rate_cut(people$A, T_hat)))
)
gof_stat <- function(o, E) sum((o - E)^2 / pmax(E, 0.5))
gof <- lapply(c(a = "a", rate = "rate"), function(key) {
  S <- sim_mat(key)
  E <- colMeans(S)
  stat_obs <- gof_stat(obs[[key]], E)
  stat_rep <- apply(S, 1, gof_stat, E = E)
  list(
    E = E,
    lo = apply(S, 2, quantile, 0.025),
    hi = apply(S, 2, quantile, 0.975),
    p = mean(stat_rep >= stat_obs),
    stat = stat_obs
  )
})
labs <- list(a = levels(a_cut(0)), rate = levels(rate_cut(0, 1)))
gof_df <- do.call(
  rbind,
  lapply(names(gof), function(key) {
    data.frame(
      panel = c(a = "Scans attended (A)", rate = "Attendance-rate proxy A*m/T")[[key]],
      bin = factor(labs[[key]], levels = labs[[key]]),
      observed = obs[[key]],
      expected = gof[[key]]$E,
      lo = gof[[key]]$lo,
      hi = gof[[key]]$hi
    )
  })
)
print(gof_df, digits = 3, row.names = FALSE)
message(sprintf("GOF p (A bins) = %.3f; p (rate bins) = %.3f", gof$a$p, gof$rate$p))
write.csv(gof_df, paper_path("output", "adherence_gof.csv"), row.names = FALSE)
p_gof <- ggplot(gof_df, aes(bin)) +
  geom_col(aes(y = observed), fill = "grey80", colour = "black") +
  geom_pointrange(aes(y = expected, ymin = lo, ymax = hi)) +
  facet_wrap(~panel, scales = "free_x") +
  labs(
    x = NULL,
    y = "People",
    title = "Observed (bars) vs fitted-model simulation (mean, 95% range)"
  ) +
  theme_hcc()
save_figure(p_gof, "adherence_gof", 8, 4, "png", dpi = 300)

# ----- Sensitivity to data preparation ------------------------------------------
variants <- c("base", "withcovid", "diagscan", "noshift", "floor")
sens <- do.call(
  rbind,
  lapply(variants, function(v) {
    f <- if (v == "base") {
      "adherence_counts.csv"
    } else {
      paste0("adherence_counts_", v, ".csv")
    }
    tv <- read.csv(paper_path("data", f))
    fv <- fit_adherence(tv, fix = FIX)
    data.frame(variant = v, N = sum(tv$n), t(fv$par), loglik = fv$loglik)
  })
)
# The free 3-class model on the base table, as a sensitivity row.
sens <- rbind(
  sens,
  data.frame(variant = "free Att_Low", N = N, t(fit3$par), loglik = fit3$loglik)
)
print(sens, digits = 4, row.names = FALSE)
write.csv(sens, paper_path("output", "adherence_sensitivity.csv"), row.names = FALSE)

# ----- Paper-2 Parameters CSVs --------------------------------------------------
# Paper 1's files are copied with the adherence rows replaced; values come from
# the fit above (bootstrap percentile 95% CI), never typed by hand. Profile CIs
# are still computed and recorded in the ledger.
fmt <- function(x) sprintf("%.12g", x)
p1 <- read.csv(
  here::here("papers", "2024-cancers", "data", "Parameters.csv"),
  colClasses = "character",
  fileEncoding = "UTF-8-BOM"
)
ci_lo <- boot_ci[, "LOWER"]
ci_hi <- boot_ci[, "UPPER"]
new_rows <- data.frame(
  PARAMETER = ADH_PAR_NAMES,
  POINT = fmt(est[ADH_PAR_NAMES]),
  LOWER = fmt(ci_lo[ADH_PAR_NAMES]),
  UPPER = fmt(ci_hi[ADH_PAR_NAMES]),
  stringsAsFactors = FALSE
)
# Paper 2 drops paper 1's p_NA, p_CA, p_Att so the two adherence models cannot
# be mixed; the new rows take the place of the old block.
i_old <- match(c("p_NA", "p_CA", "p_Att"), p1$PARAMETER)
stopifnot(!anyNA(i_old), identical(i_old, min(i_old):max(i_old)))
p2 <- rbind(p1[seq_len(min(i_old) - 1), ], new_rows, p1[-seq_len(max(i_old)), ])
csv_field <- function(x) {
  if (grepl("[\",\n]", x)) paste0("\"", gsub("\"", "\"\"", x), "\"") else x
}
write_csv_bom <- function(df, path, bom) {
  rows <- apply(df, 1, function(r) {
    paste(vapply(r, csv_field, character(1)), collapse = ",")
  })
  txt <- paste(c(paste(names(df), collapse = ","), rows), collapse = "\n")
  cat(if (bom) "\ufeff", txt, "\n", file = path, sep = "")
}
write_csv_bom(p2, paper_path("data", "Parameters.csv"), bom = TRUE)

src1 <- read.csv(
  here::here("papers", "2024-cancers", "data", "Parameters_sources.csv"),
  colClasses = "character"
)
src_text <- paste0(
  "fitted (MLE, 3-class due-until-attend model with Att_Low fixed at 0; EM; ",
  "nonparametric bootstrap percentile 95% CI, B = 1000) ",
  "to unpublished tertiary-hospital attendance data ",
  "(aggregate table data/adherence_counts.csv)"
)
ref <- c(
  Pi_Low = paste(
    "Latent-class share (low attenders); replaces paper 1 p_NA (an attendance-rate",
    "cut). LOWER/UPPER: bootstrap percentile"
  ),
  Pi_Med = paste(
    "Latent-class share (medium attenders) = 1 - Pi_Low - Pi_High; replaces the",
    "derived paper 1 partially-adherent share. LOWER/UPPER: bootstrap percentile"
  ),
  Pi_High = paste(
    "Latent-class share (high attenders); replaces paper 1 p_CA (an attendance-rate",
    "cut). LOWER/UPPER: bootstrap percentile"
  ),
  Att_Low = paste(
    "fixed at 0 (never attends at the surveillance hospital): preferred by BIC,",
    "AICc and ICL; boundary LRT p = 0.14, parametric bootstrap LRT p = 0.095;",
    "see papers/stochastic/README.md"
  ),
  Att_Med = paste(
    "MONTHLY probability of attending while due, medium class; replaces paper 1",
    "p_Att (per scheduled scan)"
  ),
  Att_High = paste(
    "MONTHLY probability of attending while due, high class; new in paper 2"
  )
)
new_src <- data.frame(
  PARAMETER = names(ref),
  SOURCE = src_text,
  TABLE_REF = unname(ref),
  stringsAsFactors = FALSE
)
new_src <- new_src[match(ADH_PAR_NAMES, new_src$PARAMETER), ]
j_old <- match(c("p_NA", "p_CA", "p_Att"), src1$PARAMETER)
stopifnot(!anyNA(j_old), identical(j_old, min(j_old):max(j_old)))
src2 <- rbind(src1[seq_len(min(j_old) - 1), ], new_src, src1[-seq_len(max(j_old)), ])
write_csv_bom(src2, paper_path("data", "Parameters_sources.csv"), bom = FALSE)
chk <- unlist(load_params(paper_path("data", "Parameters.csv"))[ADH_PAR_NAMES])
stopifnot(all(chk >= 0), abs(sum(chk[ADH_PI_NAMES]) - 1) < 1e-8)

# ----- Ledger ---------------------------------------------------------------------
vals <- c(
  adh_N = N,
  adh_lrt_stat = lr,
  adh_lrt_p = lrt_p,
  adh_boot_B = boot$B,
  adh_gof_p_A = gof$a$p,
  adh_gof_p_rate = gof$rate$p
)
vals["adh_lrt_p_boot"] <- lrt_p_boot
vals["adh_lrt_boot_B"] <- B_LRT
vals["adh_diag_n_A1_lowclass"] <- n_drv
for (j in seq_len(nrow(comp))) {
  key <- paste0("K", comp$K[j], if (comp$zero_low[j]) "z" else "")
  for (cn in c("loglik", "df", "AIC", "AICc", "BIC", "ICL")) {
    vals[paste0("adh_", cn, "_", key)] <- comp[j, cn]
    vals[paste0("adh_", cn, "_", key, "_withcovid")] <- comp_cov[j, cn]
  }
}
for (p in ADH_PAR_NAMES) {
  vals[paste0("adh_", p)] <- est[[p]]
  if (p %in% PROF_NAMES) {
    vals[paste0("adh_", p, "_prof_lo")] <- profile_ci[p, "LOWER"]
    vals[paste0("adh_", p, "_prof_hi")] <- profile_ci[p, "UPPER"]
  }
  vals[paste0("adh_", p, "_boot_lo")] <- boot_ci[p, "LOWER"]
  vals[paste0("adh_", p, "_boot_hi")] <- boot_ci[p, "UPPER"]
}
record_values(vals)
message("Adherence fit done.")
