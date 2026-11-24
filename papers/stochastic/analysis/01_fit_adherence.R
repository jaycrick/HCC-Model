# papers/stochastic/analysis/01_fit_adherence.R
#
# Fit the 3-class "due-until-attend" adherence model (R/adherence.R) to the
# committed aggregate table data/adherence_counts.csv (columns A, F, n); no raw
# data are read. Steps: model comparison (1, 2, 3 classes; Att_Low = 0),
# profile-likelihood CIs, nonparametric bootstrap, plug-in posterior-predictive
# goodness of fit, sensitivity to the data-preparation choices, and the
# paper-2 Parameters CSVs.
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
stopifnot(identical(names(tab), c("A", "F", "n")))
N <- sum(tab$n)
message("Cohort N = ", N, ", cells = ", nrow(tab))
ledger <- list()

# ----- Model comparison -------------------------------------------------------
# k = number of identifiable parameters. Fixed values of unused classes (their
# q, when the class has probability 0) are placeholders and not counted.
fit3 <- fit_adherence(tab)
fit3_0 <- fit_adherence(tab, fix = list(Att_Low = 0))
fit2 <- fit_adherence(tab, method = "optim", fix = list(Pi_High = 0, Att_High = 0.5))
fit1 <- fit_adherence(
  tab,
  method = "optim",
  fix = list(Pi_Low = 0, Pi_High = 0, Att_Low = 0.5, Att_High = 0.5)
)
comp <- data.frame(
  model = c("1 class", "2 classes", "3 classes, Att_Low = 0", "3 classes"),
  loglik = c(fit1$loglik, fit2$loglik, fit3_0$loglik, fit3$loglik),
  k = c(1, 3, 4, 5)
)
comp$AIC <- -2 * comp$loglik + 2 * comp$k
comp$BIC <- -2 * comp$loglik + comp$k * log(N)
print(comp, digits = 6, row.names = FALSE)
stopifnot(which.min(comp$BIC) == 4, which.min(comp$AIC) == 4)

# LRT for Att_Low = 0: the null is on the boundary of the parameter space, so
# the reference distribution is a 50:50 mixture of chi-square(0) and
# chi-square(1) (Self & Liang 1987, JASA 82:605-610).
lr <- 2 * (fit3$loglik - fit3_0$loglik)
lrt_p <- 0.5 * pchisq(lr, 1, lower.tail = FALSE)
message(sprintf("LRT Att_Low = 0: LR = %.3f, boundary p = %.3g", lr, lrt_p))
message(
  "Free 3-class fit at a bound: ",
  paste(names(fit3$par)[fit3$at_bound], collapse = ", ")
)
message("Distinct optima across starts: ", length(unique(round(fit3$starts$loglik, 4))))
dir.create(paper_path("output"), showWarnings = FALSE, recursive = TRUE)
write.csv(
  comp,
  paper_path("output", "adherence_model_comparison.csv"),
  row.names = FALSE
)

# ----- Profile-likelihood CIs (chosen model: free 3-class) ---------------------
est <- fit3$par
prof <- lapply(setNames(ADH_FREE_NAMES, ADH_FREE_NAMES), function(p) {
  adherence_profile_ci(fit3, tab, p)
})
profile_ci <- t(vapply(prof, function(r) c(r$lower, r$upper), numeric(2)))
colnames(profile_ci) <- c("LOWER", "UPPER")
print(cbind(POINT = est[ADH_FREE_NAMES], profile_ci), digits = 4)

# ----- Nonparametric bootstrap --------------------------------------------------
boot_path <- paper_path("data-derived", "adherence_bootstrap.rds")
boot <- if (!run_bootstrap && file.exists(boot_path)) {
  tryCatch(readRDS(boot_path), error = function(e) NULL)
} else {
  NULL
}
if (is.null(boot) || boot$B != B_BOOT || boot$N != N) {
  set.seed(SEED)
  counts <- draw_bootstrap_counts(tab, B_BOOT)
  t0 <- Sys.time()
  draws <- adherence_bootstrap(tab, B_BOOT, list(n_starts = 10), counts)
  message("Bootstrap: ", B_BOOT, " resamples in ", format(round(Sys.time() - t0, 1)))
  boot <- list(draws = draws, B = B_BOOT, N = N, seed = SEED, n_starts = 10)
  dir.create(dirname(boot_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(boot, boot_path)
}
bd <- boot$draws[, ADH_PAR_NAMES]
# Pi_Med is derived, so it has no profile CI here: its CSV limits are the
# bootstrap percentile limits.
boot_ci <- t(apply(bd, 2, quantile, probs = c(0.025, 0.975)))
colnames(boot_ci) <- c("LOWER", "UPPER")
print(cbind(POINT = est, boot_ci), digits = 4)
message(
  "Bootstrap draws with a parameter at a bound (<1e-6 or >1-1e-6): ",
  sum(apply(bd, 1, function(r) any(r <= 1e-6 | r >= 1 - 1e-6)))
)

# ----- Goodness of fit: plug-in posterior-predictive check ----------------------
# The table gives A and F but not the follow-up length. Approximation: exposure
# T_i = F_i + M_GAP * A_i (due months plus the not-due months after each
# scan). This is exact unless follow-up ended inside the not-due period after
# the last scan, in which case T_i is too large by up to M_GAP - 1 months
# (so simulated people have slightly too much exposure: biased towards more
# scans). It also conditions on a function of the data, so the check is of the
# fit given each person's approximate exposure, not a fully independent test.
# Each person is simulated once per replicate from the fitted model (plug-in
# MLE, no parameter uncertainty); the replicated cohort has the observed size.
n_rep <- 500
people <- tab[rep(seq_len(nrow(tab)), tab$n), c("A", "F")]
T_hat <- people$F + M_GAP * people$A
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
variants <- c("base", "diagscan", "noshift", "floor")
sens <- do.call(
  rbind,
  lapply(variants, function(v) {
    f <- if (v == "base") {
      "adherence_counts.csv"
    } else {
      paste0("adherence_counts_", v, ".csv")
    }
    tv <- read.csv(paper_path("data", f))
    fv <- fit_adherence(tv)
    data.frame(variant = v, N = sum(tv$n), t(fv$par), loglik = fv$loglik)
  })
)
print(sens, digits = 4, row.names = FALSE)
write.csv(sens, paper_path("output", "adherence_sensitivity.csv"), row.names = FALSE)

# ----- Paper-2 Parameters CSVs --------------------------------------------------
# Paper 1's files are copied with the adherence rows replaced; values come from
# the fit above (profile-likelihood CI), never typed by hand.
fmt <- function(x) sprintf("%.12g", x)
p1 <- read.csv(
  here::here("papers", "2024-cancers", "data", "Parameters.csv"),
  colClasses = "character",
  fileEncoding = "UTF-8-BOM"
)
ci_lo <- c(profile_ci[, "LOWER"], Pi_Med = boot_ci["Pi_Med", "LOWER"])
ci_hi <- c(profile_ci[, "UPPER"], Pi_Med = boot_ci["Pi_Med", "UPPER"])
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
  "fitted (MLE, 3-class due-until-attend model; EM, profile-likelihood 95% CI) ",
  "to unpublished tertiary-hospital attendance data ",
  "(aggregate table data/adherence_counts.csv)"
)
ref <- c(
  Pi_Low = paste(
    "Latent-class share (low attenders); replaces paper 1 p_NA (an attendance-rate",
    "cut). LOWER/UPPER: profile likelihood"
  ),
  Pi_Med = paste(
    "Latent-class share (medium attenders) = 1 - Pi_Low - Pi_High; replaces the",
    "derived paper 1 partially-adherent share. LOWER/UPPER: bootstrap percentile",
    "(derived, no profile)"
  ),
  Pi_High = paste(
    "Latent-class share (high attenders); replaces paper 1 p_CA (an attendance-rate",
    "cut). LOWER/UPPER: profile likelihood"
  ),
  Att_Low = paste(
    "MONTHLY probability of attending while due, low class; new in paper 2"
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
for (j in seq_len(nrow(comp))) {
  key <- c("1c", "2c", "3c_pna0", "3c")[j]
  vals[paste0("adh_", c("loglik", "AIC", "BIC"), "_", key)] <- unlist(comp[
    j,
    c("loglik", "AIC", "BIC")
  ])
}
for (p in ADH_PAR_NAMES) {
  vals[paste0("adh_", p)] <- est[[p]]
  if (p %in% ADH_FREE_NAMES) {
    vals[paste0("adh_", p, "_prof_lo")] <- profile_ci[p, "LOWER"]
    vals[paste0("adh_", p, "_prof_hi")] <- profile_ci[p, "UPPER"]
  }
  vals[paste0("adh_", p, "_boot_lo")] <- boot_ci[p, "LOWER"]
  vals[paste0("adh_", p, "_boot_hi")] <- boot_ci[p, "UPPER"]
}
record_values(vals)
message("Adherence fit done.")
