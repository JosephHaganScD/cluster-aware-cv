###############################################################################
# Full Factorial Simulation Study: Cluster-Aware Cross-Validation (v5)
# =============================================================================
# CHANGES FROM v4 (Diagnostic and Prognostic Research resubmission).
#   Addresses code-review findings 1-4 and 6, plus the reviewer's additional
#   design changes (fixed predictor correlation, population-level prevalence
#   calibration, common random numbers for the paired sensitivity arm).
#
#   F1  Exactly N_SIM valid datasets per condition. The >=3 events AND
#       >=3 non-events requirement is enforced INSIDE generation by regeneration
#       (attempt-capped, hard error if exceeded). The simulated population is
#       therefore explicitly conditional on at least three events and three
#       non-events. Attempts, regenerations, invalid folds, and fitting errors
#       are logged per condition. Headline summaries are the equal-weight mean
#       of the condition means (pooled means also reported).
#   F2  Data-generating mechanism now matches the Methods exactly:
#         - ONE fixed, unit-diagonal predictor correlation matrix R, built once
#           and reused across every dataset, condition, and the test sample
#           (no per-dataset redraw). Cov_b = ICC * R, so each predictor has
#           between-subject variance exactly ICC and total marginal variance
#           exactly 1.0, homogeneous across predictors.
#         - Four active coefficients of COMMON magnitude (signs +,+,-,+); the
#           magnitude is frozen by the pilot (Section P) at the reference cell.
#   F3  No silent lambda fallback. Every fit uses lambda = LAMBDA. A failed fit
#       is recorded (fit_status, message) and the dataset flagged, never
#       silently refit at a different lambda. Run reports zero-fallback status.
#   F4  Outcome-stratified folds (subjects stratified by subject-level outcome
#       for cluster folds; observations stratified by their subject-level
#       outcome for naive folds). All N_CV_REPS replicates are required to
#       yield complete out-of-fold predictions; replicate-validity is recorded,
#       and a dataset that cannot is flagged (graceful, not a global abort).
#   F6  Signal-variance-controlled sensitivity arm (SIGNAL_MODE = "controlled"):
#       within each ICC the coefficient vector is scaled by sqrt(ICC_REF/ICC),
#       holding the latent linear-predictor variance constant across ICC. Paired
#       to the baseline arm by common random numbers (same subject effects,
#       within-subject innovations, outcome uniforms, fold assignments, and
#       test-sample draws; only the coefficient scale and intercept differ).
#   ADD Population-level prevalence calibration: the intercept solves
#       E[logistic(a0 + eta)] = prevalence under the POPULATION distribution of
#       the subject linear predictor eta ~ N(0, Var_eta), via Gauss-Hermite
#       quadrature, computed once per (ICC, prevalence, signal mode). The
#       development and independent test samples are then draws from exactly the
#       same prespecified population.
#
# WORKFLOW (run in this order):
#   1. Set RUN_MODE <- "pilot" and source the script. It freezes |beta| to hit
#      the reference-cell oracle AUROC target, smoke-tests every stage on a
#      couple of cells, and writes results/frozen_magnitude.rds. Inspect the
#      printed oracle-AUROC band and regeneration counts.
#   2. Set RUN_MODE <- "main"  -> full 162-cell baseline factorial.
#   3. Set RUN_MODE <- "sensitivity" -> paired 18-cell (optionally +6) arm.
#   The main/sensitivity runs refuse to start until the pilot has frozen |beta|.
#
# Expected runtime (sequential, checkpoint/RESUME enabled): main ~24-42 h at
#   N_TEST = 2000; pilot ~minutes; sensitivity a few hours.
#
# Author: Joseph L. Hagan, ScD, MSPH
# Date:   August 2026
###############################################################################

suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
})

# ── USER SETTINGS ────────────────────────────────────────────────────────────

# Set to "pilot", "main", or "sensitivity". Pilot must be run first.
RUN_MODE <- "sensitivity"

# Output directory. IMPORTANT: use a LOCAL, non-OneDrive path for long runs to
# avoid write-collisions during checkpointing (e.g. "C:/simulations/cawcv_v5").
# OUT_DIR <- file.path("results", "simulation_v5")
OUT_DIR <- "C:/simulations/cluster aware cv/results/simulation_v5"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

N_SIM     <- 500          # VALID datasets required per condition (exact)
N_CV_REPS <- 10           # CV replicates per k-fold strategy (all must be valid)
N_TEST    <- 2000         # new subjects for the true-performance reference
RESUME    <- TRUE         # skip completed conditions on restart
LAMBDA    <- 0.01         # fixed ridge penalty (no fallback)

FIXED_N_DAYS <- 57
FIXED_N_PRED <- 8         # 4 active (+,+,-,+), 4 inactive

# Reference cell and target for the frozen common coefficient magnitude (Pilot).
ICC_REF   <- 0.3
PREV_REF  <- 0.20
AR1_REF   <- 0.5
AUROC_TARGET <- 0.65      # population ORACLE AUROC at the reference cell
AUROC_TOL    <- 0.01      # operational tolerance (0.64-0.66); pilot warns if missed

# Sign pattern of the four active coefficients (magnitude supplied by pilot).
ACTIVE_SIGNS <- c(1, 1, -1, 1)

# Regeneration cap: hard error if a single valid dataset needs more attempts.
MAX_REGEN <- 500L

# Master seed and Gauss-Hermite nodes for population-intercept calibration.
MASTER_SEED <- 20260410L
GH_NODES    <- 64L

# Optional AR1-robustness slice for the sensitivity arm (3 ICC x 2 AR1 at one N).
SENS_INCLUDE_AR1_SLICE <- TRUE
SENS_SLICE_N           <- 50

set.seed(MASTER_SEED)

# ── 1. FIXED PREDICTOR CORRELATION MATRIX R (built once) ─────────────────────
# A single positive-definite, unit-diagonal correlation matrix, reused for every
# dataset and the test sample. Built from a mild 2-factor structure with a fixed
# seed, then converted to a correlation matrix. Deterministic given R_SEED.

R_SEED <- 20260811L
build_fixed_R <- function(n_pred = FIXED_N_PRED, n_fac = 2, seed = R_SEED,
                          load_scale = 0.30, unique_var = 0.70) {
  g_state <- .Random.seed                      # preserve caller RNG stream
  set.seed(seed)
  L  <- matrix(rnorm(n_pred * n_fac), n_pred, n_fac) * load_scale
  C  <- L %*% t(L) + diag(unique_var, n_pred)
  R  <- cov2cor(C)                             # force unit diagonal
  assign(".Random.seed", g_state, envir = .GlobalEnv)
  R
}
R_PRED     <- build_fixed_R()
CHOL_R     <- chol(R_PRED)                     # upper-triangular; b = z %*% chol(Cov_b)
MAXABS_OFF <- max(abs(R_PRED[upper.tri(R_PRED)]))

# Quadratic form for the ACTIVE pattern (per unit |beta|^2): Var_eta = ICC*|b|^2*qform
active_pattern <- c(ACTIVE_SIGNS, rep(0, FIXED_N_PRED - length(ACTIVE_SIGNS)))
QFORM_ACTIVE   <- as.numeric(t(active_pattern) %*% R_PRED %*% active_pattern)

var_eta <- function(icc, bmag) icc * (bmag^2) * QFORM_ACTIVE

# ── 2. POPULATION-LEVEL PREVALENCE CALIBRATION (Gauss-Hermite) ───────────────
# Solve E[logistic(a0 + eta)] = prev, eta ~ N(0, v), for the intercept a0.

.gh <- (function(n) {
  # Probabilists' Hermite nodes/weights (weight exp(-x^2/2)); E[g(eta)] with
  # eta ~ N(0,v) == sum(w * g(sqrt(v) * x)) / sqrt(2*pi).
  # Golub-Welsch via the symmetric tridiagonal companion matrix.
  i <- 1:(n - 1)
  J <- matrix(0, n, n); J[cbind(i, i + 1)] <- sqrt(i); J[cbind(i + 1, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  x <- e$values
  w <- (e$vectors[1, ]^2) * sqrt(2 * pi)       # so sum(w)/sqrt(2*pi) = 1
  list(x = x, w = w)
})(GH_NODES)

pop_prevalence <- function(a0, v) {
  eta <- sqrt(v) * .gh$x
  sum(.gh$w * plogis(a0 + eta)) / sqrt(2 * pi)
}
calibrate_intercept <- function(prev, v) {
  if (v <= 0) return(qlogis(prev))
  uniroot(function(a0) pop_prevalence(a0, v) - prev,
          lower = -25, upper = 25, tol = 1e-10)$root
}

# ── 3. COEFFICIENT VECTOR BY SIGNAL MODE ─────────────────────────────────────
# baseline:   |beta| = bmag_ref (frozen), NOT rescaled across ICC.
# controlled: |beta| = bmag_ref * sqrt(ICC_REF / icc)  -> Var_eta constant.

make_true_beta <- function(bmag_ref, icc, signal_mode) {
  scale <- if (signal_mode == "controlled") sqrt(ICC_REF / icc) else 1.0
  bmag  <- bmag_ref * scale
  beta  <- numeric(FIXED_N_PRED)
  beta[seq_along(ACTIVE_SIGNS)] <- ACTIVE_SIGNS * bmag
  beta
}

# ── 4. DATA-GENERATING MECHANISM ─────────────────────────────────────────────
# Draws latent quantities in a FIXED order so that, for a given seed, baseline
# and signal-controlled arms share subject effects, innovations, and outcome
# uniforms (common random numbers). Only true_beta and intercept differ.
#
# Returns exactly one dataset satisfying >=3 events and >=3 non-events, using
# regeneration that advances the seed; regen count is returned. For the paired
# sensitivity arm, generate_pair() (Section 11) keeps the two arms in lockstep.

generate_core <- function(n_subjects, n_days, n_pred, icc, ar1_rho,
                          true_beta, intercept, seed) {
  set.seed(seed)
  sigma_e     <- sqrt(1.0 - icc)
  sigma_innov <- sqrt((1.0 - icc) * (1.0 - ar1_rho^2))
  Cov_b_chol  <- CHOL_R * sqrt(icc)            # chol(ICC * R) = sqrt(ICC) * chol(R)

  # (a) subject random effects; (b) subject linear predictor; (c) outcome uniforms
  z_b   <- matrix(rnorm(n_subjects * n_pred), n_subjects, n_pred)
  b_mat <- z_b %*% Cov_b_chol
  eta   <- as.numeric(b_mat %*% true_beta)
  probs <- plogis(intercept + eta)
  u_y   <- runif(n_subjects)                   # shared across arms (CRN)
  y     <- as.integer(u_y < probs)

  # (d) within-subject AR(1) predictor trajectories
  total_rows <- n_subjects * n_days
  X_mat    <- matrix(NA_real_, total_rows, n_pred)
  case_ids <- integer(total_rows)
  dol_vec  <- integer(total_rows)
  y_obs    <- integer(total_rows)
  row_idx  <- 0L
  for (i in seq_len(n_subjects)) {
    e_prev <- rnorm(n_pred) * sigma_e
    for (t in seq_len(n_days)) {
      innov  <- rnorm(n_pred) * sigma_innov
      e_curr <- ar1_rho * e_prev + innov
      row_idx <- row_idx + 1L
      X_mat[row_idx, ] <- b_mat[i, ] + e_curr
      case_ids[row_idx] <- i
      dol_vec[row_idx]  <- t
      y_obs[row_idx]    <- y[i]
      e_prev <- e_curr
    }
  }
  X_full <- cbind(DOL = dol_vec, X_mat)
  colnames(X_full) <- c("DOL", paste0("X", seq_len(n_pred)))

  struct <- list(chol_b = Cov_b_chol, true_beta = true_beta, intercept = intercept,
                 sigma_e = sigma_e, sigma_innov = sigma_innov,
                 ar1_rho = ar1_rho, n_pred = n_pred)

  list(X = X_full, y_obs = y_obs, case_ids = case_ids, y_subject = y,
       n_subjects = n_subjects, n_events = sum(y), struct = struct,
       seed_used = seed)
}

# Single-arm generator with regeneration to satisfy the event constraint.
generate_dataset <- function(n_subjects, n_days, n_pred, icc, ar1_rho,
                             true_beta, intercept, seed) {
  regen <- 0L
  repeat {
    dat <- generate_core(n_subjects, n_days, n_pred, icc, ar1_rho,
                         true_beta, intercept, seed + regen * 100003L)
    if (dat$n_events >= 3 && dat$n_events <= n_subjects - 3) {
      dat$regen <- regen
      return(dat)
    }
    regen <- regen + 1L
    if (regen > MAX_REGEN)
      stop(sprintf("MAX_REGEN exceeded (N=%d, ICC=%.1f, prev target unmet)",
                   n_subjects, icc))
  }
}

# ── 5. TRUE-PERFORMANCE REFERENCE ────────────────────────────────────────────
# Full-data (deployed) model fit, and its true out-of-sample AUROC on a fresh
# matched sample of N_TEST new subjects from the SAME prespecified population.
# No lambda fallback: on error, return status so the caller flags the dataset.

fit_ridge_full_coef <- function(X, y_obs, lambda = LAMBDA) {
  n_pos <- sum(y_obs == 1); n_neg <- sum(y_obs == 0); n_total <- length(y_obs)
  w <- ifelse(y_obs == 1, n_total / (2 * n_pos), n_total / (2 * n_neg))
  out <- tryCatch(
    list(coef = as.numeric(coef(glmnet(X, y_obs, family = "binomial", alpha = 0,
                                       weights = w, lambda = lambda,
                                       standardize = TRUE), s = lambda)),
         status = "ok", msg = ""),
    error = function(e) list(coef = NULL, status = "error", msg = conditionMessage(e)))
  out
}

compute_true_auroc <- function(struct, coefs, n_test, n_days, seed) {
  set.seed(seed)
  n_pred <- struct$n_pred
  b_test <- matrix(rnorm(n_test * n_pred), n_test, n_pred) %*% struct$chol_b
  lin    <- as.numeric(b_test %*% struct$true_beta)
  probs  <- plogis(struct$intercept + lin)
  y_test <- as.integer(runif(n_test) < probs)
  if (length(unique(y_test)) < 2) return(NA_real_)

  a0 <- coefs[1]; beta <- coefs[-1]            # [DOL, X1..Xn_pred]
  sum_p <- numeric(n_test)
  e <- matrix(rnorm(n_test * n_pred), n_test, n_pred) * struct$sigma_e
  for (t in seq_len(n_days)) {
    innov <- matrix(rnorm(n_test * n_pred), n_test, n_pred) * struct$sigma_innov
    e     <- struct$ar1_rho * e + innov
    X_row <- cbind(t, b_test + e)
    sum_p <- sum_p + plogis(a0 + as.numeric(X_row %*% beta))
  }
  subj_pred <- sum_p / n_days
  as.numeric(pROC::auc(pROC::roc(y_test, subj_pred, quiet = TRUE)))
}

# ── 6. EMPIRICAL DATA CHARACTERISTICS (descriptive; unequal-size ICC) ────────
compute_empirical_chars <- function(dat) {
  X <- dat$X; ids <- dat$case_ids
  unique_ids <- unique(ids); n_obs <- nrow(X)
  iccs <- numeric(4)
  for (p in 1:4) {
    vals <- X[, p + 1]
    gm <- tapply(vals, ids, mean); gs <- tapply(vals, ids, length)
    k <- length(gm); grand <- mean(vals)
    ms_b <- sum(gs * (gm - grand)^2) / (k - 1)
    wr   <- vals - gm[as.character(ids)]
    ms_w <- sum(wr^2) / (n_obs - k)
    n0   <- (n_obs - sum(gs^2) / n_obs) / (k - 1)
    s2b  <- max((ms_b - ms_w) / n0, 0)
    iccs[p] <- s2b / (s2b + ms_w)
  }
  # lag-1 autocorrelation on DOL-ordered rows (explicit ordering)
  ar1s <- numeric(0)
  for (s in unique_ids[1:min(50, length(unique_ids))]) {
    idx <- which(ids == s); idx <- idx[order(X[idx, "DOL"])]
    for (p in 1:2) {
      v <- X[idx, p + 1]
      if (length(v) > 5) { ac <- cor(v[-length(v)], v[-1]); if (!is.na(ac)) ar1s <- c(ar1s, ac) }
    }
  }
  list(empirical_icc = mean(iccs),
       empirical_ar1 = if (length(ar1s)) mean(ar1s) else NA_real_)
}

# ── 7. CV STRATEGY FUNCTIONS (stratified folds; no lambda fallback) ──────────

fit_predict_ridge <- function(X_train, y_train, X_test) {
  n_pos <- sum(y_train == 1); n_neg <- sum(y_train == 0); n_total <- length(y_train)
  w <- ifelse(y_train == 1, n_total / (2 * n_pos), n_total / (2 * n_neg))
  tryCatch(
    list(pred = as.numeric(predict(
           glmnet(X_train, y_train, family = "binomial", alpha = 0, weights = w,
                  lambda = LAMBDA, standardize = TRUE),
           newx = X_test, s = LAMBDA, type = "response")),
         status = "ok"),
    error = function(e) list(pred = NULL, status = "error"))
}

subject_auroc <- function(case_ids, pred_probs, y_subject) {
  sp <- tapply(pred_probs, case_ids, mean)
  st <- y_subject[as.integer(names(sp))]
  if (length(unique(st)) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(st, sp, quiet = TRUE)))
}

# Stratified assignment of items to k folds within each outcome class.
strat_folds <- function(class_labels, k, seed) {
  set.seed(seed)
  f <- integer(length(class_labels))
  for (cl in unique(class_labels)) {
    who <- which(class_labels == cl)
    f[who] <- sample(rep(seq_len(k), length.out = length(who)))
  }
  f
}

# naive: observations stratified by their subject-level outcome.
naive_kfold <- function(dat, k, seed) {
  fold <- strat_folds(dat$y_obs, k, seed)
  preds <- rep(NA_real_, nrow(dat$X)); ok <- TRUE
  for (fld in seq_len(k)) {
    tr <- which(fold != fld); te <- which(fold == fld)
    if (length(unique(dat$y_obs[tr])) < 2) { ok <- FALSE; next }
    fp <- fit_predict_ridge(dat$X[tr, , drop = FALSE], dat$y_obs[tr],
                            dat$X[te, , drop = FALSE])
    if (fp$status != "ok") { ok <- FALSE; next }
    preds[te] <- fp$pred
  }
  if (any(is.na(preds))) ok <- FALSE
  list(auc = if (all(!is.na(preds))) subject_auroc(dat$case_ids, preds, dat$y_subject) else NA_real_,
       valid = ok)
}

# subject-level: subjects stratified by subject-level outcome, then mapped to obs.
cluster_kfold <- function(dat, k, seed) {
  uid <- sort(unique(dat$case_ids))
  y_subj <- dat$y_subject[uid]
  actual_k <- min(k, length(uid))
  sfold <- strat_folds(y_subj, actual_k, seed)
  fold_map <- setNames(sfold, uid)
  obs_fold <- fold_map[as.character(dat$case_ids)]
  preds <- rep(NA_real_, nrow(dat$X)); ok <- TRUE
  for (fld in seq_len(actual_k)) {
    tr <- which(obs_fold != fld); te <- which(obs_fold == fld)
    if (length(te) == 0 || length(unique(dat$y_obs[tr])) < 2) { ok <- FALSE; next }
    fp <- fit_predict_ridge(dat$X[tr, , drop = FALSE], dat$y_obs[tr],
                            dat$X[te, , drop = FALSE])
    if (fp$status != "ok") { ok <- FALSE; next }
    preds[te] <- fp$pred
  }
  if (any(is.na(preds))) ok <- FALSE
  list(auc = if (all(!is.na(preds))) subject_auroc(dat$case_ids, preds, dat$y_subject) else NA_real_,
       valid = ok)
}

loco_cv <- function(dat) {
  uid <- sort(unique(dat$case_ids))
  sp  <- numeric(length(uid)); names(sp) <- uid; ok <- TRUE
  for (i in seq_along(uid)) {
    te <- which(dat$case_ids == uid[i]); tr <- which(dat$case_ids != uid[i])
    if (length(unique(dat$y_obs[tr])) < 2) { sp[i] <- 0.5; next }
    fp <- fit_predict_ridge(dat$X[tr, , drop = FALSE], dat$y_obs[tr],
                            dat$X[te, , drop = FALSE])
    if (fp$status != "ok") { ok <- FALSE; sp[i] <- NA_real_; next }
    sp[i] <- mean(fp$pred)
  }
  st <- dat$y_subject[as.integer(names(sp))]
  if (!ok || any(is.na(sp)) || length(unique(st)) < 2)
    return(list(auc = NA_real_, valid = FALSE))
  list(auc = as.numeric(pROC::auc(pROC::roc(st, sp, quiet = TRUE))), valid = TRUE)
}

# Average over N_CV_REPS replicates; require ALL valid, else flag the dataset.
cv_mean <- function(dat, fn, k, n_reps) {
  a <- numeric(n_reps); v <- logical(n_reps)
  for (r in seq_len(n_reps)) { out <- fn(dat, k, r); a[r] <- out$auc; v[r] <- out$valid }
  list(mean = if (all(v)) mean(a) else NA_real_, sd = if (all(v)) sd(a) else NA_real_,
       n_valid = sum(v))
}

# ── 8. SINGLE-DATASET ANALYSIS ───────────────────────────────────────────────
# Returns one row of per-dataset results, or a flag row if fitting/folds failed.

analyze_dataset <- function(dat, params) {
  chars <- compute_empirical_chars(dat)
  naive <- cv_mean(dat, naive_kfold,   10, N_CV_REPS)
  s10   <- cv_mean(dat, cluster_kfold, 10, N_CV_REPS)
  s5    <- cv_mean(dat, cluster_kfold,  5, N_CV_REPS)
  loco  <- loco_cv(dat)

  cf <- fit_ridge_full_coef(dat$X, dat$y_obs, LAMBDA)
  fit_ok <- cf$status == "ok"
  true_auroc <- if (fit_ok)
    compute_true_auroc(dat$struct, cf$coef, N_TEST, FIXED_N_DAYS,
                       dat$seed_used + 500000L) else NA_real_

  folds_ok <- naive$n_valid == N_CV_REPS && s10$n_valid == N_CV_REPS &&
              s5$n_valid == N_CV_REPS && loco$valid
  flagged  <- !(fit_ok && folds_ok && !is.na(true_auroc))

  data.frame(
    true_N = params$N, true_ICC = params$ICC, true_AR1 = params$AR1,
    true_event_rate = params$event_rate, signal_mode = params$signal_mode,
    n_events = dat$n_events, actual_event_rate = dat$n_events / dat$n_subjects,
    regen = dat$regen,
    empirical_icc = chars$empirical_icc, empirical_ar1 = chars$empirical_ar1,
    naive_mean = naive$mean, naive_sd = naive$sd,
    s10_mean = s10$mean, s10_sd = s10$sd, s5_mean = s5$mean, loco_auroc = loco$auc,
    true_auroc = true_auroc,
    # optimism / residual discrepancy vs true (primary) and vs LOCO (retained)
    optimism_naive_true = naive$mean - true_auroc,
    resid_s10_true = s10$mean  - true_auroc,   # residual internal-validation discrepancy
    resid_s5_true  = s5$mean   - true_auroc,
    resid_loco_true = loco$auc - true_auroc,
    optimism_naive_loco = naive$mean - loco$auc,
    dev_s10_loco = s10$mean - loco$auc,
    # diagnostics
    fit_status = cf$status, folds_ok = folds_ok, flagged = flagged,
    naive_reps_valid = naive$n_valid, s10_reps_valid = s10$n_valid,
    s5_reps_valid = s5$n_valid, loco_valid = loco$valid,
    stringsAsFactors = FALSE)
}

# ── P. PILOT: freeze |beta| and smoke-test every stage ───────────────────────

oracle_auroc_ref <- function(bmag, n = 3e6) {
  v  <- var_eta(ICC_REF, bmag)
  a0 <- calibrate_intercept(PREV_REF, v)
  set.seed(MASTER_SEED)
  eta <- rnorm(n, 0, sqrt(v)); y <- as.integer(runif(n) < plogis(a0 + eta))
  as.numeric(pROC::auc(pROC::roc(y, eta, quiet = TRUE,
                                 direction = "<", levels = c(0, 1))))
}

run_pilot <- function() {
  cat(sprintf("Fixed R: max|off-diagonal| = %.3f ; QFORM_ACTIVE = %.4f\n",
              MAXABS_OFF, QFORM_ACTIVE))
  bmag <- uniroot(function(b) oracle_auroc_ref(b) - AUROC_TARGET,
                  lower = 0.05, upper = 5, tol = 1e-4)$root
  got <- oracle_auroc_ref(bmag)
  cat(sprintf("Frozen |beta| = %.4f  ->  reference oracle AUROC = %.4f (target %.2f)\n",
              bmag, got, AUROC_TARGET))
  if (abs(got - AUROC_TARGET) > AUROC_TOL)
    warning("Reference oracle AUROC outside tolerance; inspect before proceeding.")

  cat("\nBaseline oracle AUROC across ICC (expected to rise; confound acknowledged):\n")
  for (icc in c(0.1, 0.3, 0.5)) {
    v <- var_eta(icc, bmag); a0 <- calibrate_intercept(PREV_REF, v)
    set.seed(MASTER_SEED)
    eta <- rnorm(2e6, 0, sqrt(v)); y <- as.integer(runif(2e6) < plogis(a0 + eta))
    cat(sprintf("  ICC=%.1f  Var_eta=%.4f  oracle AUROC=%.4f\n", icc, v,
                as.numeric(pROC::auc(pROC::roc(y, eta, quiet = TRUE,
                                               direction = "<", levels = c(0, 1))))))
  }

  saveRDS(list(bmag_ref = bmag, R_PRED = R_PRED, QFORM_ACTIVE = QFORM_ACTIVE,
               ICC_REF = ICC_REF, PREV_REF = PREV_REF, AUROC_TARGET = AUROC_TARGET,
               reference_oracle_auroc = got, timestamp = Sys.time()),
          file.path(OUT_DIR, "frozen_magnitude.rds"))

  cat("\nEnd-to-end smoke test on 2 cells (small N_SIM, small N_TEST):\n")
  old_ntest <- N_TEST; assign("N_TEST", 300, envir = .GlobalEnv)
  smoke <- expand.grid(N = c(20, 100), ICC = 0.3, AR1 = 0.5, event_rate = 0.20)
  for (i in seq_len(nrow(smoke))) {
    p <- smoke[i, ]; p$signal_mode <- "baseline"
    beta <- make_true_beta(bmag, p$ICC, "baseline")
    a0   <- calibrate_intercept(p$event_rate, var_eta(p$ICC, bmag))
    rows <- lapply(1:5, function(s) {
      dat <- generate_dataset(p$N, FIXED_N_DAYS, FIXED_N_PRED, p$ICC, p$AR1,
                              beta, a0, seed = 999000L + s)
      analyze_dataset(dat, p)
    })
    df <- do.call(rbind, rows)
    cat(sprintf("  N=%3d: true=%.3f naive=%.3f s10=%.3f loco=%.3f | optT=%+.3f residT=%+.3f | flagged=%d regen(mean)=%.1f\n",
                p$N, mean(df$true_auroc), mean(df$naive_mean), mean(df$s10_mean),
                mean(df$loco_auroc), mean(df$optimism_naive_true),
                mean(df$resid_s10_true), sum(df$flagged), mean(df$regen)))
  }
  assign("N_TEST", old_ntest, envir = .GlobalEnv)
  cat("\nPilot complete. frozen_magnitude.rds written. Set RUN_MODE and re-source.\n")
}

# ── 9. LOAD FROZEN MAGNITUDE (main / sensitivity) ────────────────────────────
load_frozen <- function() {
  f <- file.path(OUT_DIR, "frozen_magnitude.rds")
  if (!file.exists(f)) stop("Run RUN_MODE='pilot' first; frozen_magnitude.rds not found.")
  readRDS(f)$bmag_ref
}

# ── 10. MAIN FACTORIAL (baseline arm) ────────────────────────────────────────

param_grid_main <- expand.grid(
  N = c(20, 30, 50, 75, 100, 150), ICC = c(0.1, 0.3, 0.5),
  AR1 = c(0.2, 0.5, 0.8), event_rate = c(0.10, 0.20, 0.35),
  stringsAsFactors = FALSE)

cond_key_of <- function(p, mode) sprintf("%s_N%d_ICC%.1f_AR%.1f_ER%.2f",
                                         mode, p$N, p$ICC, p$AR1, p$event_rate)

run_condition_baseline <- function(p, bmag) {
  beta <- make_true_beta(bmag, p$ICC, "baseline")
  a0   <- calibrate_intercept(p$event_rate, var_eta(p$ICC, bmag))
  p$signal_mode <- "baseline"
  base_seed <- p$N * 100000L + as.integer(p$ICC * 10) * 10000L +
               as.integer(p$AR1 * 10) * 1000L + as.integer(p$event_rate * 100) * 10L
  rows <- vector("list", N_SIM)
  regen_total <- 0L; flagged_total <- 0L; fit_err <- 0L
  for (s in seq_len(N_SIM)) {
    dat <- generate_dataset(p$N, FIXED_N_DAYS, FIXED_N_PRED, p$ICC, p$AR1,
                            beta, a0, seed = base_seed + s)
    r <- analyze_dataset(dat, p)
    regen_total <- regen_total + r$regen
    flagged_total <- flagged_total + r$flagged
    if (r$fit_status != "ok") fit_err <- fit_err + 1L
    rows[[s]] <- r
  }
  list(df = do.call(rbind, rows),
       log = data.frame(cond_key = cond_key_of(p, "baseline"),
                        n_valid = N_SIM, regen_total = regen_total,
                        flagged = flagged_total, fit_errors = fit_err,
                        stringsAsFactors = FALSE))
}

# ── 11. SENSITIVITY ARM (paired baseline + controlled by common random numbers)
# For each (cell, sim) a single set of latent draws and outcome uniforms is used
# to produce BOTH arms: only the coefficient scale and the (precomputed) intercept
# differ. Regeneration advances BOTH arms together so pairing is preserved.

param_grid_sens <- expand.grid(
  N = c(20, 30, 50, 75, 100, 150), ICC = c(0.1, 0.3, 0.5),
  AR1 = AR1_REF, event_rate = PREV_REF, stringsAsFactors = FALSE)
if (SENS_INCLUDE_AR1_SLICE) {
  slice <- expand.grid(N = SENS_SLICE_N, ICC = c(0.1, 0.3, 0.5),
                       AR1 = c(0.2, 0.8), event_rate = PREV_REF,
                       stringsAsFactors = FALSE)
  param_grid_sens <- unique(rbind(param_grid_sens, slice))
}

generate_pair <- function(p, bmag, seed) {
  beta_b <- make_true_beta(bmag, p$ICC, "baseline")
  beta_c <- make_true_beta(bmag, p$ICC, "controlled")
  a0_b <- calibrate_intercept(p$event_rate, var_eta(p$ICC, bmag))
  a0_c <- calibrate_intercept(p$event_rate, var_eta(ICC_REF, bmag))  # constant Var_eta
  regen <- 0L
  repeat {
    sd <- seed + regen * 100003L
    db <- generate_core(p$N, FIXED_N_DAYS, FIXED_N_PRED, p$ICC, p$AR1, beta_b, a0_b, sd)
    dc <- generate_core(p$N, FIXED_N_DAYS, FIXED_N_PRED, p$ICC, p$AR1, beta_c, a0_c, sd)
    okb <- db$n_events >= 3 && db$n_events <= p$N - 3
    okc <- dc$n_events >= 3 && dc$n_events <= p$N - 3
    if (okb && okc) { db$regen <- regen; dc$regen <- regen; return(list(baseline = db, controlled = dc)) }
    regen <- regen + 1L
    if (regen > MAX_REGEN) stop("MAX_REGEN exceeded in paired sensitivity generation.")
  }
}

run_condition_sens <- function(p, bmag) {
  base_seed <- 700000000L + p$N * 100000L + as.integer(p$ICC * 10) * 10000L +
               as.integer(p$AR1 * 10) * 1000L + as.integer(p$event_rate * 100) * 10L
  rows <- vector("list", N_SIM)
  for (s in seq_len(N_SIM)) {
    pr <- generate_pair(p, bmag, base_seed + s)
    pb <- p; pb$signal_mode <- "baseline"
    pc <- p; pc$signal_mode <- "controlled"
    rows[[s]] <- rbind(analyze_dataset(pr$baseline, pb),
                       analyze_dataset(pr$controlled, pc))
  }
  do.call(rbind, rows)
}

# ── 12. DRIVER: dispatch on RUN_MODE, with checkpoint / RESUME ────────────────

run_factorial <- function(grid, runner, tag, bmag) {
  res_file  <- file.path(OUT_DIR, sprintf("sim_results_%s.csv", tag))
  log_file  <- file.path(OUT_DIR, sprintf("sim_cond_log_%s.csv", tag))
  prog_file <- file.path(OUT_DIR, sprintf("sim_progress_%s.rds", tag))

  if (RESUME && file.exists(res_file)) {
    all_results <- read.csv(res_file, stringsAsFactors = FALSE)
    completed <- if (file.exists(prog_file)) readRDS(prog_file) else character(0)
    cond_log <- if (file.exists(log_file)) read.csv(log_file, stringsAsFactors = FALSE) else NULL
    cat(sprintf("Resuming %s: %d rows, %d conditions done\n", tag, nrow(all_results), length(completed)))
  } else {
    all_results <- data.frame(); completed <- character(0); cond_log <- NULL
  }

  t0 <- Sys.time()
  for (i in seq_len(nrow(grid))) {
    p <- grid[i, ]; ck <- cond_key_of(p, tag)
    if (ck %in% completed) { cat(sprintf("[%d/%d] %s SKIPPED\n", i, nrow(grid), ck)); next }
    tc <- Sys.time()
    out <- runner(p, bmag)
    df  <- if (is.list(out) && !is.null(out$df)) out$df else out
    all_results <- rbind(all_results, df)
    if (is.list(out) && !is.null(out$log)) cond_log <- rbind(cond_log, out$log)
    completed <- c(completed, ck)
    cat(sprintf("[%d/%d] %s  %d rows  %.0fs\n", i, nrow(grid), ck, nrow(df),
                as.numeric(difftime(Sys.time(), tc, units = "secs"))))
    if (TRUE) {  # checkpoint after every condition (was: i %% 5 == 0 || i == nrow(grid))
      write.csv(all_results, res_file, row.names = FALSE)
      if (!is.null(cond_log)) write.csv(cond_log, log_file, row.names = FALSE)
      saveRDS(completed, prog_file)
      cat(sprintf("  >> checkpoint: %d rows saved\n", nrow(all_results)))
    }
  }
  write.csv(all_results, res_file, row.names = FALSE)
  if (!is.null(cond_log)) write.csv(cond_log, log_file, row.names = FALSE)
  saveRDS(completed, prog_file)
  cat(sprintf("%s complete: %d rows in %.1f min\n", tag, nrow(all_results),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  list(results = all_results, res_file = res_file, log_file = log_file)
}

# ── 13. SUMMARIES (condition-weighted primary; scheme-invariance; SD ratio) ──

summarize_run <- function(all_results, tag) {
  keys <- c("true_N", "true_ICC", "true_AR1", "true_event_rate", "signal_mode")
  cond <- aggregate(cbind(true_auroc, naive_mean, s10_mean, s5_mean, loco_auroc,
                          optimism_naive_true, resid_s10_true, resid_s5_true,
                          resid_loco_true, optimism_naive_loco, naive_sd, s10_sd) ~
                    true_N + true_ICC + true_AR1 + true_event_rate + signal_mode,
                    data = all_results, FUN = function(x) mean(x, na.rm = TRUE))
  # scheme-invariance of the residual discrepancy (absolute differences)
  cond$absdiff_s5_s10   <- abs(cond$resid_s5_true  - cond$resid_s10_true)
  cond$absdiff_loco_s10 <- abs(cond$resid_loco_true - cond$resid_s10_true)
  cond$precision_ratio  <- cond$naive_sd / cond$s10_sd    # precision illusion

  write.csv(cond, file.path(OUT_DIR, sprintf("sim_condsummary_%s.csv", tag)),
            row.names = FALSE)

  eqw <- function(v, mode) mean(cond[cond$signal_mode == mode, v], na.rm = TRUE)
  cat(sprintf("\n== HEADLINE (condition-weighted equal mean), tag=%s ==\n", tag))
  for (mode in unique(cond$signal_mode)) {
    cat(sprintf(" [%s] total opt(naive-true)=%+.4f  partition(naive-s10)=%+.4f  resid(s10-true)=%+.4f\n",
                mode, eqw("optimism_naive_true", mode),
                mean(cond[cond$signal_mode == mode, "naive_mean"] -
                     cond[cond$signal_mode == mode, "s10_mean"], na.rm = TRUE),
                eqw("resid_s10_true", mode)))
    cat(sprintf("      scheme-invariance: mean|s5-s10|=%.4f max|s5-s10|=%.4f  mean|loco-s10|=%.4f max|loco-s10|=%.4f\n",
                mean(cond$absdiff_s5_s10[cond$signal_mode == mode], na.rm = TRUE),
                max(cond$absdiff_s5_s10[cond$signal_mode == mode], na.rm = TRUE),
                mean(cond$absdiff_loco_s10[cond$signal_mode == mode], na.rm = TRUE),
                max(cond$absdiff_loco_s10[cond$signal_mode == mode], na.rm = TRUE)))
  }
  # residual by N and by ICC (for the reframed "governed by N" claim)
  cat(" residual(s10-true) by N (baseline):\n")
  bl <- cond[cond$signal_mode == "baseline", ]
  print(round(tapply(bl$resid_s10_true, bl$true_N, mean, na.rm = TRUE), 4))
  cat(" residual(s10-true) by ICC (baseline vs controlled if present):\n")
  print(round(tapply(cond$resid_s10_true, list(cond$true_ICC, cond$signal_mode), mean, na.rm = TRUE), 4))
  cond
}

# ── 14. REPRODUCIBILITY ARTIFACT ─────────────────────────────────────────────

write_repro_artifact <- function(files, tag, bmag) {
  meta <- list(
    tag = tag, timestamp = Sys.time(), master_seed = MASTER_SEED, R_seed = R_SEED,
    frozen_beta_magnitude = bmag, active_signs = ACTIVE_SIGNS,
    icc_ref = ICC_REF, prev_ref = PREV_REF, auroc_target = AUROC_TARGET,
    N_SIM = N_SIM, N_CV_REPS = N_CV_REPS, N_TEST = N_TEST, LAMBDA = LAMBDA,
    max_abs_off_diagonal_R = MAXABS_OFF, qform_active = QFORM_ACTIVE,
    glmnet = as.character(packageVersion("glmnet")),
    pROC = as.character(packageVersion("pROC")),
    R_version = R.version.string)
  writeLines(vapply(names(meta), function(k)
    sprintf("%s: %s", k, paste(meta[[k]], collapse = ", ")), character(1)),
    file.path(OUT_DIR, sprintf("run_metadata_%s.txt", tag)))
  present <- files[file.exists(files)]
  if (length(present))
    write.csv(data.frame(file = basename(present),
                         md5 = tools::md5sum(present), row.names = NULL),
              file.path(OUT_DIR, sprintf("checksums_%s.csv", tag)), row.names = FALSE)
  writeLines(capture.output(sessionInfo()),
             file.path(OUT_DIR, sprintf("sessionInfo_%s.txt", tag)))
  cat(sprintf("Reproducibility artifact written for tag=%s\n", tag))
}

# ── 15. EXECUTE ──────────────────────────────────────────────────────────────

if (RUN_MODE == "pilot") {
  run_pilot()
} else if (RUN_MODE == "main") {
  bmag <- load_frozen()
  cat(sprintf("MAIN factorial (baseline): |beta|=%.4f, %d conditions\n",
              bmag, nrow(param_grid_main)))
  out  <- run_factorial(param_grid_main, run_condition_baseline, "main", bmag)
  cond <- summarize_run(out$results, "main")
  write_repro_artifact(c(out$res_file, out$log_file,
                         file.path(OUT_DIR, "sim_condsummary_main.csv")), "main", bmag)
} else if (RUN_MODE == "sensitivity") {
  bmag <- load_frozen()
  cat(sprintf("SENSITIVITY (paired baseline+controlled): |beta|=%.4f, %d cells\n",
              bmag, nrow(param_grid_sens)))
  out  <- run_factorial(param_grid_sens, function(p, b) run_condition_sens(p, b),
                        "sensitivity", bmag)
  cond <- summarize_run(out$results, "sensitivity")
  write_repro_artifact(c(out$res_file,
                         file.path(OUT_DIR, "sim_condsummary_sensitivity.csv")),
                       "sensitivity", bmag)
} else {
  stop("RUN_MODE must be 'pilot', 'main', or 'sensitivity'.")
}