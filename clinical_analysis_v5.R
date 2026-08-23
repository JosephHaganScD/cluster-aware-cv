###############################################################################
# Clinical Illustration: Cluster-Aware Cross-Validation for ROP Prediction
# =========================================================================
# UNIFIED PIPELINE (v5) — supersedes clinical_analysis.R (April 2026) and
# clinical_addons.R (August 2026), which are merged here.
#
# Why this file exists:
#   The April 2026 clinical_analysis.R still reported DeLong CIs for LOCO and
#   percentile intervals on the clinical optimism, both of which were changed
#   in the manuscript (tracker T2, T3) but never propagated to the code. The
#   bootstrap CI was instead computed in the standalone clinical_addons.R,
#   which used a DIFFERENT standardization convention (manual, unweighted,
#   1/(n-1) denominator, standardize = FALSE) from the main script (glmnet
#   internal, weighted, 1/n denominator, standardize = TRUE). Ridge penalties
#   are scale-dependent, so the two conventions fit slightly different models:
#   that is the source of the severe-ROP LOCO 0.606 (add-on) vs 0.608 (main)
#   discrepancy.
#
#   This script computes EVERYTHING from one convention and one set of fits:
#   Tables 1-2 (AUROC, optimism, intervals), Table 3 (fold-wise coefficient
#   stability), Section 3.3 (precision illusion), Section 3.4 (summary models),
#   Section 3.5 (regularization sensitivity), and the Appendix dependence table.
#
# Convention adopted: glmnet(standardize = TRUE).
#   Chosen because it produced every AUROC currently reported in the
#   manuscript, so no reported discrimination value changes. Standardized-scale
#   coefficients for Table 3 are recovered analytically from the same fits
#   (see standardized_coefs() and the verification block), rather than by
#   refitting under a second convention.
#
# Patch 2026-08-12: standardized_coefs() intercept correction (see comment at
#   the function). Slopes, and therefore every value in manuscript Table 3,
#   were unaffected; Table 3 does not report the intercept.
#
# Patch 2026-08-15 (reviewer verification items; no reported value changes):
#   (a) Observations are now explicitly ordered by Case_ID and then DOL after
#       the complete-case restriction. The supplied final_daily_data.csv is
#       already in this order and no subject has a non-monotone DOL sequence,
#       so every result is numerically unchanged; the sort makes the lag-1
#       autocorrelations in Appendix B correct by construction rather than by
#       reliance on the input file order.
#   (b) The one-way ICC is now reported under BOTH the arithmetic mean cluster
#       size and the effective cluster size for an unbalanced design,
#       n0_eff = (N - sum(n_i^2)/N)/(k - 1). Cluster sizes are nearly balanced
#       (44 to 60), so the two differ by 0.003 observations and every ICC by
#       less than 1e-4. Appendix B continues to report the arithmetic-mean
#       version, which is what manuscript Section 2.3 now states explicitly.
#
# Author: Joseph L. Hagan, ScD, MSPH
# Date:   August 2026 (patched 2026-08-12; verification patch 2026-08-15)
###############################################################################

library(glmnet)
library(pROC)

# ── USER SETTINGS ────────────────────────────────────────────────────────────

# DATA_FILE  <- file.path("data", "final_daily_data.csv")
# OUTPUT_DIR <- file.path("results", "rop_analysis")

DATA_FILE  <- "C:/Users/jlhagan/OneDrive - Texas Children's/WSS-onedrive-migration/Desktop/MY NOTES/cluster aware CV/3-13-2026 idea/3-13-26 claude material/final_daily_data.csv"
OUTPUT_DIR <- "C:/simulations/cluster aware cv/results/rop_analysis_v5"

N_CV_REPS    <- 50     # replicate fold assignments, Tables 1-2 and Section 3.3
N_STAB_REPS  <- 10     # replicate assignments for coefficient stability (Table 3)
N_BOOT       <- 2000   # subject-level bootstrap resamples for the LOCO CI
LAMBDA       <- 0.01
BOOT_SEED    <- 1

VERIFY_STANDARDIZATION <- TRUE  # run the one-off check that the analytic
                                # standardized coefficients match a refit

if (!file.exists(DATA_FILE)) {
  stop(sprintf("Data file not found: %s", DATA_FILE))
}
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

output_file <- file.path(
  OUTPUT_DIR, paste0("ROP_results_v5_", format(Sys.Date(), "%Y-%m-%d"), ".txt"))
sink(output_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat(sprintf("ROP Analysis Output (unified v5 pipeline) - %s\n", Sys.time()))
cat(sprintf("R %s; glmnet %s; pROC %s\n\n",
            getRversion(),
            as.character(packageVersion("glmnet")),
            as.character(packageVersion("pROC"))))

# ── 1. DATA PREPARATION ─────────────────────────────────────────────────────

raw <- read.csv(DATA_FILE)

pred_vars <- c("DOL", "Avg_FiO2", "Avg_SpO2", "Severe_Hypoxemia__",
               "Amb_Hyperox_", "Iatr_HyperOx_", "Swings", "Tit_Index")
oxy_vars  <- pred_vars[-1]

cc  <- complete.cases(raw[, pred_vars])
dat <- raw[cc, ]

# Explicit within-subject time ordering. Required for the lag-1 autocorrelations
# in Appendix B to be computed along the time axis; harmless if already sorted.
n_unsorted <- sum(tapply(dat$DOL, dat$Case_ID, function(z) is.unsorted(z)))
dat <- dat[order(dat$Case_ID, dat$DOL), ]
cat(sprintf("Complete cases: %d rows (%d excluded)\n", nrow(dat), sum(!cc)))
cat(sprintf("Subjects with non-monotone DOL in the input file order: %d\n",
            n_unsorted))
cl_sizes <- as.integer(table(dat$Case_ID))
cat(sprintf("Observations per subject: %d to %d (median %.0f)\n",
            min(cl_sizes), max(cl_sizes), median(cl_sizes)))

dat$severe_rop <- as.integer(dat$ROP == 2)
dat$any_rop    <- as.integer(dat$ROP >= 1)

for (nm in c("severe_rop", "any_rop")) {
  cat(sprintf("%-11s %d events among %d subjects (%.1f%%)\n", nm,
              sum(tapply(dat[[nm]], dat$Case_ID, max)),
              length(unique(dat$Case_ID)),
              100 * mean(tapply(dat[[nm]], dat$Case_ID, max))))
}

# ── 2. DEPENDENCE STRUCTURE (Appendices A and B) ────────────────────────────

# Cluster-size term for the one-way random-effects ICC. With unequal cluster
# sizes the standard unbalanced-design term is the effective cluster size
# n0_eff = (N - sum(n_i^2)/N)/(k - 1) rather than the arithmetic mean. Both are
# computed here so that the choice is documented rather than implicit.
n_i     <- as.integer(table(dat$Case_ID))
N_obs   <- sum(n_i); k_cl <- length(n_i)
n0      <- mean(n_i)                                        # arithmetic mean
n0_eff  <- (N_obs - sum(n_i^2) / N_obs) / (k_cl - 1)        # unbalanced design

cat(sprintf("\nCluster-size term: arithmetic mean n0 = %.4f; effective n0 = %.4f\n",
            n0, n0_eff))

compute_icc <- function(x, group, n0_use = n0) {
  ms   <- summary(aov(x ~ factor(group)))[[1]]
  ms_b <- ms["Mean Sq"][1, 1]; ms_w <- ms["Mean Sq"][2, 1]
  (ms_b - ms_w) / (ms_b + (n0_use - 1) * ms_w)
}

compute_lag1 <- function(var_name) {
  # dat is sorted by Case_ID then DOL above, so x is in time order within subject
  lag1 <- sapply(unique(dat$Case_ID), function(id) {
    sub <- dat[dat$Case_ID == id, ]
    x   <- sub[[var_name]][order(sub$DOL)]
    if (length(x) < 3) return(NA)
    cor(x[-length(x)], x[-1])
  })
  c(median = median(lag1, na.rm = TRUE),
    q25 = unname(quantile(lag1, .25, na.rm = TRUE)),
    q75 = unname(quantile(lag1, .75, na.rm = TRUE)))
}

dep_table <- data.frame(
  Predictor = oxy_vars,
  Mean = sapply(oxy_vars, function(v) mean(dat[[v]], na.rm = TRUE)),
  SD   = sapply(oxy_vars, function(v) sd(dat[[v]], na.rm = TRUE)),
  ICC  = sapply(oxy_vars, function(v) compute_icc(dat[[v]], dat$Case_ID)),
  stringsAsFactors = FALSE)
lag1_res <- t(sapply(oxy_vars, compute_lag1))
dep_table$Lag1_Median <- lag1_res[, "median"]
dep_table$Lag1_Q25    <- lag1_res[, "q25"]
dep_table$Lag1_Q75    <- lag1_res[, "q75"]

cat("\n===== APPENDIX B: DEPENDENCE STRUCTURE =====\n")
print(format(dep_table, digits = 3), row.names = FALSE)

# Sensitivity of the ICC to the cluster-size term (manuscript Section 2.3).
icc_eff  <- sapply(oxy_vars, function(v)
  compute_icc(dat[[v]], dat$Case_ID, n0_use = n0_eff))
icc_cmp  <- data.frame(Predictor = oxy_vars,
                       ICC_mean_n0 = dep_table$ICC,
                       ICC_eff_n0  = icc_eff,
                       Difference  = icc_eff - dep_table$ICC,
                       row.names = NULL)
cat("\n----- ICC under arithmetic-mean vs effective cluster size -----\n")
print(format(icc_cmp, digits = 5), row.names = FALSE)
cat(sprintf("Maximum absolute difference in ICC = %.2e\n",
            max(abs(icc_cmp$Difference))))
write.csv(icc_cmp, file.path(OUTPUT_DIR, "appendixB_icc_n0_sensitivity.csv"),
          row.names = FALSE)
write.csv(dep_table, file.path(OUTPUT_DIR, "appendixB_dependence.csv"),
          row.names = FALSE)

cor_mat <- cor(dat[, oxy_vars], use = "complete.obs")
cat("\n===== APPENDIX A: PREDICTOR CORRELATION MATRIX =====\n")
print(round(cor_mat, 2))
cat(sprintf("\nMaximum |off-diagonal r| = %.2f\n",
            max(abs(cor_mat[lower.tri(cor_mat)]))))
write.csv(round(cor_mat, 3), file.path(OUTPUT_DIR, "appendixA_correlations.csv"))

# ── 3. SINGLE FITTING CONVENTION ────────────────────────────────────────────
# One function used by every CV strategy, the summary models, the regularization
# sweep, and the coefficient-stability analysis.

.fallback_count <- 0L   # counts silent lambda fallbacks; must be 0 at the end

fit_ridge <- function(X_train, y_train, lambda = LAMBDA) {
  n_pos <- sum(y_train == 1); n_neg <- sum(y_train == 0)
  n_tot <- length(y_train)
  w <- ifelse(y_train == 1, n_tot / (2 * n_pos), n_tot / (2 * n_neg))
  fit <- tryCatch(
    glmnet(X_train, y_train, family = "binomial", alpha = 0,
           weights = w, lambda = lambda, standardize = TRUE),
    error = function(e) {
      .fallback_count <<- .fallback_count + 1L
      glmnet(X_train, y_train, family = "binomial", alpha = 0,
             weights = w, lambda = 0.1, standardize = TRUE)
    })
  attr(fit, "weights") <- w
  fit
}

predict_ridge <- function(fit, X_test) {
  as.numeric(predict(fit, newx = X_test, s = fit$lambda[1], type = "response"))
}

fit_predict_ridge <- function(X_train, y_train, X_test, lambda = LAMBDA) {
  predict_ridge(fit_ridge(X_train, y_train, lambda), X_test)
}

# Standardized-scale coefficients recovered from a standardize = TRUE fit.
# glmnet standardizes x using the WEIGHTED mean and the weighted second moment
# with a 1/n (not 1/(n-1)) denominator, fits on the standardized scale, and
# returns beta on the original scale. Therefore beta_std_j = beta_raw_j * xs_j.
# The intercept requires its own correction: on the standardized scale
#   beta0_std = beta0_raw + sum_j(beta_raw_j * xm_j),
# because centring x shifts the intercept by the centred contribution of every
# slope. Omitting this term leaves the slopes correct but the intercept wrong,
# which is what the [CHECK] block below caught on 2026-08-12.
standardized_coefs <- function(fit, X_train, w) {
  wn <- w / sum(w)
  xm <- as.numeric(crossprod(wn, X_train))
  xs <- sqrt(as.numeric(crossprod(wn, sweep(X_train, 2, xm, "-")^2)))
  xs[xs < 1e-8] <- 1e-8
  b  <- as.numeric(coef(fit, s = fit$lambda[1]))   # intercept first
  c(b[1] + sum(b[-1] * xm), b[-1] * xs)
}

# ── 4. CROSS-VALIDATION STRATEGIES ──────────────────────────────────────────

subject_auroc <- function(case_ids, pred_probs, y_outcome) {
  subj_preds <- tapply(pred_probs, case_ids, mean)
  subj_true  <- tapply(y_outcome, case_ids, function(x) x[1])[names(subj_preds)]
  if (length(unique(subj_true)) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(as.numeric(subj_true),
                                 as.numeric(subj_preds), quiet = TRUE)))
}

naive_kfold <- function(X, y_obs, case_ids, k = 10, seed = 0, lambda = LAMBDA) {
  set.seed(seed)
  n <- nrow(X)
  fold_ids <- sample(rep(1:k, length.out = n))
  preds <- rep(NA_real_, n)
  for (fold in 1:k) {
    te <- which(fold_ids == fold); tr <- which(fold_ids != fold)
    if (length(unique(y_obs[tr])) < 2) next
    preds[te] <- fit_predict_ridge(X[tr, , drop = FALSE], y_obs[tr],
                                   X[te, , drop = FALSE], lambda)
  }
  v <- !is.na(preds)
  subject_auroc(case_ids[v], preds[v], y_obs[v])
}

cluster_kfold <- function(X, y_obs, case_ids, k = 10, seed = 0, lambda = LAMBDA) {
  set.seed(seed)
  uid <- unique(case_ids); n_subj <- length(uid)
  kk  <- min(k, n_subj)
  fold_map  <- setNames(rep(1:kk, length.out = n_subj), uid[sample(n_subj)])
  obs_folds <- fold_map[as.character(case_ids)]
  preds <- rep(NA_real_, nrow(X))
  for (fold in 1:kk) {
    te <- which(obs_folds == fold); tr <- which(obs_folds != fold)
    if (length(te) == 0 || length(tr) == 0) next
    if (length(unique(y_obs[tr])) < 2) next
    preds[te] <- fit_predict_ridge(X[tr, , drop = FALSE], y_obs[tr],
                                   X[te, , drop = FALSE], lambda)
  }
  v <- !is.na(preds)
  subject_auroc(case_ids[v], preds[v], y_obs[v])
}

# LOCO returns the subject-level prediction vector itself, so that the point
# estimate, the bootstrap CI, and the DeLong CI all derive from ONE object.
loco_predictions <- function(X, y_obs, case_ids, lambda = LAMBDA) {
  uid <- unique(case_ids)
  subj_pred <- numeric(length(uid)); names(subj_pred) <- uid
  subj_true <- numeric(length(uid)); names(subj_true) <- uid
  for (i in seq_along(uid)) {
    te <- which(case_ids == uid[i]); tr <- which(case_ids != uid[i])
    subj_true[i] <- y_obs[te][1]
    if (length(unique(y_obs[tr])) < 2) { subj_pred[i] <- 0.5; next }
    subj_pred[i] <- mean(fit_predict_ridge(X[tr, , drop = FALSE], y_obs[tr],
                                           X[te, , drop = FALSE], lambda))
  }
  list(pred = subj_pred, true = subj_true)
}

loco_summary <- function(lp, n_boot = N_BOOT, seed = BOOT_SEED) {
  roc0  <- pROC::roc(lp$true, lp$pred, quiet = TRUE)
  point <- as.numeric(pROC::auc(roc0))
  set.seed(seed)
  n <- length(lp$true); boot <- numeric(0)
  for (b in 1:n_boot) {
    idx <- sample.int(n, n, replace = TRUE)
    if (length(unique(lp$true[idx])) < 2) next
    boot <- c(boot, as.numeric(pROC::auc(
      pROC::roc(lp$true[idx], lp$pred[idx], quiet = TRUE))))
  }
  ci     <- as.numeric(quantile(boot, c(.025, .975), na.rm = TRUE))
  delong <- as.numeric(pROC::ci.auc(roc0, method = "delong"))
  list(auc = point, boot_lo = ci[1], boot_hi = ci[2],
       n_boot_valid = length(boot), delong_lo = delong[1], delong_hi = delong[3])
}

# ── 5. FOLD-WISE COEFFICIENT STABILITY (Table 3) ────────────────────────────

coef_stability <- function(dat, outcome_var, pred_vars, k = 10,
                           n_reps = N_STAB_REPS) {
  X <- as.matrix(dat[, pred_vars]); y <- dat[[outcome_var]]
  case_ids <- dat$Case_ID; uid <- unique(case_ids); n_subj <- length(uid)
  collected <- list()
  for (r in 1:n_reps) {
    set.seed(r)
    kk <- min(k, n_subj)
    fold_map  <- setNames(rep(1:kk, length.out = n_subj), uid[sample(n_subj)])
    obs_folds <- fold_map[as.character(case_ids)]
    for (fold in 1:kk) {
      tr <- which(obs_folds != fold)
      if (length(unique(y[tr])) < 2) next
      Xtr <- X[tr, , drop = FALSE]
      fit <- fit_ridge(Xtr, y[tr])
      collected[[length(collected) + 1]] <-
        standardized_coefs(fit, Xtr, attr(fit, "weights"))
    }
  }
  M <- do.call(rbind, collected)
  colnames(M) <- c("(Intercept)", pred_vars)
  data.frame(predictor = colnames(M),
             mean_coef = colMeans(M),
             sd_coef   = apply(M, 2, sd),
             n_folds   = nrow(M),
             row.names = NULL, stringsAsFactors = FALSE)
}

# One-off numerical check that standardized_coefs() inverts glmnet correctly:
# refit with standardize = FALSE on manually weighted-standardized X and
# compare. Differences should be at solver tolerance (< 1e-6).
if (VERIFY_STANDARDIZATION) {
  Xv <- as.matrix(dat[, pred_vars]); yv <- dat$severe_rop
  f1 <- fit_ridge(Xv, yv); wv <- attr(f1, "weights")
  b_analytic <- standardized_coefs(f1, Xv, wv)
  wn <- wv / sum(wv)
  xm <- as.numeric(crossprod(wn, Xv))
  xs <- sqrt(as.numeric(crossprod(wn, sweep(Xv, 2, xm, "-")^2)))
  Xs <- scale(Xv, center = xm, scale = xs)
  f2 <- glmnet(Xs, yv, family = "binomial", alpha = 0, weights = wv,
               lambda = LAMBDA, standardize = FALSE)
  b_refit <- as.numeric(coef(f2, s = LAMBDA))
  cat(sprintf(
    "\n[CHECK] standardized-coefficient recovery: max abs difference = %.3e\n",
    max(abs(b_analytic - b_refit))))
  cat("        (should be < 1e-6; if not, report Table 3 from the refit path)\n")
}

# ── 6. SUBJECT-LEVEL SUMMARY MODELS (Section 3.4) ───────────────────────────

summary_cv <- function(dat, outcome_var, use_slopes = FALSE, k = 10,
                       n_reps = N_CV_REPS, pred_vars) {
  uid <- unique(dat$Case_ID); n_subj <- length(uid)
  subj_means <- do.call(rbind, lapply(uid, function(id)
    colMeans(dat[dat$Case_ID == id, pred_vars, drop = FALSE], na.rm = TRUE)))
  colnames(subj_means) <- paste0("mean_", pred_vars)
  features <- subj_means
  if (use_slopes) {
    subj_slopes <- do.call(rbind, lapply(uid, function(id) {
      sub <- dat[dat$Case_ID == id, ]
      sapply(pred_vars[-1], function(v) {
        if (sd(sub$DOL) == 0) return(0)
        unname(coef(lm(sub[[v]] ~ sub$DOL))[2])
      })
    }))
    colnames(subj_slopes) <- paste0("slope_", pred_vars[-1])
    features <- cbind(subj_means, subj_slopes)
  }
  subj_y <- sapply(uid, function(id) dat[[outcome_var]][dat$Case_ID == id][1])
  Xs <- as.matrix(features)

  kfold_aucs <- sapply(1:n_reps, function(r) {
    set.seed(r)
    fold_ids <- sample(rep(1:k, length.out = n_subj))
    preds <- rep(NA_real_, n_subj)
    for (fold in 1:k) {
      te <- which(fold_ids == fold); tr <- which(fold_ids != fold)
      if (length(unique(subj_y[tr])) < 2) next
      preds[te] <- fit_predict_ridge(Xs[tr, , drop = FALSE], subj_y[tr],
                                     Xs[te, , drop = FALSE])
    }
    v <- !is.na(preds)
    if (length(unique(subj_y[v])) < 2) return(NA_real_)
    as.numeric(pROC::auc(pROC::roc(subj_y[v], preds[v], quiet = TRUE)))
  })

  preds_loco <- rep(NA_real_, n_subj)
  for (i in 1:n_subj) {
    tr <- setdiff(1:n_subj, i)
    if (length(unique(subj_y[tr])) < 2) { preds_loco[i] <- 0.5; next }
    preds_loco[i] <- fit_predict_ridge(Xs[tr, , drop = FALSE], subj_y[tr],
                                       Xs[i, , drop = FALSE])
  }
  lsum <- loco_summary(list(pred = preds_loco, true = subj_y))

  list(kfold_mean = mean(kfold_aucs, na.rm = TRUE),
       kfold_sd   = sd(kfold_aucs, na.rm = TRUE),
       kfold_pi   = quantile(kfold_aucs, c(.025, .975), na.rm = TRUE),
       loco       = lsum$auc, loco_boot = c(lsum$boot_lo, lsum$boot_hi),
       n_features = ncol(Xs))
}

# ── 7. RUN MAIN ANALYSIS ────────────────────────────────────────────────────

run_analysis <- function(dat, outcome_var, outcome_label, pred_vars, n_reps) {
  cat(sprintf("\n\n========== %s ==========\n", outcome_label))
  X <- as.matrix(dat[, pred_vars]); y <- dat[[outcome_var]]
  case_ids <- dat$Case_ID

  cat("Naive 10-fold CV ...\n")
  naive_aucs  <- sapply(1:n_reps, function(r) naive_kfold(X, y, case_ids, 10, r))
  cat("Subject-level 10-fold CV ...\n")
  subj10_aucs <- sapply(1:n_reps, function(r) cluster_kfold(X, y, case_ids, 10, r))
  cat("LOCO CV + subject-level bootstrap CI ...\n")
  lp   <- loco_predictions(X, y, case_ids)
  lsum <- loco_summary(lp)

  naive_pi  <- quantile(naive_aucs,  c(.025, .975))
  subj10_pi <- quantile(subj10_aucs, c(.025, .975))

  cat("\n--- Observation-level results (Tables 1-2) ---\n")
  cat(sprintf("Naive 10-fold:    mean = %.3f  95%% fold-assignment interval (%.3f, %.3f)  SD = %.3f  optimism = %+.3f\n",
              mean(naive_aucs), naive_pi[1], naive_pi[2], sd(naive_aucs),
              mean(naive_aucs) - lsum$auc))
  cat(sprintf("Subject 10-fold:  mean = %.3f  95%% fold-assignment interval (%.3f, %.3f)  SD = %.3f  optimism = %+.3f\n",
              mean(subj10_aucs), subj10_pi[1], subj10_pi[2], sd(subj10_aucs),
              mean(subj10_aucs) - lsum$auc))
  cat(sprintf("LOCO (reference): %.4f  subject bootstrap 95%% CI (%.4f, %.4f) [%d valid resamples]\n",
              lsum$auc, lsum$boot_lo, lsum$boot_hi, lsum$n_boot_valid))
  cat(sprintf("                  DeLong 95%% CI (comparison only): (%.4f, %.4f)\n",
              lsum$delong_lo, lsum$delong_hi))

  cat(sprintf("\n--- Precision illusion (Section 3.3) ---\n"))
  cat(sprintf("SD naive = %.4f; SD subject-level = %.4f; ratio = %.3f (%.1f-fold)\n",
              sd(naive_aucs), sd(subj10_aucs),
              sd(naive_aucs) / sd(subj10_aucs),
              sd(subj10_aucs) / sd(naive_aucs)))

  cat("\nSummary models ...\n")
  sm <- summary_cv(dat, outcome_var, FALSE, 10, n_reps, pred_vars)
  ss <- summary_cv(dat, outcome_var, TRUE,  10, n_reps, pred_vars)

  cat("\n--- Subject-level summary models (Section 3.4) ---\n")
  for (z in list(list(sm, "means only"), list(ss, "means + slopes"))) {
    o <- z[[1]]; lab <- z[[2]]
    cat(sprintf("%-16s (%2d features): 10-fold %.3f (%.3f, %.3f), optimism %+.3f | LOCO %.3f, boot 95%% CI (%.3f, %.3f)\n",
                lab, o$n_features, o$kfold_mean, o$kfold_pi[1], o$kfold_pi[2],
                o$kfold_mean - lsum$auc, o$loco, o$loco_boot[1], o$loco_boot[2]))
  }

  list(naive_aucs = naive_aucs, subj10_aucs = subj10_aucs,
       loco = lsum, summ_means = sm, summ_slopes = ss)
}

res_severe <- run_analysis(dat, "severe_rop", "SEVERE ROP (15 events)",
                           pred_vars, N_CV_REPS)
res_any    <- run_analysis(dat, "any_rop", "ANY ROP (48 events)",
                           pred_vars, N_CV_REPS)

# ── 8. COEFFICIENT STABILITY (Table 3) ──────────────────────────────────────

for (oc in list(c("severe_rop", "Severe ROP"), c("any_rop", "Any ROP"))) {
  cat(sprintf("\n\n===== TABLE 3: fold-wise standardized coefficients, %s =====\n",
              oc[2]))
  cs <- coef_stability(dat, oc[1], pred_vars)
  print(format(cs, digits = 3), row.names = FALSE)
  write.csv(cs, file.path(OUTPUT_DIR, sprintf("coef_stability_%s.csv", oc[1])),
            row.names = FALSE)
}

# ── 9. REGULARIZATION SENSITIVITY (Section 3.5) ─────────────────────────────

cat("\n\n===== REGULARIZATION SENSITIVITY (Section 3.5) =====\n")
lambda_grid   <- c(1.0, 0.01, 0.00001)
lambda_labels <- c("strong (1.0)", "moderate (0.01)", "weak (0.00001)")

for (outcome_var in c("severe_rop", "any_rop")) {
  cat(sprintf("\n--- %s ---\n", outcome_var))
  X <- as.matrix(dat[, pred_vars]); y <- dat[[outcome_var]]
  case_ids <- dat$Case_ID
  for (j in seq_along(lambda_grid)) {
    lam <- lambda_grid[j]
    naive_lam <- sapply(1:N_CV_REPS, function(r)
      naive_kfold(X, y, case_ids, 10, r, lambda = lam))
    lp_lam <- loco_predictions(X, y, case_ids, lambda = lam)
    loco_lam <- as.numeric(pROC::auc(
      pROC::roc(lp_lam$true, lp_lam$pred, quiet = TRUE)))
    cat(sprintf("  lambda = %-16s naive = %.3f  LOCO = %.3f  optimism = %+.3f\n",
                lambda_labels[j], mean(naive_lam), loco_lam,
                mean(naive_lam) - loco_lam))
  }
}

# ── 10. SESSION AND INTEGRITY LOG ───────────────────────────────────────────

cat(sprintf("\n\n[INTEGRITY] silent lambda fallbacks during run: %d (must be 0)\n",
            .fallback_count))
cat("\n===== sessionInfo() =====\n")
print(sessionInfo())
cat(sprintf("\nResults written to: %s\n", output_file))
