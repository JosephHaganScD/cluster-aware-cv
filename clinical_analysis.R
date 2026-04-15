###############################################################################
# Clinical Illustration: Cluster-Aware Cross-Validation for ROP Prediction
# =========================================================================
# Applies the same ridge logistic regression + CV framework used in the
# simulation study to the motivating ROP dataset (Srivatsa et al., 2022).
#
# Produces Tables 1–2, Appendix dependence structure table, and
# supplementary results reported in Sections 3.1–3.4.
#
# Updates (April 2026):
#   - Appendix table: mean, SD, ICC, median lag-1 autocorrelation [IQR]
#     for all 7 oxygenation predictors + full 7x7 correlation matrix
#   - 95% percentile intervals for k-fold AUROC and mean optimism
#   - DeLong 95% CI for LOCO AUROC
#
# Author: Joseph L. Hagan, ScD, MSPH
# Date:   April 2026
###############################################################################

library(glmnet)
library(pROC)

# ── USER SETTINGS ────────────────────────────────────────────────────────────

# Expected repository layout:
#   data/final_daily_data.csv
#   results/rop_analysis/
DATA_FILE  <- file.path("data", "final_daily_data.csv")
OUTPUT_DIR <- file.path("results", "rop_analysis")

N_CV_REPS <- 50   # replicate fold assignments for k-fold strategies

if (!file.exists(DATA_FILE)) {
  stop(sprintf("Data file not found: %s\nPlace final_daily_data.csv in the data/ directory.", DATA_FILE))
}
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

output_file <- file.path(
  OUTPUT_DIR,
  paste0("ROP_results_", format(Sys.Date(), "%Y-%m-%d"), ".txt")
)
sink(output_file, split = TRUE)   # split = TRUE echoes to console too
on.exit(sink(), add = TRUE)
cat(sprintf("ROP Analysis Output — %s\n\n", Sys.time()))

rop <- read.csv(DATA_FILE)


# ── 1. DATA PREPARATION ─────────────────────────────────────────────────────

raw <- rop

# 8 predictors: DOL + 7 oxygenation variables
pred_vars <- c("DOL", "Avg_FiO2", "Avg_SpO2", "Severe_Hypoxemia__",
               "Amb_Hyperox_", "Iatr_HyperOx_", "Swings", "Tit_Index")

# 7 oxygenation predictors (excluding DOL) for dependence structure
oxy_vars <- pred_vars[-1]

# Complete-case restriction (drop rows with any missing predictor)
cc <- complete.cases(raw[, pred_vars])
dat <- raw[cc, ]
cat(sprintf("Complete cases: %d rows (%d excluded)\n", nrow(dat), sum(!cc)))

# Outcome dichotomizations (ROP: 0=none, 1=mild, 2=severe)
dat$severe_rop <- as.integer(dat$ROP == 2)   # 15 events
dat$any_rop    <- as.integer(dat$ROP >= 1)    # 48 events

cat(sprintf("Severe ROP: %d events among %d subjects (%.1f%%)\n",
            sum(tapply(dat$severe_rop, dat$Case_ID, max)),
            length(unique(dat$Case_ID)),
            100 * mean(tapply(dat$severe_rop, dat$Case_ID, max))))
cat(sprintf("Any ROP:    %d events among %d subjects (%.1f%%)\n",
            sum(tapply(dat$any_rop, dat$Case_ID, max)),
            length(unique(dat$Case_ID)),
            100 * mean(tapply(dat$any_rop, dat$Case_ID, max))))

# ── 2. DEPENDENCE STRUCTURE (Appendix Table) ────────────────────────────────

cat("\n========== DEPENDENCE STRUCTURE: ALL 7 OXYGENATION PREDICTORS ==========\n")

n0 <- mean(table(dat$Case_ID))  # average cluster size

# ICC via one-way random effects ANOVA for each predictor
compute_icc <- function(x, group) {
  aov_fit <- aov(x ~ factor(group))
  ms <- summary(aov_fit)[[1]]
  ms_b <- ms["Mean Sq"][1, 1]
  ms_w <- ms["Mean Sq"][2, 1]
  (ms_b - ms_w) / (ms_b + (n0 - 1) * ms_w)
}

# Lag-1 autocorrelation (median and IQR) for each predictor
compute_lag1 <- function(var_name) {
  lag1 <- sapply(unique(dat$Case_ID), function(id) {
    x <- dat[[var_name]][dat$Case_ID == id]
    if (length(x) < 3) return(NA)
    cor(x[-length(x)], x[-1])
  })
  c(median = median(lag1, na.rm = TRUE),
    q25    = unname(quantile(lag1, 0.25, na.rm = TRUE)),
    q75    = unname(quantile(lag1, 0.75, na.rm = TRUE)))
}

# Build summary table
dep_table <- data.frame(
  Predictor = oxy_vars,
  Mean = sapply(oxy_vars, function(v) mean(dat[[v]], na.rm = TRUE)),
  SD   = sapply(oxy_vars, function(v) sd(dat[[v]], na.rm = TRUE)),
  ICC  = sapply(oxy_vars, function(v) compute_icc(dat[[v]], dat$Case_ID)),
  stringsAsFactors = FALSE
)

lag1_results <- t(sapply(oxy_vars, compute_lag1))
dep_table$Lag1_Median <- lag1_results[, "median"]
dep_table$Lag1_Q25    <- lag1_results[, "q25"]
dep_table$Lag1_Q75    <- lag1_results[, "q75"]

# Print dependence structure table
cat("\nPredictor Summaries:\n")
cat(sprintf("%-22s  %10s  %10s  %6s  %s\n",
            "Predictor", "Mean", "SD", "ICC", "Lag-1 rho [IQR]"))
cat(paste(rep("-", 80), collapse = ""), "\n")
for (i in 1:nrow(dep_table)) {
  cat(sprintf("%-22s  %10.3f  %10.3f  %6.2f  %.2f [%.2f, %.2f]\n",
              dep_table$Predictor[i],
              dep_table$Mean[i], dep_table$SD[i], dep_table$ICC[i],
              dep_table$Lag1_Median[i], dep_table$Lag1_Q25[i], dep_table$Lag1_Q75[i]))
}

# Full 7x7 correlation matrix (daily observations)
cor_mat <- cor(dat[, oxy_vars], use = "complete.obs")

cat("\n7x7 Pairwise Correlation Matrix (daily observations):\n\n")
# Print header
cat(sprintf("%-22s", ""))
for (v in oxy_vars) cat(sprintf("  %8s", substr(v, 1, 8)))
cat("\n")
for (i in seq_along(oxy_vars)) {
  cat(sprintf("%-22s", oxy_vars[i]))
  for (j in seq_along(oxy_vars)) {
    if (j <= i) {
      cat(sprintf("  %8.2f", cor_mat[i, j]))
    } else {
      cat(sprintf("  %8s", ""))
    }
  }
  cat("\n")
}

# ── 3. MODEL FITTING FUNCTIONS ──────────────────────────────────────────────
# Mirrors simulation_factorial_v2.R exactly

fit_predict_ridge <- function(X_train, y_train, X_test, lambda = 0.01) {
  n_pos <- sum(y_train == 1)
  n_neg <- sum(y_train == 0)
  n_total <- length(y_train)
  w <- ifelse(y_train == 1, n_total / (2 * n_pos), n_total / (2 * n_neg))

  fit <- tryCatch({
    glmnet(X_train, y_train, family = "binomial", alpha = 0,
           weights = w, lambda = lambda, standardize = TRUE)
  }, error = function(e) {
    glmnet(X_train, y_train, family = "binomial", alpha = 0,
           weights = w, lambda = 0.1, standardize = TRUE)
  })

  as.numeric(predict(fit, newx = X_test, s = fit$lambda[1], type = "response"))
}

subject_auroc <- function(case_ids, pred_probs, y_outcome) {
  subj_preds <- tapply(pred_probs, case_ids, mean)
  subj_true  <- tapply(y_outcome, case_ids, function(x) x[1])
  subj_true  <- subj_true[names(subj_preds)]

  if (length(unique(subj_true)) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(as.numeric(subj_true),
                                   as.numeric(subj_preds), quiet = TRUE)))
}

# Return pROC::roc object for DeLong CI computation
subject_roc <- function(case_ids, pred_probs, y_outcome) {
  subj_preds <- tapply(pred_probs, case_ids, mean)
  subj_true  <- tapply(y_outcome, case_ids, function(x) x[1])
  subj_true  <- subj_true[names(subj_preds)]

  if (length(unique(subj_true)) < 2) return(NULL)
  pROC::roc(as.numeric(subj_true), as.numeric(subj_preds), quiet = TRUE)
}

# ── 4. CROSS-VALIDATION STRATEGIES ──────────────────────────────────────────

naive_kfold <- function(X, y_obs, case_ids, k = 10, seed = 0) {
  set.seed(seed)
  n <- nrow(X)
  fold_ids <- sample(rep(1:k, length.out = n))
  preds <- rep(NA_real_, n)

  for (fold in 1:k) {
    test_idx  <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)
    if (length(unique(y_obs[train_idx])) < 2) next
    preds[test_idx] <- fit_predict_ridge(
      X[train_idx, , drop = FALSE], y_obs[train_idx],
      X[test_idx, , drop = FALSE])
  }

  valid <- !is.na(preds)
  subject_auroc(case_ids[valid], preds[valid], y_obs)
}

cluster_kfold <- function(X, y_obs, case_ids, k = 10, seed = 0) {
  set.seed(seed)
  unique_ids <- unique(case_ids)
  n_subj <- length(unique_ids)
  actual_k <- min(k, n_subj)
  fold_map <- setNames(rep(1:actual_k, length.out = n_subj),
                        unique_ids[sample(n_subj)])
  obs_folds <- fold_map[as.character(case_ids)]
  preds <- rep(NA_real_, nrow(X))

  for (fold in 1:actual_k) {
    test_idx  <- which(obs_folds == fold)
    train_idx <- which(obs_folds != fold)
    if (length(test_idx) == 0 || length(train_idx) == 0) next
    if (length(unique(y_obs[train_idx])) < 2) next
    preds[test_idx] <- fit_predict_ridge(
      X[train_idx, , drop = FALSE], y_obs[train_idx],
      X[test_idx, , drop = FALSE])
  }

  valid <- !is.na(preds)
  subject_auroc(case_ids[valid], preds[valid], y_obs)
}

loco_cv <- function(X, y_obs, case_ids) {
  unique_ids <- unique(case_ids)
  subj_preds <- numeric(length(unique_ids))
  names(subj_preds) <- unique_ids

  for (i in seq_along(unique_ids)) {
    held_out <- unique_ids[i]
    test_idx  <- which(case_ids == held_out)
    train_idx <- which(case_ids != held_out)
    if (length(unique(y_obs[train_idx])) < 2) {
      subj_preds[i] <- 0.5; next
    }
    p <- fit_predict_ridge(
      X[train_idx, , drop = FALSE], y_obs[train_idx],
      X[test_idx, , drop = FALSE])
    subj_preds[i] <- mean(p)
  }

  subj_true <- tapply(y_obs, case_ids, function(x) x[1])
  subj_true <- subj_true[names(subj_preds)]
  if (length(unique(subj_true)) < 2) return(list(auc = NA_real_, roc = NULL))

  roc_obj <- pROC::roc(as.numeric(subj_true),
                        as.numeric(subj_preds), quiet = TRUE)
  list(auc = as.numeric(pROC::auc(roc_obj)), roc = roc_obj)
}

# ── 5. SUBJECT-LEVEL SUMMARY MODELS ─────────────────────────────────────────

summary_cv <- function(dat, outcome_var, use_slopes = FALSE, k = 10,
                       n_reps = 50, pred_vars) {
  unique_ids <- unique(dat$Case_ID)
  n_subj <- length(unique_ids)

  # Means
  subj_means <- do.call(rbind, lapply(unique_ids, function(id) {
    sub <- dat[dat$Case_ID == id, pred_vars, drop = FALSE]
    colMeans(sub, na.rm = TRUE)
  }))
  colnames(subj_means) <- paste0("mean_", pred_vars)
  features <- subj_means

  if (use_slopes) {
    subj_slopes <- do.call(rbind, lapply(unique_ids, function(id) {
      sub <- dat[dat$Case_ID == id, ]
      sapply(pred_vars[-1], function(v) {
        if (sd(sub$DOL) == 0) return(0)
        coef(lm(sub[[v]] ~ sub$DOL))[2]
      })
    }))
    colnames(subj_slopes) <- paste0("slope_", pred_vars[-1])
    features <- cbind(subj_means, subj_slopes)
  }

  subj_y <- sapply(unique_ids, function(id) dat[[outcome_var]][dat$Case_ID == id][1])
  X_subj <- as.matrix(features)

  # k-fold CV on subject-level data
  kfold_aucs <- sapply(1:n_reps, function(r) {
    set.seed(r)
    fold_ids <- sample(rep(1:k, length.out = n_subj))
    preds <- rep(NA_real_, n_subj)
    for (fold in 1:k) {
      test_idx  <- which(fold_ids == fold)
      train_idx <- which(fold_ids != fold)
      if (length(unique(subj_y[train_idx])) < 2) next
      preds[test_idx] <- fit_predict_ridge(
        X_subj[train_idx, , drop = FALSE], subj_y[train_idx],
        X_subj[test_idx, , drop = FALSE])
    }
    valid <- !is.na(preds)
    if (length(unique(subj_y[valid])) < 2) return(NA_real_)
    as.numeric(pROC::auc(pROC::roc(subj_y[valid], preds[valid], quiet = TRUE)))
  })

  # LOCO on subject-level data
  preds_loco <- rep(NA_real_, n_subj)
  for (i in 1:n_subj) {
    train_idx <- setdiff(1:n_subj, i)
    if (length(unique(subj_y[train_idx])) < 2) { preds_loco[i] <- 0.5; next }
    preds_loco[i] <- fit_predict_ridge(
      X_subj[train_idx, , drop = FALSE], subj_y[train_idx],
      X_subj[i, , drop = FALSE])
  }
  loco_roc <- pROC::roc(subj_y, preds_loco, quiet = TRUE)
  loco_auc <- as.numeric(pROC::auc(loco_roc))
  loco_ci  <- as.numeric(pROC::ci.auc(loco_roc, method = "delong"))

  list(kfold_mean = mean(kfold_aucs, na.rm = TRUE),
       kfold_sd   = sd(kfold_aucs, na.rm = TRUE),
       loco       = loco_auc,
       loco_ci    = loco_ci,
       n_features = ncol(X_subj))
}

# ── 6. RUN ANALYSIS ─────────────────────────────────────────────────────────

run_analysis <- function(dat, outcome_var, outcome_label, pred_vars, n_reps) {
  cat(sprintf("\n========== %s ==========\n", outcome_label))

  X <- as.matrix(dat[, pred_vars])
  y <- dat[[outcome_var]]
  case_ids <- dat$Case_ID

  # --- Observation-level CV strategies ---
  cat("Running naive 10-fold CV...\n")
  naive_aucs <- sapply(1:n_reps, function(r) {
    if (r %% 10 == 0) cat(sprintf("  rep %d/%d\n", r, n_reps))
    naive_kfold(X, y, case_ids, k = 10, seed = r)
  })

  cat("Running subject-level 10-fold CV...\n")
  subj10_aucs <- sapply(1:n_reps, function(r) {
    if (r %% 10 == 0) cat(sprintf("  rep %d/%d\n", r, n_reps))
    cluster_kfold(X, y, case_ids, k = 10, seed = r)
  })

  cat("Running LOCO CV...\n")
  loco_result <- loco_cv(X, y, case_ids)
  loco_auc <- loco_result$auc
  loco_roc <- loco_result$roc

  # --- CIs ---
  # DeLong 95% CI for LOCO
  loco_ci <- as.numeric(pROC::ci.auc(loco_roc, method = "delong"))
  # ci.auc returns: lower, auc, upper

  # Percentile intervals for k-fold AUROCs
  naive_pi  <- quantile(naive_aucs, c(0.025, 0.975))
  subj10_pi <- quantile(subj10_aucs, c(0.025, 0.975))

  # Percentile intervals for optimism
  naive_opt  <- naive_aucs - loco_auc
  subj10_dev <- subj10_aucs - loco_auc
  naive_opt_pi  <- quantile(naive_opt, c(0.025, 0.975))
  subj10_dev_pi <- quantile(subj10_dev, c(0.025, 0.975))

  cat("\n--- Observation-Level Results ---\n")
  cat(sprintf("Naive 10-fold:   mean = %.3f [%.3f, %.3f], SD = %.3f\n",
              mean(naive_aucs), naive_pi[1], naive_pi[2], sd(naive_aucs)))
  cat(sprintf("  Optimism: %+.3f [%+.3f, %+.3f]\n",
              mean(naive_opt), naive_opt_pi[1], naive_opt_pi[2]))
  cat(sprintf("Subject 10-fold: mean = %.3f [%.3f, %.3f], SD = %.3f\n",
              mean(subj10_aucs), subj10_pi[1], subj10_pi[2], sd(subj10_aucs)))
  cat(sprintf("  Deviation: %+.3f [%+.3f, %+.3f]\n",
              mean(subj10_dev), subj10_dev_pi[1], subj10_dev_pi[2]))
  cat(sprintf("LOCO:            %.3f (DeLong 95%% CI: %.3f–%.3f)\n",
              loco_auc, loco_ci[1], loco_ci[3]))

  # --- Subject-level summary models ---
  cat("\nRunning summary models (means only)...\n")
  summ_means <- summary_cv(dat, outcome_var, use_slopes = FALSE,
                           k = 10, n_reps = n_reps, pred_vars = pred_vars)

  cat("Running summary models (means + slopes)...\n")
  summ_slopes <- summary_cv(dat, outcome_var, use_slopes = TRUE,
                            k = 10, n_reps = n_reps, pred_vars = pred_vars)

  cat("\n--- Summary Model Results ---\n")
  cat(sprintf("Means only (%d features):\n", summ_means$n_features))
  cat(sprintf("  10-fold: mean = %.3f, SD = %.3f, optimism = %+.3f\n",
              summ_means$kfold_mean, summ_means$kfold_sd,
              summ_means$kfold_mean - loco_auc))
  cat(sprintf("  LOCO:    %.3f (DeLong 95%% CI: %.3f–%.3f), optimism = %+.3f\n",
              summ_means$loco, summ_means$loco_ci[1], summ_means$loco_ci[3],
              summ_means$loco - loco_auc))

  cat(sprintf("Means + slopes (%d features):\n", summ_slopes$n_features))
  cat(sprintf("  10-fold: mean = %.3f, SD = %.3f, optimism = %+.3f\n",
              summ_slopes$kfold_mean, summ_slopes$kfold_sd,
              summ_slopes$kfold_mean - loco_auc))
  cat(sprintf("  LOCO:    %.3f (DeLong 95%% CI: %.3f–%.3f), optimism = %+.3f\n",
              summ_slopes$loco, summ_slopes$loco_ci[1], summ_slopes$loco_ci[3],
              summ_slopes$loco - loco_auc))

  # Return all results
  list(naive_aucs = naive_aucs, subj10_aucs = subj10_aucs,
       loco = loco_auc, loco_ci = loco_ci,
       naive_pi = naive_pi, subj10_pi = subj10_pi,
       naive_opt_pi = naive_opt_pi, subj10_dev_pi = subj10_dev_pi,
       summ_means = summ_means, summ_slopes = summ_slopes)
}

# Run for both outcomes
res_severe <- run_analysis(dat, "severe_rop", "SEVERE ROP (15 events)",
                           pred_vars, N_CV_REPS)
res_any    <- run_analysis(dat, "any_rop", "ANY ROP (48 events)",
                           pred_vars, N_CV_REPS)

# ── 7. REGULARIZATION SENSITIVITY (Section 3.4) ─────────────────────────────

cat("\n========== REGULARIZATION SENSITIVITY ==========\n")

lambda_grid <- c(1.0, 0.01, 0.00001)  # strong, moderate, weak
lambda_labels <- c("strong (1.0)", "moderate (0.01)", "weak (0.00001)")

for (outcome_var in c("severe_rop", "any_rop")) {
  cat(sprintf("\n--- %s ---\n", outcome_var))
  X <- as.matrix(dat[, pred_vars])
  y <- dat[[outcome_var]]
  case_ids <- dat$Case_ID

  for (j in seq_along(lambda_grid)) {
    lam <- lambda_grid[j]

    naive_aucs_lam <- sapply(1:10, function(r) {
      set.seed(r)
      n <- nrow(X)
      fold_ids <- sample(rep(1:10, length.out = n))
      preds <- rep(NA_real_, n)
      for (fold in 1:10) {
        test_idx  <- which(fold_ids == fold)
        train_idx <- which(fold_ids != fold)
        if (length(unique(y[train_idx])) < 2) next
        preds[test_idx] <- fit_predict_ridge(
          X[train_idx, , drop = FALSE], y[train_idx],
          X[test_idx, , drop = FALSE], lambda = lam)
      }
      valid <- !is.na(preds)
      subject_auroc(case_ids[valid], preds[valid], y)
    })

    # LOCO with this lambda
    unique_ids <- unique(case_ids)
    subj_preds <- numeric(length(unique_ids))
    names(subj_preds) <- unique_ids
    for (i in seq_along(unique_ids)) {
      held_out <- unique_ids[i]
      test_idx  <- which(case_ids == held_out)
      train_idx <- which(case_ids != held_out)
      p <- fit_predict_ridge(
        X[train_idx, , drop = FALSE], y[train_idx],
        X[test_idx, , drop = FALSE], lambda = lam)
      subj_preds[i] <- mean(p)
    }
    subj_true <- tapply(y, case_ids, function(x) x[1])
    subj_true <- subj_true[names(subj_preds)]
    loco_roc_lam <- pROC::roc(as.numeric(subj_true),
                               as.numeric(subj_preds), quiet = TRUE)
    loco_lam <- as.numeric(pROC::auc(loco_roc_lam))
    loco_ci_lam <- as.numeric(pROC::ci.auc(loco_roc_lam, method = "delong"))

    cat(sprintf("  lambda = %-14s  naive = %.3f  LOCO = %.3f [%.3f, %.3f]  optimism = %+.3f\n",
                lambda_labels[j], mean(naive_aucs_lam), loco_lam,
                loco_ci_lam[1], loco_ci_lam[3],
                mean(naive_aucs_lam) - loco_lam))
  }
}

# ── 8. PRINT FORMATTED TABLES ───────────────────────────────────────────────

print_table <- function(res, label) {
  cat(sprintf("\n\n===== TABLE: %s =====\n", label))
  cat(sprintf("%-30s  %8s  %18s  %8s  %10s  %18s\n",
              "Strategy", "Mean AUC", "95% PI / CI", "SD", "Optimism", "95% PI"))
  cat(paste(rep("-", 110), collapse = ""), "\n")

  naive_opt <- mean(res$naive_aucs) - res$loco
  subj10_dev <- mean(res$subj10_aucs) - res$loco

  cat(sprintf("%-30s  %8.3f  [%5.3f, %5.3f]  %8.3f  %+10.3f  [%+.3f, %+.3f]\n",
              "Naive 10-fold",
              mean(res$naive_aucs), res$naive_pi[1], res$naive_pi[2],
              sd(res$naive_aucs), naive_opt,
              res$naive_opt_pi[1], res$naive_opt_pi[2]))
  cat(sprintf("%-30s  %8.3f  [%5.3f, %5.3f]  %8.3f  %+10.3f  [%+.3f, %+.3f]\n",
              "Subject 10-fold",
              mean(res$subj10_aucs), res$subj10_pi[1], res$subj10_pi[2],
              sd(res$subj10_aucs), subj10_dev,
              res$subj10_dev_pi[1], res$subj10_dev_pi[2]))
  cat(sprintf("%-30s  %8.3f  [%5.3f, %5.3f]  %8s  %10s\n",
              "LOCO (reference)", res$loco,
              res$loco_ci[1], res$loco_ci[3],
              "---", "ref"))

  cat(sprintf("%-30s  %8.3f  %18s  %8.3f  %+10.3f\n",
              "Summary 10f (means+slopes)", res$summ_slopes$kfold_mean,
              "", res$summ_slopes$kfold_sd,
              res$summ_slopes$kfold_mean - res$loco))
  cat(sprintf("%-30s  %8.3f  [%5.3f, %5.3f]  %8s  %+10.3f\n",
              "Summary LOCO (means+slopes)", res$summ_slopes$loco,
              res$summ_slopes$loco_ci[1], res$summ_slopes$loco_ci[3],
              "---", res$summ_slopes$loco - res$loco))

  cat(sprintf("%-30s  %8.3f  %18s  %8.3f  %+10.3f\n",
              "Summary 10f (means only)", res$summ_means$kfold_mean,
              "", res$summ_means$kfold_sd,
              res$summ_means$kfold_mean - res$loco))
  cat(sprintf("%-30s  %8.3f  [%5.3f, %5.3f]  %8s  %+10.3f\n",
              "Summary LOCO (means only)", res$summ_means$loco,
              res$summ_means$loco_ci[1], res$summ_means$loco_ci[3],
              "---", res$summ_means$loco - res$loco))
}

print_table(res_severe, "Severe ROP (N=101, 15 events)")
print_table(res_any,    "Any ROP (N=101, 48 events)")

cat("\n\nDone. Compare output with Tables 1-2 in manuscript.\n")

# Close output file
cat(sprintf("Results saved to:\n  %s\n", output_file))
