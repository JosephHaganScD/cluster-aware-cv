###############################################################################
# Full Factorial Simulation Study: Cluster-Aware Cross-Validation (v3)
# ====================================================================
# CHANGES FROM v2:
#   1. N_SIM increased from 100 to 500 datasets per condition
#   2. Monte Carlo 95% CIs computed for mean optimism and mean AUROC
#   3. Median [IQR] reported alongside mean for optimism
#   4. Proportion of failed/null datasets tracked per condition
#   5. Master RNG seed set for full reproducibility
#   6. Per-dataset results saved (not just condition summaries)
#   7. Output path updated
#
# Optimization: uses glmnet with fixed lambda (no inner CV)
# Expected runtime: ~20-40 hours (500 datasets x 162 conditions)
#
# Author: Joseph L. Hagan, ScD, MSPH
# Date:   April 2026
###############################################################################

library(glmnet)
library(pROC)

# ── USER SETTINGS ────────────────────────────────────────────────────────────

# Expected repository layout:
#   results/simulation/
OUT_DIR <- file.path("results", "simulation")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

N_SIM     <- 500   # datasets per condition
N_CV_REPS <- 10    # CV replicates per k-fold strategy
RESUME    <- TRUE  # skip completed conditions on restart

# Master seed for reproducibility
set.seed(20260410)

# ── 1. FACTORIAL DESIGN ─────────────────────────────────────────────────────

param_N          <- c(20, 30, 50, 75, 100, 150)
param_ICC        <- c(0.1, 0.3, 0.5)
param_AR1        <- c(0.2, 0.5, 0.8)
param_event_rate <- c(0.10, 0.20, 0.35)

FIXED_N_DAYS     <- 57
FIXED_N_PRED     <- 8
FIXED_BETA_SCALE <- 0.5

param_grid <- expand.grid(
  N          = param_N,
  ICC        = param_ICC,
  AR1        = param_AR1,
  event_rate = param_event_rate,
  stringsAsFactors = FALSE
)

cat(sprintf("Factorial design: %d conditions x %d sims = %d total datasets\n",
            nrow(param_grid), N_SIM, nrow(param_grid) * N_SIM))

# ── 2. DATA-GENERATING MECHANISM ────────────────────────────────────────────

generate_dataset <- function(n_subjects, n_days = 57, n_pred = 8,
                              icc = 0.4, ar1_rho = 0.55,
                              event_rate = 0.15, beta_scale = 0.5,
                              seed = 42) {
  set.seed(seed)

  sigma_b2 <- icc
  sigma_e2 <- 1.0 - icc
  sigma_b  <- sqrt(sigma_b2)
  sigma_e  <- sqrt(sigma_e2)
  sigma_innov <- sqrt(sigma_e2 * (1.0 - ar1_rho^2))

  n_fac <- min(3, n_pred)
  L <- matrix(rnorm(n_pred * n_fac), nrow = n_pred, ncol = n_fac) * 0.5
  Cov_b <- L %*% t(L) + diag(0.3, n_pred)
  sc <- sigma_b / sqrt(mean(diag(Cov_b)))
  Cov_b <- Cov_b * sc^2
  chol_b <- chol(Cov_b)

  true_beta <- rep(0, n_pred)
  true_beta[1] <- beta_scale
  true_beta[2] <- beta_scale * 0.7
  true_beta[3] <- -beta_scale * 0.5
  true_beta[4] <- beta_scale * 0.4

  b_mat <- matrix(rnorm(n_subjects * n_pred), nrow = n_subjects) %*% chol_b
  linear_preds <- as.numeric(b_mat %*% true_beta)

  lo <- -10; hi <- 10
  for (iter in 1:100) {
    mid <- (lo + hi) / 2
    if (mean(1 / (1 + exp(-(mid + linear_preds)))) < event_rate) lo <- mid else hi <- mid
  }
  intercept <- mid

  probs <- 1 / (1 + exp(-(intercept + linear_preds)))
  y <- rbinom(n_subjects, 1, probs)

  if (sum(y) < 2 || sum(y) > n_subjects - 2) {
    return(generate_dataset(n_subjects, n_days, n_pred, icc, ar1_rho,
                            event_rate, beta_scale, seed + 10000))
  }

  total_rows <- n_subjects * n_days
  X_mat <- matrix(NA_real_, nrow = total_rows, ncol = n_pred)
  case_ids <- integer(total_rows)
  dol_vec <- integer(total_rows)
  y_obs <- integer(total_rows)

  row_idx <- 0
  for (i in 1:n_subjects) {
    e_prev <- rnorm(n_pred) * sigma_e
    for (t in 1:n_days) {
      innov <- rnorm(n_pred) * sigma_innov
      e_curr <- ar1_rho * e_prev + innov
      row_idx <- row_idx + 1
      X_mat[row_idx, ] <- b_mat[i, ] + e_curr
      case_ids[row_idx] <- i
      dol_vec[row_idx] <- t
      y_obs[row_idx] <- y[i]
      e_prev <- e_curr
    }
  }

  X_full <- cbind(dol_vec, X_mat)
  colnames(X_full) <- c("DOL", paste0("X", 1:n_pred))

  list(X = X_full, y_obs = y_obs, case_ids = case_ids,
       y_subject = y, n_subjects = n_subjects, n_events = sum(y))
}

# ── 3. EMPIRICAL DATA CHARACTERISTICS ───────────────────────────────────────

compute_empirical_chars <- function(dat) {
  X <- dat$X
  ids <- dat$case_ids
  unique_ids <- unique(ids)
  n_subj <- length(unique_ids)
  n_obs <- nrow(X)

  iccs <- numeric(4)
  for (p in 1:4) {
    col_idx <- p + 1
    vals <- X[, col_idx]
    grand_mean <- mean(vals)
    group_means <- tapply(vals, ids, mean)
    group_sizes <- tapply(vals, ids, length)
    k <- length(group_means)
    ms_between <- sum(group_sizes * (group_means - grand_mean)^2) / (k - 1)
    within_resid <- vals - group_means[as.character(ids)]
    ms_within <- sum(within_resid^2) / (n_obs - k)
    n0 <- (n_obs - sum(group_sizes^2) / n_obs) / (k - 1)
    sigma_b2_hat <- (ms_between - ms_within) / n0
    if (sigma_b2_hat < 0) sigma_b2_hat <- 0
    iccs[p] <- sigma_b2_hat / (sigma_b2_hat + ms_within)
  }
  mean_icc <- mean(iccs)

  sample_ids <- unique_ids[1:min(50, n_subj)]
  ar1s <- numeric(0)
  for (s in sample_ids) {
    idx <- which(ids == s)
    idx <- idx[order(X[idx, 1])]
    for (p in 1:2) {
      col_idx <- p + 1
      vals <- X[idx, col_idx]
      if (length(vals) > 5) {
        ac <- cor(vals[-length(vals)], vals[-1])
        if (!is.na(ac)) ar1s <- c(ar1s, ac)
      }
    }
  }
  mean_ar1 <- if (length(ar1s) > 0) mean(ar1s) else 0

  list(empirical_icc = mean_icc, empirical_ar1 = mean_ar1)
}

# ── 4. CV STRATEGY FUNCTIONS ────────────────────────────────────────────────

fit_predict_ridge <- function(X_train, y_train, X_test) {
  n_pos <- sum(y_train == 1)
  n_neg <- sum(y_train == 0)
  n_total <- length(y_train)
  w <- ifelse(y_train == 1, n_total / (2 * n_pos), n_total / (2 * n_neg))

  fit <- tryCatch({
    glmnet(X_train, y_train, family = "binomial", alpha = 0,
           weights = w, lambda = 0.01, standardize = TRUE)
  }, error = function(e) {
    glmnet(X_train, y_train, family = "binomial", alpha = 0,
           weights = w, lambda = 0.1, standardize = TRUE)
  })

  as.numeric(predict(fit, newx = X_test, s = fit$lambda[1], type = "response"))
}

subject_auroc <- function(case_ids, pred_probs, y_subject, unique_ids) {
  subj_preds <- tapply(pred_probs, case_ids, mean)
  subj_true <- y_subject[as.integer(names(subj_preds))]
  if (length(unique(subj_true)) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(subj_true, subj_preds, quiet = TRUE)))
}

naive_kfold <- function(dat, k = 10, seed = 0) {
  set.seed(seed)
  n <- nrow(dat$X)
  fold_ids <- sample(rep(1:k, length.out = n))
  preds <- rep(NA_real_, n)
  for (fold in 1:k) {
    test_idx  <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)
    if (length(unique(dat$y_obs[train_idx])) < 2) next
    preds[test_idx] <- fit_predict_ridge(
      dat$X[train_idx, , drop = FALSE], dat$y_obs[train_idx],
      dat$X[test_idx, , drop = FALSE])
  }
  valid <- !is.na(preds)
  if (sum(valid) == 0) return(NA_real_)
  subject_auroc(dat$case_ids[valid], preds[valid], dat$y_subject, unique(dat$case_ids))
}

cluster_kfold <- function(dat, k = 10, seed = 0) {
  set.seed(seed)
  unique_ids <- unique(dat$case_ids)
  n_subj <- length(unique_ids)
  actual_k <- min(k, n_subj)
  fold_map <- setNames(rep(1:actual_k, length.out = n_subj),
                        unique_ids[sample(n_subj)])
  obs_folds <- fold_map[as.character(dat$case_ids)]
  preds <- rep(NA_real_, nrow(dat$X))
  for (fold in 1:actual_k) {
    test_idx  <- which(obs_folds == fold)
    train_idx <- which(obs_folds != fold)
    if (length(test_idx) == 0 || length(train_idx) == 0) next
    if (length(unique(dat$y_obs[train_idx])) < 2) next
    preds[test_idx] <- fit_predict_ridge(
      dat$X[train_idx, , drop = FALSE], dat$y_obs[train_idx],
      dat$X[test_idx, , drop = FALSE])
  }
  valid <- !is.na(preds)
  if (sum(valid) == 0) return(NA_real_)
  subject_auroc(dat$case_ids[valid], preds[valid], dat$y_subject, unique(dat$case_ids))
}

loco_cv <- function(dat) {
  unique_ids <- unique(dat$case_ids)
  subj_preds <- numeric(length(unique_ids))
  names(subj_preds) <- unique_ids
  for (i in seq_along(unique_ids)) {
    held_out <- unique_ids[i]
    test_idx  <- which(dat$case_ids == held_out)
    train_idx <- which(dat$case_ids != held_out)
    if (length(unique(dat$y_obs[train_idx])) < 2) {
      subj_preds[i] <- 0.5; next
    }
    p <- fit_predict_ridge(
      dat$X[train_idx, , drop = FALSE], dat$y_obs[train_idx],
      dat$X[test_idx, , drop = FALSE])
    subj_preds[i] <- mean(p)
  }
  subj_true <- dat$y_subject[as.integer(names(subj_preds))]
  if (length(unique(subj_true)) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(subj_true, subj_preds, quiet = TRUE)))
}

# ── 5. SINGLE-DATASET ANALYSIS ──────────────────────────────────────────────

analyze_one_dataset <- function(params, sim_idx, n_cv_reps = 10) {
  seed <- params$N * 100000 +
          as.integer(params$ICC * 10) * 10000 +
          as.integer(params$AR1 * 10) * 1000 +
          as.integer(params$event_rate * 100) * 10 +
          sim_idx

  dat <- generate_dataset(
    n_subjects = params$N, n_days = FIXED_N_DAYS, n_pred = FIXED_N_PRED,
    icc = params$ICC, ar1_rho = params$AR1, event_rate = params$event_rate,
    beta_scale = FIXED_BETA_SCALE, seed = seed)

  if (dat$n_events < 3 || dat$n_events > dat$n_subjects - 3) return(NULL)

  chars <- compute_empirical_chars(dat)

  naive_aucs  <- sapply(1:n_cv_reps, function(r) naive_kfold(dat, k = 10, seed = r))
  subj10_aucs <- sapply(1:n_cv_reps, function(r) cluster_kfold(dat, k = 10, seed = r))
  subj5_aucs  <- sapply(1:n_cv_reps, function(r) cluster_kfold(dat, k = 5, seed = r))
  loco_auc    <- loco_cv(dat)

  naive_aucs  <- naive_aucs[!is.na(naive_aucs)]
  subj10_aucs <- subj10_aucs[!is.na(subj10_aucs)]
  subj5_aucs  <- subj5_aucs[!is.na(subj5_aucs)]

  if (length(naive_aucs) == 0 || is.na(loco_auc)) return(NULL)

  data.frame(
    true_N = params$N, true_ICC = params$ICC, true_AR1 = params$AR1,
    true_event_rate = params$event_rate, sim_idx = sim_idx,
    n_events = dat$n_events,
    actual_event_rate = dat$n_events / dat$n_subjects,
    empirical_icc = chars$empirical_icc,
    empirical_ar1 = chars$empirical_ar1,
    naive_10f_mean = mean(naive_aucs),
    naive_10f_sd   = if (length(naive_aucs) > 1) sd(naive_aucs) else NA_real_,
    subj_10f_mean  = if (length(subj10_aucs) > 0) mean(subj10_aucs) else NA_real_,
    subj_10f_sd    = if (length(subj10_aucs) > 1) sd(subj10_aucs) else NA_real_,
    subj_5f_mean   = if (length(subj5_aucs) > 0) mean(subj5_aucs) else NA_real_,
    subj_5f_sd     = if (length(subj5_aucs) > 1) sd(subj5_aucs) else NA_real_,
    loco_auroc     = loco_auc,
    optimism_naive = mean(naive_aucs) - loco_auc,
    deviation_s10  = if (length(subj10_aucs) > 0) mean(subj10_aucs) - loco_auc else NA_real_,
    deviation_s5   = if (length(subj5_aucs) > 0) mean(subj5_aucs) - loco_auc else NA_real_,
    stringsAsFactors = FALSE)
}

# ── 6. MAIN SIMULATION LOOP ─────────────────────────────────────────────────

results_file  <- file.path(OUT_DIR, "sim_results_v3.csv")
progress_file <- file.path(OUT_DIR, "sim_progress_v3.rds")

if (RESUME && file.exists(results_file)) {
  all_results <- read.csv(results_file, stringsAsFactors = FALSE)
  cat(sprintf("Resuming: found %d existing results\n", nrow(all_results)))
  completed <- if (file.exists(progress_file)) readRDS(progress_file) else character(0)
} else {
  all_results <- data.frame()
  completed <- character(0)
}

# Track failures per condition
failure_log <- data.frame(
  cond_key = character(0), n_attempted = integer(0),
  n_valid = integer(0), n_failed = integer(0),
  pct_failed = numeric(0), stringsAsFactors = FALSE)

total_conditions <- nrow(param_grid)
t_global_start <- Sys.time()

cat(sprintf("\nStarting simulation at %s\n", format(t_global_start)))
cat(sprintf("Conditions: %d | Datasets/condition: %d | CV reps: %d\n\n",
            total_conditions, N_SIM, N_CV_REPS))

for (cond_idx in 1:total_conditions) {
  params <- param_grid[cond_idx, ]
  cond_key <- sprintf("N%d_ICC%.1f_AR%.1f_ER%.2f",
                       params$N, params$ICC, params$AR1, params$event_rate)

  if (cond_key %in% completed) {
    cat(sprintf("[%d/%d] %s — SKIPPED\n", cond_idx, total_conditions, cond_key))
    next
  }

  t_cond_start <- Sys.time()
  cond_results <- list()
  n_valid <- 0

  for (sim in 1:N_SIM) {
    result <- tryCatch(
      analyze_one_dataset(params, sim, n_cv_reps = N_CV_REPS),
      error = function(e) NULL)
    if (!is.null(result)) {
      cond_results[[length(cond_results) + 1]] <- result
      n_valid <- n_valid + 1
    }
  }

  # Log failures for this condition
  n_failed <- N_SIM - n_valid
  failure_log <- rbind(failure_log, data.frame(
    cond_key = cond_key, n_attempted = N_SIM,
    n_valid = n_valid, n_failed = n_failed,
    pct_failed = round(100 * n_failed / N_SIM, 1),
    stringsAsFactors = FALSE))

  if (length(cond_results) > 0) {
    cond_df <- do.call(rbind, cond_results)
    all_results <- rbind(all_results, cond_df)
    elapsed <- as.numeric(difftime(Sys.time(), t_cond_start, units = "secs"))
    cat(sprintf("[%d/%d] %s — %d/%d valid (%d failed, %.0fs) | Opt=%+.3f Dev(s10)=%+.3f LOCO=%.3f\n",
                cond_idx, total_conditions, cond_key, n_valid, N_SIM, n_failed, elapsed,
                mean(cond_df$optimism_naive, na.rm = TRUE),
                mean(cond_df$deviation_s10, na.rm = TRUE),
                mean(cond_df$loco_auroc, na.rm = TRUE)))
  } else {
    cat(sprintf("[%d/%d] %s — 0/%d valid (ALL FAILED)\n",
                cond_idx, total_conditions, cond_key, N_SIM))
  }

  completed <- c(completed, cond_key)

  if (cond_idx %% 5 == 0 || cond_idx == total_conditions) {
    write.csv(all_results, results_file, row.names = FALSE)
    saveRDS(completed, progress_file)
    write.csv(failure_log, file.path(OUT_DIR, "sim_failures_v3.csv"), row.names = FALSE)
    cat(sprintf("  >> Checkpoint: %d results saved\n", nrow(all_results)))
  }
}

write.csv(all_results, results_file, row.names = FALSE)
saveRDS(completed, progress_file)
write.csv(failure_log, file.path(OUT_DIR, "sim_failures_v3.csv"), row.names = FALSE)
total_time <- as.numeric(difftime(Sys.time(), t_global_start, units = "mins"))
cat(sprintf("\nSimulation complete: %d results in %.1f minutes\n",
            nrow(all_results), total_time))

# ── 7. HELPER: Monte Carlo 95% CI ───────────────────────────────────────────

mc_ci <- function(x) {
  # Returns formatted (lower, upper) string
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 2) return("(NA, NA)")
  m <- mean(x)
  se <- sd(x) / sqrt(n)
  sprintf("(%.4f, %.4f)", m - 1.96 * se, m + 1.96 * se)
}

mc_ci_num <- function(x) {
  # Returns numeric lower and upper
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 2) return(c(NA_real_, NA_real_))
  m <- mean(x)
  se <- sd(x) / sqrt(n)
  c(m - 1.96 * se, m + 1.96 * se)
}

# ── 8. SUMMARY TABLES WITH CIs ──────────────────────────────────────────────

cat("\n", paste(rep("=", 120), collapse = ""), "\n")
cat("  SUMMARY: Optimism by Cluster Count (with Monte Carlo 95% CIs)\n")
cat(paste(rep("=", 120), collapse = ""), "\n\n")

for (n_val in param_N) {
  sub <- all_results[all_results$true_N == n_val, ]
  n_data <- nrow(sub)
  opt_m  <- mean(sub$optimism_naive, na.rm = TRUE)
  opt_ci <- mc_ci(sub$optimism_naive)
  opt_med <- median(sub$optimism_naive, na.rm = TRUE)
  opt_iqr <- quantile(sub$optimism_naive, c(0.25, 0.75), na.rm = TRUE)
  dev10_m  <- mean(sub$deviation_s10, na.rm = TRUE)
  dev10_ci <- mc_ci(sub$deviation_s10)
  loco_m  <- mean(sub$loco_auroc, na.rm = TRUE)
  loco_ci <- mc_ci(sub$loco_auroc)
  naive_m <- mean(sub$naive_10f_mean, na.rm = TRUE)
  naive_ci <- mc_ci(sub$naive_10f_mean)

  cat(sprintf("N = %d (n = %d datasets)\n", n_val, n_data))
  cat(sprintf("  LOCO AUROC:       mean = %.4f  95%% CI %s\n", loco_m, loco_ci))
  cat(sprintf("  Naive AUROC:      mean = %.4f  95%% CI %s\n", naive_m, naive_ci))
  cat(sprintf("  Optimism(naive):  mean = %+.4f  95%% CI %s\n", opt_m, opt_ci))
  cat(sprintf("                    median = %+.4f  IQR [%+.4f, %+.4f]\n",
              opt_med, opt_iqr[1], opt_iqr[2]))
  cat(sprintf("  Dev(subj-10f):    mean = %+.4f  95%% CI %s\n", dev10_m, dev10_ci))
  cat("\n")
}

cat("\n", paste(rep("=", 120), collapse = ""), "\n")
cat("  SUMMARY: Optimism by ICC\n")
cat(paste(rep("=", 120), collapse = ""), "\n\n")

for (icc_val in param_ICC) {
  sub <- all_results[all_results$true_ICC == icc_val, ]
  opt_m <- mean(sub$optimism_naive, na.rm = TRUE)
  opt_ci <- mc_ci(sub$optimism_naive)
  opt_med <- median(sub$optimism_naive, na.rm = TRUE)
  opt_iqr <- quantile(sub$optimism_naive, c(0.25, 0.75), na.rm = TRUE)
  cat(sprintf("  ICC = %.1f | Opt(naive): mean = %+.4f  95%% CI %s  median = %+.4f [%+.4f, %+.4f]\n",
              icc_val, opt_m, opt_ci, opt_med, opt_iqr[1], opt_iqr[2]))
}

cat("\n", paste(rep("=", 120), collapse = ""), "\n")
cat("  SUMMARY: Optimism by AR(1)\n")
cat(paste(rep("=", 120), collapse = ""), "\n\n")

for (ar_val in param_AR1) {
  sub <- all_results[all_results$true_AR1 == ar_val, ]
  opt_m <- mean(sub$optimism_naive, na.rm = TRUE)
  opt_ci <- mc_ci(sub$optimism_naive)
  opt_med <- median(sub$optimism_naive, na.rm = TRUE)
  opt_iqr <- quantile(sub$optimism_naive, c(0.25, 0.75), na.rm = TRUE)
  cat(sprintf("  AR1 = %.1f | Opt(naive): mean = %+.4f  95%% CI %s  median = %+.4f [%+.4f, %+.4f]\n",
              ar_val, opt_m, opt_ci, opt_med, opt_iqr[1], opt_iqr[2]))
}

cat("\n", paste(rep("=", 120), collapse = ""), "\n")
cat("  SUMMARY: Optimism by Event Rate\n")
cat(paste(rep("=", 120), collapse = ""), "\n\n")

for (er_val in param_event_rate) {
  sub <- all_results[all_results$true_event_rate == er_val, ]
  opt_m <- mean(sub$optimism_naive, na.rm = TRUE)
  opt_ci <- mc_ci(sub$optimism_naive)
  opt_med <- median(sub$optimism_naive, na.rm = TRUE)
  opt_iqr <- quantile(sub$optimism_naive, c(0.25, 0.75), na.rm = TRUE)
  cat(sprintf("  ER = %.2f | Opt(naive): mean = %+.4f  95%% CI %s  median = %+.4f [%+.4f, %+.4f]\n",
              er_val, opt_m, opt_ci, opt_med, opt_iqr[1], opt_iqr[2]))
}

# N x ICC cross-tabulation
cat("\n", paste(rep("=", 120), collapse = ""), "\n")
cat("  Optimism by N x ICC (mean [95% CI])\n")
cat(paste(rep("=", 120), collapse = ""), "\n\n")

NxICC <- aggregate(optimism_naive ~ true_N + true_ICC,
                    data = all_results, FUN = mean, na.rm = TRUE)
print(reshape(NxICC, idvar = "true_N", timevar = "true_ICC", direction = "wide"))

cat("\n  With Monte Carlo 95% CIs:\n\n")
cat(sprintf("  %-5s", "N"))
for (icc_val in param_ICC) cat(sprintf(" | ICC=%.1f                         ", icc_val))
cat("\n")
for (n_val in param_N) {
  cat(sprintf("  %-5d", n_val))
  for (icc_val in param_ICC) {
    sub <- all_results[all_results$true_N == n_val & all_results$true_ICC == icc_val, ]
    m <- mean(sub$optimism_naive, na.rm = TRUE)
    ci <- mc_ci(sub$optimism_naive)
    cat(sprintf(" | %+.4f %s", m, ci))
  }
  cat("\n")
}

# N x AR1 cross-tabulation
cat("\n  Optimism by N x AR(1) (mean [95% CI])\n\n")
cat(sprintf("  %-5s", "N"))
for (ar_val in param_AR1) cat(sprintf(" | AR1=%.1f                         ", ar_val))
cat("\n")
for (n_val in param_N) {
  cat(sprintf("  %-5d", n_val))
  for (ar_val in param_AR1) {
    sub <- all_results[all_results$true_N == n_val & all_results$true_AR1 == ar_val, ]
    m <- mean(sub$optimism_naive, na.rm = TRUE)
    ci <- mc_ci(sub$optimism_naive)
    cat(sprintf(" | %+.4f %s", m, ci))
  }
  cat("\n")
}

# ICC x AR1 cross-tabulation
cat("\n  Optimism by ICC x AR(1) (mean [95% CI])\n\n")
cat(sprintf("  %-7s", "ICC"))
for (ar_val in param_AR1) cat(sprintf(" | AR1=%.1f                         ", ar_val))
cat("\n")
for (icc_val in param_ICC) {
  cat(sprintf("  %-7.1f", icc_val))
  for (ar_val in param_AR1) {
    sub <- all_results[all_results$true_ICC == icc_val & all_results$true_AR1 == ar_val, ]
    m <- mean(sub$optimism_naive, na.rm = TRUE)
    ci <- mc_ci(sub$optimism_naive)
    cat(sprintf(" | %+.4f %s", m, ci))
  }
  cat("\n")
}

# ── 9. FAILURE REPORT ───────────────────────────────────────────────────────

cat("\n", paste(rep("=", 120), collapse = ""), "\n")
cat("  FAILURE REPORT\n")
cat(paste(rep("=", 120), collapse = ""), "\n\n")

# Overall
total_attempted <- sum(failure_log$n_attempted)
total_failed <- sum(failure_log$n_failed)
cat(sprintf("  Overall: %d/%d datasets failed (%.1f%%)\n\n",
            total_failed, total_attempted, 100 * total_failed / total_attempted))

# By N
cat("  Failure rate by N:\n")
for (n_val in param_N) {
  fl_sub <- failure_log[grepl(sprintf("^N%d_", n_val), failure_log$cond_key), ]
  cat(sprintf("    N = %3d: %d/%d failed (%.1f%%)\n",
              n_val, sum(fl_sub$n_failed), sum(fl_sub$n_attempted),
              100 * sum(fl_sub$n_failed) / sum(fl_sub$n_attempted)))
}

# By event rate
cat("\n  Failure rate by event rate:\n")
for (er_val in param_event_rate) {
  er_str <- sprintf("ER%.2f", er_val)
  fl_sub <- failure_log[grepl(er_str, failure_log$cond_key), ]
  cat(sprintf("    ER = %.2f: %d/%d failed (%.1f%%)\n",
              er_val, sum(fl_sub$n_failed), sum(fl_sub$n_attempted),
              100 * sum(fl_sub$n_failed) / sum(fl_sub$n_attempted)))
}

# Conditions with >10% failure
high_fail <- failure_log[failure_log$pct_failed > 10, ]
if (nrow(high_fail) > 0) {
  cat(sprintf("\n  Conditions with >10%% failure rate: %d\n", nrow(high_fail)))
  for (i in 1:nrow(high_fail)) {
    cat(sprintf("    %s: %d/%d failed (%.1f%%)\n",
                high_fail$cond_key[i], high_fail$n_failed[i],
                high_fail$n_attempted[i], high_fail$pct_failed[i]))
  }
} else {
  cat("\n  No conditions with >10% failure rate.\n")
}

# ── 10. SAVE COMPREHENSIVE SUMMARY CSV ──────────────────────────────────────

# Per-condition summary with CIs, median, IQR, and failure info
cond_summary <- do.call(rbind, lapply(1:nrow(param_grid), function(i) {
  params <- param_grid[i, ]
  sub <- all_results[all_results$true_N == params$N &
                     all_results$true_ICC == params$ICC &
                     all_results$true_AR1 == params$AR1 &
                     all_results$true_event_rate == params$event_rate, ]

  cond_key <- sprintf("N%d_ICC%.1f_AR%.1f_ER%.2f",
                       params$N, params$ICC, params$AR1, params$event_rate)
  fl_row <- failure_log[failure_log$cond_key == cond_key, ]

  n_valid <- nrow(sub)
  if (n_valid == 0) return(NULL)

  opt_ci <- mc_ci_num(sub$optimism_naive)
  loco_ci <- mc_ci_num(sub$loco_auroc)
  naive_ci <- mc_ci_num(sub$naive_10f_mean)
  s10_ci <- mc_ci_num(sub$subj_10f_mean)
  dev10_ci <- mc_ci_num(sub$deviation_s10)

  data.frame(
    true_N = params$N, true_ICC = params$ICC,
    true_AR1 = params$AR1, true_event_rate = params$event_rate,
    n_valid = n_valid,
    n_failed = if (nrow(fl_row) > 0) fl_row$n_failed else NA_integer_,
    pct_failed = if (nrow(fl_row) > 0) fl_row$pct_failed else NA_real_,
    # LOCO AUROC
    loco_mean = mean(sub$loco_auroc, na.rm = TRUE),
    loco_ci_lo = loco_ci[1], loco_ci_hi = loco_ci[2],
    # Naive AUROC
    naive_mean = mean(sub$naive_10f_mean, na.rm = TRUE),
    naive_ci_lo = naive_ci[1], naive_ci_hi = naive_ci[2],
    # Subject-level 10-fold AUROC
    subj10_mean = mean(sub$subj_10f_mean, na.rm = TRUE),
    subj10_ci_lo = s10_ci[1], subj10_ci_hi = s10_ci[2],
    # Subject-level 5-fold AUROC
    subj5_mean = mean(sub$subj_5f_mean, na.rm = TRUE),
    # Optimism (naive)
    opt_mean = mean(sub$optimism_naive, na.rm = TRUE),
    opt_ci_lo = opt_ci[1], opt_ci_hi = opt_ci[2],
    opt_median = median(sub$optimism_naive, na.rm = TRUE),
    opt_q25 = as.numeric(quantile(sub$optimism_naive, 0.25, na.rm = TRUE)),
    opt_q75 = as.numeric(quantile(sub$optimism_naive, 0.75, na.rm = TRUE)),
    # Deviation (subject-level 10-fold)
    dev_s10_mean = mean(sub$deviation_s10, na.rm = TRUE),
    dev_s10_ci_lo = dev10_ci[1], dev_s10_ci_hi = dev10_ci[2],
    # Deviation (subject-level 5-fold)
    dev_s5_mean = mean(sub$deviation_s5, na.rm = TRUE),
    # Within-dataset SDs (averaged)
    naive_sd_mean = mean(sub$naive_10f_sd, na.rm = TRUE),
    subj10_sd_mean = mean(sub$subj_10f_sd, na.rm = TRUE),
    # Empirical characteristics
    emp_icc_mean = mean(sub$empirical_icc, na.rm = TRUE),
    emp_ar1_mean = mean(sub$empirical_ar1, na.rm = TRUE),
    stringsAsFactors = FALSE)
}))

summary_file <- file.path(OUT_DIR, "sim_summary_v3.csv")
write.csv(cond_summary, summary_file, row.names = FALSE)

cat(sprintf("\n\nOutput files:\n"))
cat(sprintf("  Per-dataset results: %s\n", results_file))
cat(sprintf("  Condition summary:   %s\n", summary_file))
cat(sprintf("  Failure log:         %s\n", file.path(OUT_DIR, "sim_failures_v3.csv")))
cat(sprintf("  Total: %d results in %.1f minutes\n", nrow(all_results), total_time))
cat("\nDone. Share sim_summary_v3.csv and sim_results_v3.csv with Claude.\n")



