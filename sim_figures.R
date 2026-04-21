###############################################################################
# Simulation Figures and Updated Summary Statistics
# =========================================================================
# Reads sim_results_v3.csv and produces:
#   - Figure 1: Dot-and-line plot (3x3 panel, AR1 x event rate)
#   - Supplementary Figure: Heat map (3x3 panel, AR1 x event rate)
#   - Updated marginal summary statistics for manuscript Sections 3.5–3.8
#
# Author: Joseph L. Hagan, ScD, MSPH
# Date:   April 2026
###############################################################################

library(ggplot2)

# ── INPUT / OUTPUT PATHS ─────────────────────────────────────────────────────

# Expected repository layout:
#   results/simulation/sim_results_v3.csv
#   results/figures/
SIM_RESULTS_FILE <- file.path("results", "simulation", "sim_results_v3.csv")
OUT_DIR          <- file.path("results", "figures")

if (!file.exists(SIM_RESULTS_FILE)) {
  stop(sprintf("Simulation results file not found: %s\nRun simulation_factorial_v3_github.R first.", SIM_RESULTS_FILE))
}
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ── 1. READ AND SUMMARIZE ────────────────────────────────────────────────────

dat <- read.csv(SIM_RESULTS_FILE)
cat(sprintf("Loaded %d per-dataset rows from %d unique conditions\n",
            nrow(dat), nrow(unique(dat[, c("true_N","true_ICC","true_AR1","true_event_rate")]))))

# Condition-level summaries
cond <- aggregate(
  cbind(optimism_naive, deviation_s10, deviation_s5) ~ true_N + true_ICC + true_AR1 + true_event_rate,
  data = dat, FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), n = length(x))
)

# aggregate() with multi-value FUN creates matrix columns; unpack them
summary_df <- data.frame(
  N          = cond$true_N,
  ICC        = cond$true_ICC,
  AR1        = cond$true_AR1,
  event_rate = cond$true_event_rate,
  mean_opt   = cond$optimism_naive[, "mean"],
  se_opt     = cond$optimism_naive[, "se"],
  n_reps     = cond$optimism_naive[, "n"],
  mean_dev10 = cond$deviation_s10[, "mean"],
  se_dev10   = cond$deviation_s10[, "se"],
  mean_dev5  = cond$deviation_s5[, "mean"],
  se_dev5    = cond$deviation_s5[, "se"]
)

# Monte Carlo 95% CIs
summary_df$opt_lo <- summary_df$mean_opt - 1.96 * summary_df$se_opt
summary_df$opt_hi <- summary_df$mean_opt + 1.96 * summary_df$se_opt
summary_df$dev10_lo <- summary_df$mean_dev10 - 1.96 * summary_df$se_dev10
summary_df$dev10_hi <- summary_df$mean_dev10 + 1.96 * summary_df$se_dev10

# Readable factor labels for panels
summary_df$AR1_label <- factor(
  paste0("AR(1) = ", summary_df$AR1),
  levels = c("AR(1) = 0.2", "AR(1) = 0.5", "AR(1) = 0.8")
)
summary_df$ER_label <- factor(
  paste0("Event rate = ", sprintf("%.0f%%", 100 * summary_df$event_rate)),
  levels = c("Event rate = 10%", "Event rate = 20%", "Event rate = 35%")
)
summary_df$ICC_f <- factor(summary_df$ICC, levels = c(0.1, 0.3, 0.5),
                           labels = c("ICC = 0.1", "ICC = 0.3", "ICC = 0.5"))

# ── 2. FIGURE 1: DOT-AND-LINE PLOT ──────────────────────────────────────────
# Mean naive CV optimism (with MC 95% CI) on y-axis, N on x-axis,
# separate lines for ICC levels, subject-level 10-fold deviation overlaid.
# Paneled: columns = AR(1), rows = event rate.

# Reshape for overlaying naive optimism and subject-level deviation
plot_naive <- data.frame(
  N       = summary_df$N,
  ICC_f   = summary_df$ICC_f,
  AR1_label = summary_df$AR1_label,
  ER_label  = summary_df$ER_label,
  y       = summary_df$mean_opt,
  ylo     = summary_df$opt_lo,
  yhi     = summary_df$opt_hi,
  Strategy = "Naive 10-fold"
)

plot_subj <- data.frame(
  N       = summary_df$N,
  ICC_f   = summary_df$ICC_f,
  AR1_label = summary_df$AR1_label,
  ER_label  = summary_df$ER_label,
  y       = summary_df$mean_dev10,
  ylo     = summary_df$dev10_lo,
  yhi     = summary_df$dev10_hi,
  Strategy = "Subject-level 10-fold"
)

plot_df <- rbind(plot_naive, plot_subj)
plot_df$Strategy <- factor(plot_df$Strategy,
                           levels = c("Naive 10-fold", "Subject-level 10-fold"))

# Color palette: 3 ICC levels for naive (blue gradient), grey for subject-level
icc_colors <- c("ICC = 0.1" = "#4292C6",    # light blue
                "ICC = 0.3" = "#2171B5",    # medium blue
                "ICC = 0.5" = "#084594")    # dark blue

fig1 <- ggplot() +
  # Subject-level 10-fold: single grey band showing all ICC levels collapsed
  geom_ribbon(data = plot_subj,
              aes(x = N, ymin = ylo, ymax = yhi, group = ICC_f),
              fill = "grey80", alpha = 0.3) +
  geom_line(data = plot_subj,
            aes(x = N, y = y, group = ICC_f),
            color = "grey55", linewidth = 0.4, linetype = "dashed") +
  # Naive optimism: lines with CI ribbons by ICC
  geom_ribbon(data = plot_naive,
              aes(x = N, ymin = ylo, ymax = yhi, fill = ICC_f),
              alpha = 0.15) +
  geom_line(data = plot_naive,
            aes(x = N, y = y, color = ICC_f, group = ICC_f),
            linewidth = 0.7) +
  geom_point(data = plot_naive,
             aes(x = N, y = y, color = ICC_f),
             size = 1.8) +
  # Reference line at zero
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
  # Scales
  scale_color_manual(values = icc_colors, name = NULL) +
  scale_fill_manual(values = icc_colors, guide = "none") +
  scale_x_continuous(breaks = c(20, 50, 75, 100, 150)) +
  scale_y_continuous(limits = c(-0.03, 0.23),
                     breaks = seq(0, 0.20, by = 0.05)) +
  # Faceting
  facet_grid(ER_label ~ AR1_label) +
  # Labels
  labs(x = "Number of clusters (N)",
       y = "Mean optimism (AUROC units)",
       caption = "Solid lines: naive 10-fold CV optimism (with 95% Monte Carlo CI bands) by ICC level.\nDashed grey lines: subject-level 10-fold CV deviation from LOCO by ICC level (lines overlap near zero across all conditions,\nindicating ICC independence of subject-level CV deviation).") +
  # Theme
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.margin     = margin(t = -5),
    strip.background  = element_rect(fill = "grey95"),
    strip.text        = element_text(size = 10),
    panel.grid.minor  = element_blank(),
    plot.caption      = element_text(size = 8.5, hjust = 0, margin = margin(t = 8)),
    axis.title        = element_text(size = 10.5)
  )

ggsave(file.path(OUT_DIR, "Figure1_dot_line_optimism.pdf"),
       fig1, width = 9, height = 8, units = "in")
cat("Saved Figure 1 (PDF) to:", file.path(OUT_DIR, "Figure1_dot_line_optimism.pdf"), "\n")

# High-resolution PNG for sharing/submission (avoids screen-capture tools)
ggsave(file.path(OUT_DIR, "Figure1_dot_line_optimism.png"),
       fig1, width = 9, height = 8, units = "in", dpi = 300)
cat("Saved Figure 1 (PNG) to:", file.path(OUT_DIR, "Figure1_dot_line_optimism.png"), "\n")

# ── 3. SUPPLEMENTARY FIGURE: HEAT MAP ───────────────────────────────────────
# 6 x 3 grid (N x ICC) with color = mean optimism, paneled by AR1 x event rate.
# Annotate cells with mean optimism value and rep count where < 450.

# Tile positions use factors so spacing is uniform
hm_df <- summary_df
hm_df$N_f <- factor(hm_df$N, levels = sort(unique(hm_df$N)))
hm_df$ICC_num <- factor(hm_df$ICC)

# Text label: optimism value; add rep count only when substantially below 500
hm_df$cell_label <- sprintf("%.3f", hm_df$mean_opt)
hm_df$rep_label  <- ifelse(hm_df$n_reps < 450,
                           sprintf("%.3f\n(n=%d)", hm_df$mean_opt, hm_df$n_reps),
                           sprintf("%.3f", hm_df$mean_opt))

fig_hm <- ggplot(hm_df, aes(x = ICC_num, y = N_f, fill = mean_opt)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = rep_label), size = 2.7, lineheight = 0.85) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1,
                       name = "Mean\noptimism",
                       limits = c(0.03, 0.21),
                       breaks = seq(0.05, 0.20, by = 0.05)) +
  scale_y_discrete(limits = rev(levels(hm_df$N_f))) +
  facet_grid(ER_label ~ AR1_label) +
  labs(x = "ICC", y = "Number of clusters (N)",
       caption = "Cell values show mean naive CV optimism (AUROC units). Cells with (n=...) had fewer than 450 valid simulation replicates (of 500 attempted).") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position   = "right",
    strip.background  = element_rect(fill = "grey95", color = NA),
    strip.text        = element_text(size = 10),
    panel.grid        = element_blank(),
    plot.caption      = element_text(size = 8.5, hjust = 0, margin = margin(t = 8)),
    axis.title        = element_text(size = 10.5)
  )

ggsave(file.path(OUT_DIR, "SuppFigure_heatmap_optimism.pdf"),
       fig_hm, width = 9, height = 8.5, units = "in")
cat("Saved Supplementary Figure (PDF) to:", file.path(OUT_DIR, "SuppFigure_heatmap_optimism.pdf"), "\n")

# High-resolution PNG for sharing/submission
ggsave(file.path(OUT_DIR, "SuppFigure_heatmap_optimism.png"),
       fig_hm, width = 9, height = 8.5, units = "in", dpi = 300)
cat("Saved Supplementary Figure (PNG) to:", file.path(OUT_DIR, "SuppFigure_heatmap_optimism.png"), "\n")

# ── 4. UPDATED MANUSCRIPT SUMMARY STATISTICS ─────────────────────────────────
# Print all values needed to update Sections 3.5–3.8 in the manuscript.

cat("\n\n========== UPDATED MANUSCRIPT STATISTICS ==========\n")
cat(sprintf("Total per-dataset rows: %d\n", nrow(dat)))
cat(sprintf("Unique conditions: %d\n", nrow(summary_df)))
cat(sprintf("Grand mean optimism (across %d conditions): %.3f\n",
            nrow(summary_df), mean(summary_df$mean_opt)))
cat(sprintf("Optimism range across conditions: +%.3f to +%.3f\n",
            min(summary_df$mean_opt), max(summary_df$mean_opt)))

# Section 3.5: By cluster count
cat("\n--- Section 3.5: Optimism by cluster count ---\n")
by_N <- aggregate(optimism_naive ~ true_N, data = dat,
                  FUN = function(x) c(mean = mean(x), n = length(x)))
by_N_df <- data.frame(N = by_N$true_N,
                      mean_opt = by_N$optimism_naive[, "mean"],
                      n_datasets = by_N$optimism_naive[, "n"])
for (i in 1:nrow(by_N_df)) {
  cat(sprintf("  N = %3d: mean optimism = +%.3f  (%d datasets)\n",
              by_N_df$N[i], by_N_df$mean_opt[i], by_N_df$n_datasets[i]))
}

# Section 3.6: By ICC, AR1, event rate (marginal means)
cat("\n--- Section 3.6: Marginal effects ---\n")
cat("  By ICC:\n")
by_ICC <- aggregate(optimism_naive ~ true_ICC, data = dat, FUN = mean)
for (i in 1:nrow(by_ICC)) {
  cat(sprintf("    ICC = %.1f: mean optimism = +%.3f\n",
              by_ICC$true_ICC[i], by_ICC$optimism_naive[i]))
}
cat("  By AR(1):\n")
by_AR1 <- aggregate(optimism_naive ~ true_AR1, data = dat, FUN = mean)
for (i in 1:nrow(by_AR1)) {
  cat(sprintf("    AR(1) = %.1f: mean optimism = +%.3f\n",
              by_AR1$true_AR1[i], by_AR1$optimism_naive[i]))
}
cat("  By event rate:\n")
by_ER <- aggregate(optimism_naive ~ true_event_rate, data = dat, FUN = mean)
for (i in 1:nrow(by_ER)) {
  cat(sprintf("    ER = %.2f: mean optimism = +%.3f\n",
              by_ER$true_event_rate[i], by_ER$optimism_naive[i]))
}

# Section 3.7: Interaction tables (N x ICC, ICC x AR1)
cat("\n--- Section 3.7: N x ICC interaction (Table 4 values) ---\n")
tbl4 <- aggregate(optimism_naive ~ true_N + true_ICC, data = dat, FUN = mean)
tbl4_wide <- reshape(tbl4, idvar = "true_N", timevar = "true_ICC",
                     direction = "wide")
colnames(tbl4_wide) <- gsub("optimism_naive\\.", "ICC_", colnames(tbl4_wide))
tbl4_wide <- tbl4_wide[order(tbl4_wide$true_N), ]
cat(sprintf("  %5s  %10s  %10s  %10s\n", "N", "ICC=0.1", "ICC=0.3", "ICC=0.5"))
for (i in 1:nrow(tbl4_wide)) {
  cat(sprintf("  %5d  %+10.3f  %+10.3f  %+10.3f\n",
              tbl4_wide$true_N[i],
              tbl4_wide[i, "ICC_0.1"],
              tbl4_wide[i, "ICC_0.3"],
              tbl4_wide[i, "ICC_0.5"]))
}

cat("\n--- Section 3.7: ICC x AR(1) interaction (Table 5 values) ---\n")
tbl5 <- aggregate(optimism_naive ~ true_ICC + true_AR1, data = dat, FUN = mean)
tbl5_wide <- reshape(tbl5, idvar = "true_ICC", timevar = "true_AR1",
                     direction = "wide")
colnames(tbl5_wide) <- gsub("optimism_naive\\.", "AR1_", colnames(tbl5_wide))
tbl5_wide <- tbl5_wide[order(tbl5_wide$true_ICC), ]
cat(sprintf("  %7s  %10s  %10s  %10s\n", "ICC", "AR1=0.2", "AR1=0.5", "AR1=0.8"))
for (i in 1:nrow(tbl5_wide)) {
  cat(sprintf("  %7.1f  %+10.3f  %+10.3f  %+10.3f\n",
              tbl5_wide$true_ICC[i],
              tbl5_wide[i, "AR1_0.2"],
              tbl5_wide[i, "AR1_0.5"],
              tbl5_wide[i, "AR1_0.8"]))
}

# Section 3.8: Subject-level CV performance
cat("\n--- Section 3.8: Subject-level CV performance ---\n")
cat(sprintf("  Subject 10-fold: mean |deviation| = %.3f, max |deviation| = %.3f\n",
            mean(abs(summary_df$mean_dev10)),
            max(abs(summary_df$mean_dev10))))
cat(sprintf("  Subject  5-fold: mean |deviation| = %.3f, max |deviation| = %.3f\n",
            mean(abs(summary_df$mean_dev5)),
            max(abs(summary_df$mean_dev5))))

# Precision illusion: ratio of naive SD to subject-level SD
sd_compare <- aggregate(
  cbind(naive_10f_sd, subj_10f_sd) ~ true_N + true_ICC + true_AR1 + true_event_rate,
  data = dat, FUN = mean
)
sd_ratio <- sd_compare$naive_10f_sd / sd_compare$subj_10f_sd
cat(sprintf("  Precision illusion (naive SD / subject SD): median = %.2f, range = %.2f–%.2f\n",
            median(sd_ratio), min(sd_ratio), max(sd_ratio)))

# Abstract values
cat("\n--- Abstract update ---\n")
cat(sprintf("  Optimism range: +%.3f to +%.3f across %d conditions\n",
            min(summary_df$mean_opt), max(summary_df$mean_opt), nrow(summary_df)))
cat(sprintf("  Grand mean optimism: +%.3f\n", mean(summary_df$mean_opt)))

# Worst-case and best-case conditions
worst <- summary_df[which.max(summary_df$mean_opt), ]
best  <- summary_df[which.min(summary_df$mean_opt), ]
cat(sprintf("  Worst case: N=%d, ICC=%.1f, AR1=%.1f, ER=%.2f -> +%.3f\n",
            worst$N, worst$ICC, worst$AR1, worst$event_rate, worst$mean_opt))
cat(sprintf("  Best case:  N=%d, ICC=%.1f, AR1=%.1f, ER=%.2f -> +%.3f\n",
            best$N, best$ICC, best$AR1, best$event_rate, best$mean_opt))

# Datasets per N level (for Table 3 header note)
cat("\n--- Datasets per cluster count (for Table 3 note) ---\n")
for (i in 1:nrow(by_N_df)) {
  cat(sprintf("  N = %3d: %d datasets\n", by_N_df$N[i], by_N_df$n_datasets[i]))
}

# Conditions with < 500 reps (for footnote)
low_rep <- summary_df[summary_df$n_reps < 500, c("N","ICC","AR1","event_rate","n_reps")]
if (nrow(low_rep) > 0) {
  cat(sprintf("\n--- Conditions with < 500 valid reps (%d conditions) ---\n", nrow(low_rep)))
  low_rep <- low_rep[order(low_rep$n_reps), ]
  for (i in 1:min(10, nrow(low_rep))) {
    cat(sprintf("  N=%d, ICC=%.1f, AR1=%.1f, ER=%.2f: n=%d\n",
                low_rep$N[i], low_rep$ICC[i], low_rep$AR1[i],
                low_rep$event_rate[i], low_rep$n_reps[i]))
  }
  if (nrow(low_rep) > 10) cat(sprintf("  ... and %d more\n", nrow(low_rep) - 10))
}

cat("\n\nDone.\n")
