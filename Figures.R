###############################################################################
# Manuscript Figures — cluster-aware CV (v5 simulation)
# =============================================================================
# Produces all four figures for submission to Diagnostic and Prognostic Research:
#
#   Figure 1  CV schematic (naive / subject-level k-fold / LOCO)
#             Source: hard-coded layout; no simulation data required.
#
#   Figure 2  Partitioning component of optimism (dot-and-line, 3x3 panel)
#             x-axis: N; lines: ICC; facets: AR1 (columns) x event rate (rows).
#             Subject-level 10-fold deviation from LOCO overlaid in grey.
#             95% Monte Carlo CI bands on naive optimism.
#
#   Figure 3  Heat map of partitioning component across all 162 conditions
#             N (rows) x ICC (columns), facets: AR1 (columns) x event rate (rows).
#
#   Figure 4  Dependence-axis dissociation (two series per N panel)
#             Partitioning component (rising) and residual discrepancy (falling)
#             plotted against ICC at AR1 = 0.5, prevalence = 0.20.
#
# Input:  sim_condsummary_main.csv  (162-row condition summary from v5 run)
# Output: PDF (vector, submission quality) + PNG (300 dpi, for review)
#
# Requires: ggplot2, dplyr, patchwork, scales
#
# Author: Joseph L. Hagan, ScD, MSPH
# Date:   August 2026
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})

# ── USER SETTINGS ─────────────────────────────────────────────────────────────

# Path to the v5 condition summary (162 rows)
COND_FILE <- "C:/simulations/cluster aware cv/results/simulation_v5/sim_condsummary_main.csv"

# Output directory (figures saved here)
OUT_DIR <- "C:/simulations/cluster aware cv/results/simulation_v5"

# ── LOAD DATA ─────────────────────────────────────────────────────────────────

cond <- read.csv(COND_FILE, stringsAsFactors = FALSE) |>
  filter(signal_mode == "baseline") |>
  mutate(
    # Readable factor labels
    AR1_label = factor(paste0("AR(1) = ", true_AR1),
                       levels = c("AR(1) = 0.2", "AR(1) = 0.5", "AR(1) = 0.8")),
    ER_label  = factor(paste0("Event rate = ",
                               sprintf("%.0f%%", 100 * true_event_rate)),
                       levels = c("Event rate = 10%",
                                  "Event rate = 20%",
                                  "Event rate = 35%")),
    ICC_f     = factor(paste0("ICC = ", true_ICC),
                       levels = c("ICC = 0.1", "ICC = 0.3", "ICC = 0.5")),
    N_f       = factor(true_N, levels = c(20, 30, 50, 75, 100, 150)),
    # 95% Monte Carlo CI for naive optimism using naive_sd and n = 500
    opt_lo    = optimism_naive_loco - 1.96 * naive_sd / sqrt(500),
    opt_hi    = optimism_naive_loco + 1.96 * naive_sd / sqrt(500),
    # |deviation| of subject-level 10-fold from LOCO
    dev_s10   = s10_mean - loco_auroc
  )

cat(sprintf("Loaded %d conditions from %s\n", nrow(cond), COND_FILE))

# ── SHARED THEME ──────────────────────────────────────────────────────────────

theme_ms <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "grey90", linewidth = 0.3),
    strip.background  = element_rect(fill = "grey95", colour = "grey70"),
    strip.text        = element_text(size = 10),
    legend.position   = "bottom",
    legend.title      = element_blank(),
    axis.title        = element_text(size = 10),
    axis.text         = element_text(size = 9),
    plot.margin       = margin(6, 6, 6, 6)
  )

icc_colors <- c("ICC = 0.1" = "#4292C6",
                 "ICC = 0.3" = "#2171B5",
                 "ICC = 0.5" = "#084594")

###############################################################################
# FIGURE 1: CV schematic
###############################################################################

col_train   <- "#4E79A7"
col_test    <- "#F28E2B"
col_heldout <- "#B07AA1"
col_leak    <- "firebrick"

n_subj <- 5
n_time <- 8
subj_levels <- paste0("S", n_subj:1)

make_tile_df <- function(split_vec) {
  data.frame(
    subject = factor(rep(paste0("S", 1:n_subj), each = n_time),
                     levels = subj_levels),
    time    = rep(1:n_time, times = n_subj),
    split   = factor(split_vec, levels = c("Train", "Test", "Held out"))
  )
}

df_a <- make_tile_df(c(
  "Train","Train","Test", "Train","Train","Test", "Train","Train",
  "Train","Test", "Train","Train","Test", "Train","Train","Train",
  "Test", "Train","Train","Train","Train","Test", "Train","Train",
  "Train","Train","Test", "Train","Train","Train","Test", "Train",
  "Train","Test", "Train","Train","Train","Train","Train","Test"
))

df_b <- make_tile_df(c(
  rep("Train", n_time),
  rep("Train", n_time),
  rep("Train", n_time),
  rep("Test",  n_time),
  rep("Test",  n_time)
))

df_c <- make_tile_df(c(
  rep("Train",    n_time),
  rep("Train",    n_time),
  rep("Held out", n_time),
  rep("Train",    n_time),
  rep("Train",    n_time)
))

fill_scale <- scale_fill_manual(
  values = c(Train = col_train, Test = col_test, "Held out" = col_heldout),
  name   = NULL, drop = FALSE, guide = "none"
)

x_scale <- scale_x_continuous(breaks = 1:n_time,
                               labels = paste0("t", 1:n_time),
                               expand = c(0.04, 0.04))

base_tile_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(size = 8, color = "gray50"),
    axis.text.y     = element_text(size = 9, face = "bold"),
    plot.title      = element_text(face = "bold", size = 10,
                                   margin = margin(b = 4)),
    legend.position = "none"
  )

y_lim <- c(0.2, 5.6)

p1a <- ggplot(df_a, aes(x = time, y = subject, fill = split)) +
  geom_tile(color = "white", linewidth = 0.8, width = 0.88, height = 0.88) +
  annotate("rect",
           xmin = 0.52, xmax = 8.48, ymin = 3.52, ymax = 4.48,
           color = col_leak, fill = NA, linewidth = 1.1, linetype = "dashed") +
  annotate("segment",
           x = 4.5, xend = 4.5, y = 0.75, yend = 3.48,
           color = col_leak, linewidth = 0.85,
           arrow = arrow(length = unit(0.10, "inches"),
                         type = "open", ends = "last")) +
  annotate("text", x = 4.5, y = 0.52,
           label = "Within-subject leakage  \u2192  Inflated apparent AUROC",
           color = col_leak, size = 3.2, hjust = 0.5, fontface = "italic") +
  fill_scale + x_scale +
  coord_cartesian(clip = "off", ylim = y_lim) +
  labs(x = NULL, y = NULL, title = "A.  Naive CV (observation-level)") +
  base_tile_theme

p1b <- ggplot(df_b, aes(x = time, y = subject, fill = split)) +
  geom_tile(color = "white", linewidth = 0.8, width = 0.88, height = 0.88) +
  geom_hline(yintercept = 2.5, color = "gray25",
             linewidth = 1.1, linetype = "solid") +
  annotate("text", x = 8.6, y = 4.0, label = "Train",
           color = col_train, size = 2.8, hjust = 0, fontface = "bold") +
  annotate("text", x = 8.6, y = 1.5, label = "Test",
           color = col_test,  size = 2.8, hjust = 0, fontface = "bold") +
  annotate("text", x = 4.5, y = 0.52,
           label = "Subjects assigned intact  \u2192  No within-subject contamination",
           color = "gray15", size = 2.9, hjust = 0.5, fontface = "italic") +
  fill_scale + x_scale +
  coord_cartesian(clip = "off", ylim = y_lim) +
  labs(x = NULL, y = NULL, title = "B.  Subject-level k-fold CV") +
  base_tile_theme

p1c <- ggplot(df_c, aes(x = time, y = subject, fill = split)) +
  geom_tile(color = "white", linewidth = 0.8, width = 0.88, height = 0.88) +
  annotate("segment",
           x = 0.52, xend = 8.48, y = 3.48, yend = 3.48,
           color = col_heldout, linewidth = 0.9, linetype = "dotted") +
  annotate("segment",
           x = 0.52, xend = 8.48, y = 2.52, yend = 2.52,
           color = col_heldout, linewidth = 0.9, linetype = "dotted") +
  annotate("text", x = 8.6, y = 3.0, label = "Held\nout",
           color = col_heldout, size = 2.7, hjust = 0, fontface = "bold") +
  annotate("text", x = 4.5, y = 0.52,
           label = "One subject fully held out  \u2192  No within-subject contamination",
           color = "gray15", size = 2.9, hjust = 0.5, fontface = "italic") +
  fill_scale + x_scale +
  coord_cartesian(clip = "off", ylim = y_lim) +
  labs(x = "Time point", y = NULL,
       title = "C.  LOCO (Leave-One-Cluster-Out)") +
  base_tile_theme +
  theme(axis.title.x = element_text(size = 9, color = "gray40",
                                    margin = margin(t = 4)))

p1_legend <- ggplot() +
  annotate("rect", xmin=1.55, xmax=1.95, ymin=0.25, ymax=0.75,
           fill=col_train, color=NA) +
  annotate("text", x=2.10, y=0.50, label="Train",
           hjust=0, vjust=0.5, size=3.2, color="gray10") +
  annotate("rect", xmin=4.05, xmax=4.45, ymin=0.25, ymax=0.75,
           fill=col_test, color=NA) +
  annotate("text", x=4.60, y=0.50, label="Test",
           hjust=0, vjust=0.5, size=3.2, color="gray10") +
  annotate("rect", xmin=6.35, xmax=6.75, ymin=0.25, ymax=0.75,
           fill=col_heldout, color=NA) +
  annotate("text", x=6.90, y=0.50, label="Held out",
           hjust=0, vjust=0.5, size=3.2, color="gray10") +
  xlim(0, 10) + ylim(0, 1) + theme_void()

fig1 <- (p1a / p1b / p1c / p1_legend) +
  plot_layout(heights = c(1, 1, 1, 0.10), guides = "keep") &
  theme(legend.position = "none",
        plot.margin     = margin(t = 6, r = 55, b = 2, l = 6))

ggsave(file.path(OUT_DIR, "Figure1_schematic.pdf"),
       fig1, width = 6.5, height = 8.5, units = "in", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Figure1_schematic.png"),
       fig1, width = 6.5, height = 8.5, units = "in", dpi = 300)
cat("Figure 1 saved.\n")

###############################################################################
# FIGURE 2: Partitioning component — dot-and-line, 3x3 facet
###############################################################################
# Naive optimism (solid colored lines by ICC, with 95% MC CI ribbons) and
# subject-level 10-fold deviation from LOCO (dashed grey lines, near zero)
# plotted against N, faceted by AR1 (columns) x event rate (rows).

fig2 <- ggplot(cond) +
  # Grey CI ribbon for s10 deviation (near-zero band, all ICC collapsed)
  geom_ribbon(aes(x = true_N, ymin = dev_s10 - 0.005, ymax = dev_s10 + 0.005,
                  group = ICC_f),
              fill = "grey70", alpha = 0.25) +
  # Dashed grey line: subject-level 10-fold deviation from LOCO
  geom_line(aes(x = true_N, y = dev_s10, group = ICC_f),
            color = "grey55", linewidth = 0.4, linetype = "dashed") +
  # MC CI ribbon: naive optimism by ICC
  geom_ribbon(aes(x = true_N, ymin = opt_lo, ymax = opt_hi,
                  fill = ICC_f, group = ICC_f),
              alpha = 0.15) +
  # Solid lines: naive optimism by ICC
  geom_line(aes(x = true_N, y = optimism_naive_loco,
                color = ICC_f, group = ICC_f),
            linewidth = 0.75) +
  geom_point(aes(x = true_N, y = optimism_naive_loco,
                 color = ICC_f),
             size = 1.8) +
  # Reference line at zero
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  scale_color_manual(values = icc_colors) +
  scale_fill_manual(values  = icc_colors, guide = "none") +
  scale_x_continuous(breaks = c(20, 50, 75, 100, 150)) +
  scale_y_continuous(limits = c(-0.02, 0.28),
                     breaks = seq(0, 0.25, by = 0.05),
                     labels = label_number(accuracy = 0.01)) +
  facet_grid(ER_label ~ AR1_label) +
  labs(x = "Number of clusters (N)",
       y = "Partitioning component (AUROC units)") +
  theme_ms +
  theme(legend.text = element_text(size = 10))

ggsave(file.path(OUT_DIR, "Figure2_dot_line.pdf"),
       fig2, width = 9, height = 8, units = "in", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Figure2_dot_line.png"),
       fig2, width = 9, height = 8, units = "in", dpi = 300)
cat("Figure 2 saved.\n")

###############################################################################
# FIGURE 3: Heat map — partitioning component across all 162 conditions
###############################################################################
# N (y-axis, descending) x ICC (x-axis), faceted AR1 (columns) x ER (rows).
# Cell text = mean partitioning component (3 decimal places).

fig3 <- ggplot(cond,
               aes(x = ICC_f, y = N_f, fill = optimism_naive_loco)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", optimism_naive_loco)),
            size = 2.7) +
  scale_fill_distiller(palette  = "YlOrRd",
                       direction = 1,
                       name      = "Mean partitioning\ncomponent",
                       limits    = c(0.06, 0.26),
                       breaks    = seq(0.07, 0.25, by = 0.05)) +
  scale_y_discrete(limits = rev(levels(cond$N_f))) +
  facet_grid(ER_label ~ AR1_label) +
  labs(x = "Intraclass correlation coefficient (ICC)",
       y = "Number of clusters (N)") +
  theme_ms +
  theme(
    panel.grid        = element_blank(),
    legend.position   = "right",
    legend.title      = element_text(size = 9),
    legend.text       = element_text(size = 9)
  )

ggsave(file.path(OUT_DIR, "Figure3_heatmap.pdf"),
       fig3, width = 9, height = 8.5, units = "in", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Figure3_heatmap.png"),
       fig3, width = 9, height = 8.5, units = "in", dpi = 300)
cat("Figure 3 saved.\n")

###############################################################################
# FIGURE 4: Dependence-axis dissociation
###############################################################################
# Two series per N panel: partitioning component (rising, blue solid) and
# residual discrepancy (falling, red dashed), plotted against ICC.
# Reference slice: AR1 = 0.5, prevalence = 0.20.

fig4_dat <- cond |>
  filter(abs(true_AR1 - 0.5) < 0.01,
         abs(true_event_rate - 0.20) < 0.01) |>
  select(true_N, true_ICC, optimism_naive_loco, resid_s10_true) |>
  rename(N   = true_N,
         ICC = true_ICC,
         partitioning = optimism_naive_loco,
         residual     = resid_s10_true)

fig4_long <- bind_rows(
  fig4_dat |> mutate(component = "Partitioning component",
                     value     = partitioning),
  fig4_dat |> mutate(component = "Residual discrepancy",
                     value     = residual)
) |>
  select(N, ICC, component, value) |>
  mutate(
    N_label   = factor(paste0("N = ", N),
                       levels = paste0("N = ", c(20, 30, 50, 75, 100, 150))),
    component = factor(component,
                       levels = c("Partitioning component",
                                  "Residual discrepancy"))
  )

fig4 <- ggplot(fig4_long,
               aes(x = ICC, y = value,
                   colour   = component,
                   linetype = component,
                   shape    = component)) +
  geom_hline(yintercept = 0, colour = "grey60",
             linewidth = 0.4, linetype = "dashed") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5, fill = "white") +
  scale_x_continuous(breaks = c(0.1, 0.3, 0.5),
                     labels = c("0.1", "0.3", "0.5")) +
  scale_y_continuous(labels = label_number(accuracy = 0.01),
                     breaks = seq(-0.05, 0.30, by = 0.05)) +
  scale_colour_manual(values = c("Partitioning component" = "#1b6ca8",
                                  "Residual discrepancy"   = "#c0392b")) +
  scale_linetype_manual(values = c("Partitioning component" = "solid",
                                    "Residual discrepancy"   = "dashed")) +
  scale_shape_manual(values  = c("Partitioning component" = 16,
                                   "Residual discrepancy"   = 17)) +
  facet_wrap(~ N_label, nrow = 2, ncol = 3) +
  labs(x = "Intraclass correlation coefficient (ICC)",
       y = "AUROC difference") +
  theme_ms

ggsave(file.path(OUT_DIR, "Figure4_dissociation.pdf"),
       fig4, width = 7, height = 5, units = "in", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Figure4_dissociation.png"),
       fig4, width = 7, height = 5, units = "in", dpi = 300)
cat("Figure 4 saved.\n")

###############################################################################
# SUMMARY
###############################################################################
cat(sprintf("\nAll four figures written to:\n  %s\n", OUT_DIR))
cat("Files produced:\n")
for (f in c("Figure1_schematic.pdf",   "Figure1_schematic.png",
            "Figure2_dot_line.pdf",    "Figure2_dot_line.png",
            "Figure3_heatmap.pdf",     "Figure3_heatmap.png",
            "Figure4_dissociation.pdf","Figure4_dissociation.png")) {
  cat(sprintf("  %s\n", f))
}

