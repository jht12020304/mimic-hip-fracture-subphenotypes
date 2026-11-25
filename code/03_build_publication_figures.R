#!/usr/bin/env Rscript
# ==============================================================================
# Paper 1 v5 (Nature Communications) — Visual-First Publication Figures
#   - Figure 1: Minimalist visual schematic modeled on ICML/Nature style
#               (rounded boxes, smooth Bezier manifold, orthogonal arrow tree,
#                custom vector icons, zero prose paragraphs)
#   - Figures 2-5 & S1-S3: Pure visual-first plots (no embedded titles,
#               subtitles, floating statistical boxes, or redundant text columns;
#               all statistics reside in captions, tables, and main text)
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# Support running from repo root or v5_nature_comm/src
args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
SCRIPT_DIR <- if (length(script_arg) > 0) dirname(normalizePath(script_arg)) else getwd()

if (dir.exists(file.path(SCRIPT_DIR, "..", "tables"))) {
  BASE_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."))
} else {
  BASE_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."))
}

V5_T    <- file.path(BASE_DIR, "tables")
OUT_PNG <- file.path(BASE_DIR, "figures", "png_600dpi")
OUT_PDF <- file.path(BASE_DIR, "figures", "pdf_vector")
dir.create(OUT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_PDF, recursive = TRUE, showWarnings = FALSE)

# Remove temporary python file if present
tmp_py <- file.path(SCRIPT_DIR, "build_fig1_icml_style.py")
if (file.exists(tmp_py)) invisible(file.remove(tmp_py))

PAL_P <- c("P1" = "#2B5C8F", "P2" = "#00897B", "P3" = "#C0392B")
LBL_P <- c("P1" = "P1: Homeostatic", "P2" = "P2: Stress-Inflammatory", "P3" = "P3: Cardiorenal-Anemic")
PAL_P_FULL <- setNames(PAL_P, LBL_P)

theme_nc <- function(base_size = 9.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(color = "#111111"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.tag = element_text(face = "bold", size = base_size + 3.5, color = "#000000"),
      axis.title = element_text(face = "bold", size = base_size - 0.3, color = "#111111"),
      axis.text = element_text(size = base_size - 0.8, color = "#222222"),
      axis.line = element_line(color = "#333333", linewidth = 0.45),
      axis.ticks = element_line(color = "#333333", linewidth = 0.4),
      panel.grid.major = element_line(color = "#F1F5F9", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold", size = base_size - 0.8),
      legend.text = element_text(size = base_size - 0.9),
      legend.key.size = unit(0.36, "cm"),
      legend.background = element_rect(fill = "transparent", color = NA),
      strip.text = element_text(face = "bold", size = base_size - 0.5, color = "#111111"),
      strip.background = element_rect(fill = "#F8FAFC", color = "#CBD5E1", linewidth = 0.4),
      plot.margin = margin(6, 10, 6, 6)
    )
}

save_fig <- function(plt, stem, w, h, also_svg = FALSE) {
  ggsave(file.path(OUT_PNG, paste0(stem, ".png")), plt, width = w, height = h, dpi = 600, bg = "white")
  ggsave(file.path(OUT_PDF, paste0(stem, ".pdf")), plt, width = w, height = h, device = cairo_pdf, bg = "white")
  if (also_svg) {
    svg(file.path(OUT_PDF, paste0(stem, ".svg")), width = w, height = h, bg = "white")
    print(plt)
    dev.off()
  }
  cat(sprintf("Saved %s (%.1f x %.1f in)\n", stem, w, h))
}

# ==============================================================================
# FIGURE 1: Minimalist ICML / Nature Visual Schematic (Grid Graphics)
# ==============================================================================
draw_fig1_grob <- function() {
  # Coordinate helper on [0, 13.6] x [0, 4.6]
  ux <- function(x) unit(x / 13.6, "npc")
  uy <- function(y) unit(y / 4.6, "npc")
  uw <- function(w) unit(w / 13.6, "npc")
  uh <- function(h) unit(h / 4.6, "npc")

  rbox <- function(x, y, w, h, fill = "#FFFFFF", col = "#2D3748", lwd = 1.2, lty = "solid", r = 0.16) {
    roundrectGrob(
      x = ux(x + w / 2), y = uy(y + h / 2), width = uw(w), height = uh(h),
      r = unit(r, "snpc"),
      gp = gpar(fill = fill, col = col, lwd = lwd, lty = lty)
    )
  }

  arr <- function(x0, y0, x1, y1, col = "#1A202C", lwd = 1.5, lty = "solid", len_mm = 2.1) {
    segmentsGrob(
      x0 = ux(x0), y0 = uy(y0), x1 = ux(x1), y1 = uy(y1),
      gp = gpar(col = col, fill = col, lwd = lwd, lty = lty, lineend = "round"),
      arrow = arrow(type = "closed", length = unit(len_mm, "mm"), angle = 22)
    )
  }

  txt <- function(x, y, label, size = 9.0, face = "plain", col = "#111111", hjust = 0.5, vjust = 0.5, lh = 1.08) {
    textGrob(
      label, x = ux(x), y = uy(y), hjust = hjust, vjust = vjust,
      gp = gpar(fontsize = size, fontface = face, col = col, fontfamily = "sans", lineheight = lh)
    )
  }

  # Smooth wavy 2D manifold polygon (like icml2026-0058)
  manifold_grob <- function(cx, cy, w, h, fill = "#EEF3FA", col = "#5A6B82") {
    t_top <- seq(0, 1, length.out = 40)
    x_top <- (cx - w * 0.48) + w * 0.96 * t_top
    y_top <- (cy + h * 0.36) + h * 0.10 * sin(2 * pi * t_top)
    t_bot <- seq(1, 0, length.out = 40)
    x_bot <- (cx - w * 0.44) + w * 0.96 * t_bot
    y_bot <- (cy - h * 0.36) + h * 0.10 * sin(2 * pi * t_bot)
    polygonGrob(
      x = ux(c(x_top, x_bot)), y = uy(c(y_top, y_bot)),
      gp = gpar(fill = fill, col = col, lwd = 1.35)
    )
  }

  ellipse_grob <- function(cx, cy, rx, ry, angle_deg = 15, fill = "#D4E6F1", col = "#2B5C8F", lwd = 1.4) {
    th <- seq(0, 2 * pi, length.out = 60)
    ang <- angle_deg * pi / 180
    ex <- rx * cos(th) * cos(ang) - ry * sin(th) * sin(ang) + cx
    ey <- rx * cos(th) * sin(ang) + ry * sin(th) * cos(ang) + cy
    polygonGrob(x = ux(ex), y = uy(ey), gp = gpar(fill = fill, col = col, lwd = lwd))
  }

  # Icon 1: Blood tube for ED Pre-Admission
  icon_tube <- function(cx, cy) {
    gList(
      rbox(cx - 0.11, cy - 0.28, 0.22, 0.52, fill = "#FFFFFF", col = "#1A202C", lwd = 1.2, r = 0.35),
      rbox(cx - 0.085, cy - 0.25, 0.17, 0.28, fill = "#E6B0AA", col = NA, lwd = 0, r = 0.30),
      rbox(cx - 0.14, cy + 0.22, 0.28, 0.11, fill = "#C0392B", col = "#1A202C", lwd = 1.1, r = 0.22)
    )
  }

  # Icon 2: Survival / calibration curve icon (Case 1 style)
  icon_surv <- function(cx, cy) {
    x0 <- cx - 0.38; y0 <- cy - 0.32
    xs <- seq(0, 1, length.out = 25)
    gList(
      linesGrob(x = ux(c(x0, x0, x0 + 0.78)), y = uy(c(y0 + 0.64, y0, y0)), gp = gpar(col = "#111111", lwd = 1.6)),
      linesGrob(x = ux(x0 + 0.04 + 0.72 * xs), y = uy(y0 + 0.58 - 0.12 * xs^1.1), gp = gpar(col = "#2B5C8F", lwd = 1.9)),
      linesGrob(x = ux(x0 + 0.04 + 0.72 * xs), y = uy(y0 + 0.54 - 0.22 * xs^1.2), gp = gpar(col = "#00897B", lwd = 1.9)),
      linesGrob(x = ux(x0 + 0.04 + 0.72 * xs), y = uy(y0 + 0.50 - 0.40 * xs^0.85), gp = gpar(col = "#C0392B", lwd = 1.9)),
      segmentsGrob(x0 = ux(x0 + 0.38), y0 = uy(y0), x1 = ux(x0 + 0.38), y1 = uy(y0 + 0.60),
                   gp = gpar(col = "#64748B", lwd = 1.0, lty = "dashed"))
    )
  }

  # Icon 3: Stylized paired kidney + landmark icon (Case 2 style)
  icon_kidney <- function(cx, cy) {
    gList(
      ellipse_grob(cx - 0.17, cy + 0.04, 0.15, 0.26, angle_deg = 16, fill = "#AED6F1", col = "#111111", lwd = 1.4),
      ellipse_grob(cx + 0.17, cy + 0.04, 0.15, 0.26, angle_deg = -16, fill = "#F5B7B1", col = "#111111", lwd = 1.4),
      ellipse_grob(cx, cy + 0.03, 0.09, 0.10, angle_deg = 0, fill = "#FFFFFF", col = NA, lwd = 0),
      segmentsGrob(x0 = ux(cx - 0.05), y0 = uy(cy + 0.02), x1 = ux(cx - 0.05), y1 = uy(cy - 0.32), gp = gpar(col = "#111111", lwd = 1.4)),
      segmentsGrob(x0 = ux(cx + 0.05), y0 = uy(cy + 0.02), x1 = ux(cx + 0.05), y1 = uy(cy - 0.32), gp = gpar(col = "#111111", lwd = 1.4))
    )
  }

  # Icon 4: Two-arm clone branching icon (Case 3 style)
  icon_ccw <- function(cx, cy) {
    gList(
      ellipse_grob(cx - 0.30, cy, 0.10, 0.10, 0, fill = "#E2E8F0", col = "#111111", lwd = 1.5),
      arr(cx - 0.19, cy + 0.05, cx + 0.14, cy + 0.22, col = "#2B5C8F", lwd = 1.5, len_mm = 1.8),
      arr(cx - 0.19, cy - 0.05, cx + 0.14, cy - 0.22, col = "#C0392B", lwd = 1.5, len_mm = 1.8),
      ellipse_grob(cx + 0.26, cy + 0.25, 0.11, 0.11, 0, fill = "#D4E6F1", col = "#2B5C8F", lwd = 1.5),
      ellipse_grob(cx + 0.26, cy - 0.25, 0.11, 0.11, 0, fill = "#FADBD8", col = "#C0392B", lwd = 1.5),
      txt(cx + 0.26, cy + 0.25, "E", size = 7.5, face = "bold", col = "#2B5C8F"),
      txt(cx + 0.26, cy - 0.25, "D", size = 7.5, face = "bold", col = "#C0392B")
    )
  }

  grobs <- gList(
    rectGrob(gp = gpar(fill = "white", col = NA)),

    # Panel headers & subtle vertical divider
    txt(0.20, 4.30, "a", size = 13.5, face = "bold", hjust = 0),
    txt(0.48, 4.30, "Time-Zero Design & Study Cohorts", size = 11.0, face = "bold", hjust = 0),
    txt(5.12, 4.30, "b", size = 13.5, face = "bold", hjust = 0),
    txt(5.42, 4.30, "Latent Phenotyping & Evaluations", size = 11.0, face = "bold", hjust = 0),
    segmentsGrob(x0 = ux(4.86), y0 = uy(0.22), x1 = ux(4.86), y1 = uy(4.36), gp = gpar(col = "#CBD5E1", lwd = 1.1)),

    # ---------------- PANEL (a): Time-Zero Gate & Compact Cohorts -------------
    rbox(0.20, 2.74, 2.52, 1.18, fill = "#EBF4FA", col = "#2B5C8F", lwd = 1.3, r = 0.12),
    icon_tube(0.56, 3.32),
    txt(0.86, 3.56, "ED Pre-Admission", size = 9.0, face = "bold", hjust = 0),
    txt(0.86, 3.26, "First-draw CBC & BMP", size = 8.0, col = "#334155", hjust = 0),
    txt(0.86, 2.98, "[-24 h, t = 0]", size = 8.0, face = "bold", col = "#2B5C8F", hjust = 0),

    # Vertical t = 0 Admission Gate
    segmentsGrob(x0 = ux(2.92), y0 = uy(2.66), x1 = ux(2.92), y1 = uy(4.00), gp = gpar(col = "#C0392B", lwd = 1.7, lty = "dashed")),
    txt(2.92, 4.10, "t = 0 Admission", size = 7.6, face = "bold", col = "#C0392B"),

    # Blocked Inpatient Care box
    rbox(3.10, 2.74, 1.56, 1.18, fill = "#F8FAFC", col = "#CBD5E1", lwd = 1.1, lty = "dashed", r = 0.12),
    txt(3.88, 3.58, "Inpatient Care", size = 8.3, face = "bold", col = "#64748B"),
    txt(3.88, 3.25, "Surgery / Fluids /\nTransfusion", size = 7.4, col = "#64748B"),
    txt(3.88, 2.90, "Excluded from LPA", size = 7.2, face = "bold", col = "#C0392B"),

    # Arrow down to Analytic Cohort
    arr(1.46, 2.72, 1.46, 2.26, col = "#1E293B", lwd = 1.5),
    rbox(0.20, 1.44, 2.52, 0.78, fill = "#FFFFFF", col = "#1E293B", lwd = 1.3, r = 0.14),
    txt(1.46, 1.95, "Hip Fracture Cohort", size = 9.0, face = "bold"),
    txt(1.46, 1.65, "Age \u2265 65 y  (n = 2,488)", size = 8.1, col = "#334155"),

    # Split into Derivation & Validation
    arr(0.92, 1.42, 0.82, 0.98, col = "#2B5C8F", lwd = 1.4),
    arr(2.00, 1.42, 2.12, 0.98, col = "#00897B", lwd = 1.4),

    rbox(0.18, 0.24, 1.28, 0.70, fill = "#EBF4FA", col = "#2B5C8F", lwd = 1.2, r = 0.14),
    txt(0.82, 0.70, "Derivation", size = 8.3, face = "bold", col = "#2B5C8F"),
    txt(0.82, 0.43, "\u22642016 (n = 1,312)", size = 7.4),

    rbox(1.54, 0.24, 1.30, 0.70, fill = "#E6F5F3", col = "#00897B", lwd = 1.2, r = 0.14),
    txt(2.19, 0.70, "Validation", size = 8.3, face = "bold", col = "#00897B"),
    txt(2.19, 0.43, "\u22652017 (n = 956)", size = 7.4),

    # Cross-Fracture Generalization Pill
    rbox(3.10, 0.24, 1.56, 1.98, fill = "#F8FAFC", col = "#475569", lwd = 1.15, r = 0.10),
    txt(3.88, 1.85, "Cross-Fracture\nCohort", size = 8.5, face = "bold"),
    txt(3.88, 1.22, "Non-hip major\nfractures \u2265 65 y", size = 7.6, col = "#334155"),
    txt(3.88, 0.72, "n = 5,508", size = 8.4, face = "bold", col = "#2B5C8F"),
    txt(3.88, 0.44, "(Frozen model)", size = 7.2, col = "#64748B"),

    # ---------------- PANEL (b): Dashed Latent Engine + 3 Pastel Cards --------
    rbox(5.12, 0.24, 4.66, 3.70, fill = "#FCFDFF", col = "#334155", lwd = 1.15, lty = "dashed", r = 0.08),

    # 3D Biomarker Axis Icon (left inside dashed box)
    arr(5.82, 2.42, 5.82, 3.08, col = "#111111", lwd = 1.4, len_mm = 1.7),
    arr(5.82, 2.42, 6.44, 2.42, col = "#111111", lwd = 1.4, len_mm = 1.7),
    arr(5.82, 2.42, 5.46, 2.02, col = "#111111", lwd = 1.4, len_mm = 1.7),
    txt(5.82, 3.22, "Hb / RDW", size = 7.5, face = "bold"),
    txt(6.30, 2.20, "Cr / K / HCO\u2083", size = 7.4, face = "bold"),
    txt(5.86, 1.80, "WBC / Glu / Plt / Na", size = 7.3, face = "bold"),
    ellipse_grob(6.08, 2.74, 0.070, 0.070, 0, fill = "#FDEBD0", col = "#111111", lwd = 1.1),
    txt(5.96, 1.42, "9 Routine ED Labs", size = 8.4, face = "bold"),

    # Separate Anticoagulation Pill
    rbox(5.30, 0.46, 1.42, 0.58, fill = "#FEF9E7", col = "#B7791F", lwd = 1.05, r = 0.16),
    txt(6.01, 0.75, "Anticoagulation\n(Separate Covariate)", size = 7.2, face = "bold", col = "#7D5A10"),

    # Arrow -> Diagonal GMM Module
    arr(6.58, 2.38, 6.94, 2.38, col = "#111111", lwd = 1.5),
    rbox(6.98, 1.92, 1.18, 0.92, fill = "#CFE2F3", col = "#3B66C4", lwd = 1.25, r = 0.15),
    txt(7.57, 2.52, "Diagonal GMM", size = 8.5, face = "bold"),
    txt(7.57, 2.22, "(VVI, K = 3)", size = 8.0, col = "#1E293B"),

    # Arrow -> Host Subphenotype Manifold
    arr(8.18, 2.38, 8.48, 2.38, col = "#111111", lwd = 1.5),
    manifold_grob(9.08, 2.38, 1.14, 1.75, fill = "#EEF3FA", col = "#5A6B82"),
    txt(9.08, 3.52, "Host Subphenotypes", size = 8.4, face = "bold"),

    # Three Gaussian ellipses on the manifold
    ellipse_grob(8.92, 2.82, 0.17, 0.11, 14, fill = "#D4E6F1", col = "#2B5C8F", lwd = 1.5),
    ellipse_grob(9.22, 2.36, 0.17, 0.11, -10, fill = "#D1F2EB", col = "#00897B", lwd = 1.5),
    ellipse_grob(8.98, 1.88, 0.17, 0.11, 12, fill = "#FADBD8", col = "#C0392B", lwd = 1.5),

    # Subphenotype legend strip inside bottom-right of dashed box
    ellipse_grob(7.10, 1.16, 0.055, 0.055, 0, fill = "#2B5C8F", col = NA, lwd = 0),
    txt(7.22, 1.16, "P1: Homeostatic (38.9%)", size = 7.6, face = "bold", col = "#2B5C8F", hjust = 0),
    ellipse_grob(7.10, 0.86, 0.055, 0.055, 0, fill = "#00897B", col = NA, lwd = 0),
    txt(7.22, 0.86, "P2: Stress-Inflammatory (35.3%)", size = 7.6, face = "bold", col = "#00897B", hjust = 0),
    ellipse_grob(7.10, 0.56, 0.055, 0.055, 0, fill = "#C0392B", col = NA, lwd = 0),
    txt(7.22, 0.56, "P3: Cardiorenal-Anemic (25.8%)", size = 7.6, face = "bold", col = "#C0392B", hjust = 0),

    # Orthogonal 3-Way Arrow Tree -> 3 Pastel Cards + Vector Icons
    segmentsGrob(x0 = ux(9.78), y0 = uy(2.09), x1 = ux(10.06), y1 = uy(2.09), gp = gpar(col = "#111111", lwd = 1.55)),
    segmentsGrob(x0 = ux(10.06), y0 = uy(0.80), x1 = ux(10.06), y1 = uy(3.38), gp = gpar(col = "#111111", lwd = 1.55)),
    arr(10.06, 3.38, 10.34, 3.38, col = "#111111", lwd = 1.55),
    arr(10.06, 2.09, 10.34, 2.09, col = "#111111", lwd = 1.55),
    arr(10.06, 0.80, 10.34, 0.80, col = "#111111", lwd = 1.55),

    # Card 1 (Warm Cream)
    rbox(10.36, 2.88, 2.12, 1.00, fill = "#FDF2D0", col = NA, lwd = 0, r = 0.16),
    txt(10.52, 3.38, "Aim 1: Prognostic &\nTemporal Validation", size = 9.0, face = "bold", hjust = 0, lh = 1.15),
    icon_surv(13.00, 3.38),

    # Card 2 (Soft Sage Green)
    rbox(10.36, 1.59, 2.12, 1.00, fill = "#E4F2E2", col = NA, lwd = 0, r = 0.16),
    txt(10.52, 2.09, "Aim 2: Incident AKI &\nDay-7 Landmark", size = 9.0, face = "bold", hjust = 0, lh = 1.15),
    icon_kidney(13.00, 2.09),

    # Card 3 (Soft Periwinkle Blue)
    rbox(10.36, 0.30, 2.12, 1.00, fill = "#DCE6F7", col = NA, lwd = 0, r = 0.16),
    txt(10.52, 0.80, "Aim 3: Surgical-Timing\nTarget Trial (CCW)", size = 9.0, face = "bold", hjust = 0, lh = 1.15),
    icon_ccw(13.00, 0.80)
  )
  gTree(children = grobs)
}

p_fig1 <- wrap_elements(full = draw_fig1_grob())
save_fig(p_fig1, "Fig1_Study_Design_and_Cohort_Flow", 13.6, 4.6, also_svg = TRUE)

# ==============================================================================
# FIGURE 2: Host Biomarker Profiles, Unmodeled Assays, Validation KM & Forest
# ==============================================================================
var_labels <- c(
  "hb" = "Hemoglobin", "log_rdw" = "RDW (log)", "log_plt" = "Platelets (log)",
  "log_wbc" = "WBC (log)", "log_cr" = "Creatinine (log)", "log_glu" = "Glucose (log)",
  "sodium" = "Sodium", "potassium" = "Potassium", "bicarbonate" = "Bicarbonate"
)
zm_df <- read.csv(file.path(V5_T, "v5_lpa_zmeans.csv"), stringsAsFactors = FALSE) %>%
  mutate(
    var_lbl = factor(var_labels[variable], levels = unname(var_labels)),
    pheno_lbl = factor(LBL_P[phenotype], levels = unname(LBL_P))
  )

p2a <- ggplot(zm_df, aes(x = var_lbl, y = mean_z, color = pheno_lbl, group = pheno_lbl)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#94A3B8", linewidth = 0.45) +
  geom_line(linewidth = 0.9) +
  geom_errorbar(aes(ymin = mean_z - 1.96 * se_z, ymax = mean_z + 1.96 * se_z), width = 0.16, linewidth = 0.55) +
  geom_point(size = 2.3) +
  scale_color_manual(values = PAL_P_FULL, name = NULL) +
  labs(x = NULL, y = "Standardized Z-Score") +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 26, hjust = 1, face = "bold"), legend.position = "top")

ext_df <- read.csv(file.path(V5_T, "SuppTable_S3_Unmodeled_Biomarkers.csv"), stringsAsFactors = FALSE)
ext_names <- c(
  "ext_albumin"    = "Albumin (g/dL)",
  "ext_lactate"    = "Lactate (mmol/L)",
  "ext_troponin_t" = "Troponin T (ng/mL)",
  "ext_ntprobnp"   = "NT-proBNP (pg/mL)"
)
ext_df$panel_lbl <- factor(ext_names[ext_df$marker], levels = unname(ext_names))
ext_df$p_f <- factor(ext_df$p_f, levels = c("P1", "P2", "P3"))

p2b <- ggplot(ext_df, aes(x = p_f, y = med, fill = p_f)) +
  geom_col(width = 0.54, alpha = 0.88) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.18, linewidth = 0.55, color = "#111111") +
  scale_fill_manual(values = PAL_P, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  facet_wrap(~ panel_lbl, scales = "free_y", nrow = 2) +
  labs(x = NULL, y = "Median [IQR]") +
  theme_nc()

km_df <- read.csv(file.path(V5_T, "v5_km_curves.csv"), stringsAsFactors = FALSE) %>%
  filter(split == "Post-2016 Validation") %>%
  mutate(pheno_lbl = factor(LBL_P[p_f], levels = unname(LBL_P)))
nr_df <- read.csv(file.path(V5_T, "v5_val_nrisk.csv"), stringsAsFactors = FALSE) %>%
  mutate(p_f = factor(p_f, levels = c("P3", "P2", "P1")))

p2c_top <- ggplot(km_df, aes(x = time, y = surv, color = pheno_lbl, fill = pheno_lbl)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.13, color = NA) +
  geom_step(linewidth = 0.9) +
  geom_vline(xintercept = 90, linetype = "dotted", color = "#64748B", linewidth = 0.5) +
  scale_color_manual(values = PAL_P_FULL, guide = "none") +
  scale_fill_manual(values = PAL_P_FULL, guide = "none") +
  scale_x_continuous(breaks = c(0, 90, 180, 270, 365), limits = c(0, 365)) +
  scale_y_continuous(limits = c(45, 100), labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Survival Probability") +
  theme_nc() +
  theme(plot.margin = margin(6, 10, 0, 6), axis.text.x = element_blank(), axis.ticks.x = element_blank())

p2c_bot <- ggplot(nr_df, aes(x = time, y = p_f, label = n_risk, color = p_f)) +
  geom_text(size = 2.65, fontface = "bold") +
  scale_color_manual(values = PAL_P, guide = "none") +
  scale_x_continuous(breaks = c(0, 90, 180, 270, 365), limits = c(0, 365)) +
  labs(x = "Days from Admission", y = "At risk") +
  theme_nc() +
  theme(panel.grid = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(0, 10, 6, 6), axis.title.y = element_text(size = 8))

t2_v5 <- read.csv(file.path(V5_T, "Table2_v5_Validation.csv"), stringsAsFactors = FALSE)
aki_v5 <- read.csv(file.path(V5_T, "SuppTable_S4b_KDIGO_AKI_Sensitivity.csv"), stringsAsFactors = FALSE)
for_df <- data.frame(
  endpoint = factor(rep(c("Full Hip Cohort", "Post-2016 Validation", "Non-Hip Fractures", "7-Day KDIGO AKI"), each = 2),
                    levels = rev(c("Full Hip Cohort", "Post-2016 Validation", "Non-Hip Fractures", "7-Day KDIGO AKI"))),
  contrast = rep(c("P2 vs P1", "P3 vs P1"), 4),
  est = c(1.02, 2.18, 1.41, 2.12, 1.80, 2.94, aki_v5$p2_est[1], aki_v5$p3_est[1]),
  lo  = c(0.77, 1.69, 0.92, 1.42, 1.43, 2.36, aki_v5$p2_lo[1], aki_v5$p3_lo[1]),
  hi  = c(1.34, 2.80, 2.16, 3.16, 2.26, 3.66, aki_v5$p2_hi[1], aki_v5$p3_hi[1])
)

p2d <- ggplot(for_df, aes(x = est, y = endpoint, color = contrast)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#64748B", linewidth = 0.5) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.22, position = position_dodge(width = 0.50), linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.50), size = 2.6) +
  scale_color_manual(values = c("P2 vs P1" = "#00897B", "P3 vs P1" = "#C0392B"), name = NULL) +
  scale_x_continuous(limits = c(0.65, 3.9), breaks = c(1, 2, 3)) +
  labs(x = "Adjusted HR / sHR (95% CI)", y = NULL) +
  theme_nc() +
  theme(legend.position = "top")

p2c_comb <- (p2c_top / p2c_bot) + plot_layout(heights = c(3.8, 1.05))
fig2 <- ((p2a | p2b) / (p2c_comb | p2d)) +
  plot_layout(heights = c(1.0, 1.05)) +
  plot_annotation(tag_levels = list(c("a", "b", "c", "", "d")))
save_fig(fig2, "Fig2_Host_Subphenotypes_and_Survival", 11.8, 8.2)

# ==============================================================================
# FIGURE 3: Calibration, DCA, Comorbidity Cross-Tab & Secular Shift (Visual-First)
# ==============================================================================
cal_df <- read.csv(file.path(V5_T, "v5_calibration_90d.csv"), stringsAsFactors = FALSE) %>%
  mutate(m_short = recode(model,
    "Clinical + Comorbidity LCA" = "Baseline (Clinical + LCA)",
    "Clinical + Comorbidity LCA + Subphenotypes" = "+ Host Subphenotypes"
  ))

p3a <- ggplot(cal_df, aes(x = pred * 100, y = obs * 100, color = m_short)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#94A3B8", linewidth = 0.5) +
  geom_errorbar(aes(ymin = lo * 100, ymax = hi * 100), width = 0.9, alpha = 0.55, linewidth = 0.45) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.3) +
  scale_color_manual(values = c("Baseline (Clinical + LCA)" = "#64748B", "+ Host Subphenotypes" = "#C0392B"), name = NULL) +
  coord_cartesian(xlim = c(0, 52), ylim = c(0, 52)) +
  labs(x = "Predicted 90-Day Mortality (%)", y = "Observed 90-Day Mortality (%)") +
  theme_nc() +
  theme(legend.position = "top")

dca_df <- read.csv(file.path(V5_T, "v5_dca_90d.csv"), stringsAsFactors = FALSE)
dca_long <- dca_df %>%
  pivot_longer(cols = c("Treat_All", "Treat_None", "Baseline", "Baseline_LPA"), names_to = "curve", values_to = "nb") %>%
  mutate(curve_lbl = factor(recode(curve,
    "Baseline_LPA" = "+ Host Subphenotypes",
    "Baseline"     = "Baseline (Clinical + LCA)",
    "Treat_All"    = "Intervene All", "Treat_None" = "Intervene None"
  ), levels = c("+ Host Subphenotypes", "Baseline (Clinical + LCA)", "Intervene All", "Intervene None")))

p3b <- ggplot(dca_long, aes(x = pt * 100, y = nb, color = curve_lbl, linetype = curve_lbl)) +
  geom_line(linewidth = 0.85) +
  scale_color_manual(values = c("+ Host Subphenotypes" = "#C0392B", "Baseline (Clinical + LCA)" = "#2B5C8F",
                                "Intervene All" = "#94A3B8", "Intervene None" = "#CBD5E1"), name = NULL) +
  scale_linetype_manual(values = c("solid", "solid", "dashed", "dotted"), name = NULL) +
  coord_cartesian(xlim = c(4, 40), ylim = c(-0.02, 0.15)) +
  labs(x = "Threshold Probability (%)", y = "Net Benefit") +
  theme_nc() +
  theme(legend.position = "top") +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))

ll_df <- read.csv(file.path(V5_T, "v5_lca_vs_lpa_crosstab.csv"), stringsAsFactors = FALSE) %>%
  mutate(
    lca_lbl = factor(c("LCA1" = "LCA1", "LCA2" = "LCA2", "LCA3" = "LCA3")[lca_class], levels = c("LCA3", "LCA2", "LCA1")),
    p_lbl = factor(c("P1" = "P1", "P2" = "P2", "P3" = "P3")[lpa_class], levels = c("P1", "P2", "P3"))
  )

p3c <- ggplot(ll_df, aes(x = p_lbl, y = lca_lbl, fill = mort_90d_pct)) +
  geom_tile(color = "white", linewidth = 1.4) +
  geom_text(aes(label = sprintf("%.1f%%", mort_90d_pct), color = mort_90d_pct > 22), size = 3.3, fontface = "bold") +
  scale_fill_gradient(low = "#EBF4FA", high = "#922B21", name = "90-d Mort (%)") +
  scale_color_manual(values = c("FALSE" = "#111111", "TRUE" = "white"), guide = "none") +
  labs(x = "Pre-Admission Subphenotype", y = "Comorbidity Class") +
  theme_nc() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())

sec_df <- read.csv(file.path(V5_T, "v5_strict_era_trend.csv"), stringsAsFactors = FALSE) %>%
  mutate(
    strict_era_bin = factor(strict_era_bin, levels = c("2008-2010", "2011-2013", "2014-2016", "2017-2019", "2020-2022")),
    pheno_lbl = factor(c("P1" = "P1", "P2" = "P2", "P3" = "P3")[p_f], levels = c("P1", "P2", "P3"))
  )
era_m <- read.csv(file.path(V5_T, "v5_strict_era_mortality.csv"), stringsAsFactors = FALSE) %>%
  mutate(strict_era_bin = factor(strict_era_bin, levels = c("2008-2010", "2011-2013", "2014-2016", "2017-2019", "2020-2022")))

p3d <- ggplot() +
  geom_col(data = sec_df, aes(x = strict_era_bin, y = pct, fill = pheno_lbl), width = 0.60, alpha = 0.88) +
  geom_line(data = era_m, aes(x = strict_era_bin, y = mort_1y_pct * 2.0, group = 1), color = "#111111", linewidth = 0.95) +
  geom_point(data = era_m, aes(x = strict_era_bin, y = mort_1y_pct * 2.0), color = "#111111", size = 2.4) +
  scale_fill_manual(values = PAL_P, name = NULL) +
  scale_y_continuous(name = "Subphenotype Share (%)",
                     sec.axis = sec_axis(~ . / 2.0, name = "1-Year Mortality (%)")) +
  labs(x = "Admission Period") +
  theme_nc() +
  theme(legend.position = "top", axis.text.x = element_text(angle = 18, hjust = 1))

fig3 <- ((p3a | p3b) / (p3c | p3d)) + plot_annotation(tag_levels = "a")
save_fig(fig3, "Fig3_Calibration_DCA_Comorbidity_and_Era", 11.4, 7.8)

# ==============================================================================
# FIGURE 4: Competing-Risk AKI, Day-7 Landmark, Mortality Matrix & Disposition
# ==============================================================================
aki_s <- read.csv(file.path(V5_T, "SuppTable_S4b_KDIGO_AKI_Sensitivity.csv"), stringsAsFactors = FALSE)
def_lbls <- c(
  "Primary KDIGO",
  "Relative-Only (ED Cr)",
  "Incident Relative (Prior Cr)",
  "Incident KDIGO (Prior Cr)"
)
aki_plot_df <- rbind(
  data.frame(def_short = def_lbls, contrast = "P2 vs P1", est = aki_s$p2_est, lo = aki_s$p2_lo, hi = aki_s$p2_hi),
  data.frame(def_short = def_lbls, contrast = "P3 vs P1", est = aki_s$p3_est, lo = aki_s$p3_lo, hi = aki_s$p3_hi)
)
aki_plot_df$def_short <- factor(aki_plot_df$def_short, levels = rev(def_lbls))

p4a <- ggplot(aki_plot_df, aes(x = est, y = def_short, color = contrast)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#64748B", linewidth = 0.5) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.22, position = position_dodge(width = 0.52), linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.52), size = 2.6) +
  scale_color_manual(values = c("P2 vs P1" = "#00897B", "P3 vs P1" = "#C0392B"), name = NULL) +
  scale_x_continuous(limits = c(0.65, 4.4), breaks = c(1, 2, 3, 4)) +
  labs(x = "Fine-Gray sHR (95% CI)", y = NULL) +
  theme_nc() +
  theme(legend.position = "top")

lm7_km <- read.csv(file.path(V5_T, "v5_day7_landmark_km.csv"), stringsAsFactors = FALSE) %>%
  mutate(
    pheno = substr(joint_grp, 1, 2),
    aki_st = ifelse(grepl("AKI\\+", joint_grp), "7-d AKI", "No AKI"),
    pheno_lbl = factor(c("P1" = "P1", "P2" = "P2", "P3" = "P3")[pheno], levels = c("P1", "P2", "P3"))
  )

p4b <- ggplot(lm7_km, aes(x = day_from_admit, y = surv, color = pheno_lbl, linetype = aki_st)) +
  geom_step(linewidth = 0.85) +
  scale_color_manual(values = PAL_P, name = NULL) +
  scale_linetype_manual(values = c("No AKI" = "solid", "7-d AKI" = "dashed"), name = NULL) +
  scale_x_continuous(breaks = c(7, 30, 60, 90), limits = c(7, 90)) +
  scale_y_continuous(limits = c(64, 100), labels = function(x) paste0(x, "%")) +
  labs(x = "Days from Admission (Day-7 Landmark)", y = "Survival Probability") +
  theme_nc() +
  theme(legend.position = "top")

rmat_aki <- read.csv(file.path(V5_T, "v5_nonesrd_risk_matrix.csv"), stringsAsFactors = FALSE) %>%
  mutate(
    p_lbl = factor(LBL_P[p_f], levels = rev(unname(LBL_P))),
    aki_lbl = factor(ifelse(grepl("KDIGO", aki_grp), "7-Day AKI", "No AKI"), levels = c("No AKI", "7-Day AKI"))
  )

p4c <- ggplot(rmat_aki, aes(x = aki_lbl, y = p_lbl, fill = mort_90d_pct)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = sprintf("%.1f%%", mort_90d_pct), color = mort_90d_pct > 25), size = 3.4, fontface = "bold") +
  scale_fill_gradient(low = "#EBF4FA", high = "#922B21", name = "90-d Mort (%)") +
  scale_color_manual(values = c("FALSE" = "#111111", "TRUE" = "white"), guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_nc() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())

dis_df <- read.csv(file.path(V5_T, "v5_discharge_disposition.csv"), stringsAsFactors = FALSE) %>%
  mutate(
    p_lbl = factor(c("P1" = "P1", "P2" = "P2", "P3" = "P3")[p_f], levels = c("P1", "P2", "P3")),
    disch_short = factor(recode(disch_cat,
      "In-hospital death / hospice" = "Death / hospice",
      "Skilled nursing / long-term care" = "Skilled nursing / LTC",
      "Home / acute rehabilitation" = "Home / acute rehab"
    ), levels = c("Death / hospice", "Skilled nursing / LTC", "Home / acute rehab"))
  )

p4d <- ggplot(dis_df, aes(x = p_lbl, y = pct, fill = disch_short)) +
  geom_col(width = 0.54, color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(pct >= 8, sprintf("%.1f%%", pct), "")),
            position = position_stack(vjust = 0.5), size = 2.8, fontface = "bold",
            color = ifelse(dis_df$disch_short == "Death / hospice", "white", "#111111")) +
  scale_fill_manual(values = c("Death / hospice" = "#922B21", "Skilled nursing / LTC" = "#D5DBDB", "Home / acute rehab" = "#A9CCE3"), name = NULL) +
  scale_y_continuous(limits = c(0, 101), breaks = c(0, 25, 50, 75, 100)) +
  labs(x = "Pre-Admission Subphenotype", y = "Admissions (%)") +
  theme_nc() +
  theme(legend.position = "top")

fig4 <- ((p4a | p4b) / (p4c | p4d)) + plot_annotation(tag_levels = "a")
save_fig(fig4, "Fig4_Renal_Landmark_and_Disposition", 11.4, 7.8)

# ==============================================================================
# FIGURE 5: Exploratory Target Trial Emulation (Visual-First)
# ==============================================================================
bal_df <- read.csv(file.path(V5_T, "SuppTable_S6_CCW_IPCW_Balance_PerPhenotype.csv"), stringsAsFactors = FALSE) %>%
  pivot_longer(cols = c("smd_unweighted", "smd_stabilized_ipcw"), names_to = "wt_type", values_to = "smd") %>%
  mutate(
    wt_lbl = factor(ifelse(wt_type == "smd_unweighted", "Unweighted", "IPCW Weighted"),
                    levels = c("Unweighted", "IPCW Weighted")),
    p_lbl = factor(LBL_P[stratum], levels = unname(LBL_P)),
    covariate = factor(covariate, levels = rev(unique(covariate)))
  )

p5a <- ggplot(bal_df, aes(x = smd, y = covariate, color = wt_lbl, shape = wt_lbl)) +
  geom_vline(xintercept = 0.10, linetype = "dashed", color = "#C0392B", linewidth = 0.5) +
  geom_point(size = 2.1, alpha = 0.9) +
  facet_wrap(~ p_lbl, nrow = 1) +
  scale_color_manual(values = c("Unweighted" = "#94A3B8", "IPCW Weighted" = "#2B5C8F"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  labs(x = "Absolute Standardized Mean Difference", y = NULL) +
  theme_nc() +
  theme(legend.position = "top")

ccw_tab <- read.csv(file.path(V5_T, "Table3_v5_Exploratory_CCW.csv"), stringsAsFactors = FALSE) %>%
  mutate(
    stratum_clean = factor(gsub("^Day0_All_Eligible: |^Operated_Subgroup: ", "", stratum),
                           levels = rev(c("Overall", "P1: Homeostatic", "P2: Stress-Inflammatory", "P3: Cardiorenal-Anemic"))),
    cohort_lbl = factor(ifelse(cohort == "Day0_All_Eligible", "All-Eligible (n = 2,485)", "Operated (n = 2,216)"),
                        levels = c("All-Eligible (n = 2,485)", "Operated (n = 2,216)"))
  )

p5b <- ggplot(ccw_tab, aes(x = rd_30d_pct, y = stratum_clean, color = cohort_lbl)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#64748B", linewidth = 0.5) +
  geom_errorbar(aes(xmin = rd_30d_lo, xmax = rd_30d_hi), orientation = "y", width = 0.22, position = position_dodge(width = 0.50), linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.50), size = 2.6) +
  scale_color_manual(values = c("#C0392B", "#2B5C8F"), name = NULL) +
  scale_x_continuous(limits = c(-6, 22), breaks = c(-5, 0, 5, 10, 15, 20)) +
  labs(x = "30-Day Absolute Mortality Risk Difference (%)", y = NULL) +
  theme_nc() +
  theme(legend.position = "top")

ckm_df <- read.csv(file.path(V5_T, "v5_ccw_km_curves.csv"), stringsAsFactors = FALSE) %>%
  mutate(p_lbl = factor(LBL_P[p_f], levels = unname(LBL_P)))

p5c <- ggplot(ckm_df, aes(x = time, y = surv, color = strategy, fill = strategy)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, color = NA) +
  geom_step(linewidth = 0.9) +
  geom_vline(xintercept = 30, linetype = "dotted", color = "#64748B", linewidth = 0.45) +
  facet_wrap(~ p_lbl, nrow = 1) +
  scale_color_manual(values = c("Early Surgery (Day 0-1)" = "#2B5C8F", "Delayed / Deferred (>Day 1)" = "#C0392B"), name = NULL) +
  scale_fill_manual(values = c("Early Surgery (Day 0-1)" = "#2B5C8F", "Delayed / Deferred (>Day 1)" = "#C0392B"), name = NULL) +
  scale_x_continuous(breaks = c(0, 30, 60, 90), limits = c(0, 90)) +
  scale_y_continuous(limits = c(60, 100), labels = function(x) paste0(x, "%")) +
  labs(x = "Days from Admission", y = "Weighted Survival") +
  theme_nc() +
  theme(legend.position = "top")

fig5 <- (wrap_elements(p5a) / wrap_elements(p5b) / wrap_elements(p5c)) +
  plot_layout(heights = c(1.08, 0.80, 0.98)) +
  plot_annotation(tag_levels = "a")
save_fig(fig5, "Fig5_Exploratory_Target_Trial", 11.4, 8.8)

# ==============================================================================
# SUPPLEMENTARY FIGURE S1: Correlations, 20-Model Grid & BLRT (Visual-First)
# ==============================================================================
cor_m <- read.csv(file.path(V5_T, "SuppTable_S2a_Pairwise_Correlations_9Lab.csv"), row.names = 1, check.names = FALSE)
cor_long <- as.data.frame(as.table(as.matrix(cor_m)))
colnames(cor_long) <- c("v1", "v2", "r")
cor_long$v1_lbl <- factor(var_labels[as.character(cor_long$v1)], levels = unname(var_labels))
cor_long$v2_lbl <- factor(var_labels[as.character(cor_long$v2)], levels = rev(unname(var_labels)))

ps1a <- ggplot(cor_long, aes(x = v1_lbl, y = v2_lbl, fill = r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#2B5C8F", mid = "#FFFFFF", high = "#C0392B", midpoint = 0, limits = c(-1, 1), name = "r") +
  labs(x = NULL, y = NULL) +
  theme_nc() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7.8),
        axis.text.y = element_text(size = 7.8),
        axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.position = "right")

grid_df <- read.csv(file.path(V5_T, "SuppTable_S2b_LPA_20Model_Grid.csv"), stringsAsFactors = FALSE)
bic_long <- rbind(
  data.frame(K = grid_df$K, val = grid_df$BIC / 1000, series = paste0(grid_df$model, " (BIC)")),
  data.frame(K = grid_df$K[grid_df$model == "VVI"], val = grid_df$ICL[grid_df$model == "VVI"] / 1000, series = "VVI (ICL)")
)
bic_long$series <- factor(bic_long$series, levels = c("EEI (BIC)", "EEE (BIC)", "VVI (BIC)", "VVV (BIC)", "VVI (ICL)"))

ps1b <- ggplot(bic_long, aes(x = K, y = val, color = series, linetype = series, shape = series)) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.2) +
  annotate("point", x = 3, y = 32529.4 / 1000, shape = 21, size = 4.4, stroke = 1.2, color = "#C0392B", fill = NA) +
  scale_color_manual(values = c("EEI (BIC)" = "#94A3B8", "EEE (BIC)" = "#64748B", "VVI (BIC)" = "#C0392B", "VVV (BIC)" = "#2B5C8F", "VVI (ICL)" = "#111111"), name = NULL) +
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid", "dashed"), name = NULL) +
  scale_shape_manual(values = c(16, 16, 16, 16, 17), name = NULL) +
  scale_x_continuous(breaks = 2:6) +
  labs(x = "Latent Profiles (K)", y = "Criterion (\u00d710\u00b3)") +
  theme_nc() +
  theme(legend.position = "top") +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

blrt_df <- read.csv(file.path(V5_T, "SuppTable_S2c_VVI_BLRT.csv"), stringsAsFactors = FALSE) %>%
  mutate(transition = factor(transition, levels = transition))
ps1c <- ggplot(blrt_df, aes(x = transition, y = LRTS_obs)) +
  geom_col(fill = ifelse(blrt_df$transition %in% c("1 vs 2", "2 vs 3"), "#2B5C8F", "#CBD5E1"), width = 0.54) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Sequential K", y = "Bootstrap LRT Statistic") +
  theme_nc()

figs1 <- (ps1a | ps1b | ps1c) + plot_layout(widths = c(1.12, 1.08, 0.85)) + plot_annotation(tag_levels = "a")
save_fig(figs1, "FigS1_Model_Selection_and_Correlations", 12.0, 4.2)

# ==============================================================================
# SUPPLEMENTARY FIGURE S2: Non-Hip Generalization & 15-Arm Forest (Visual-First)
# ==============================================================================
km_nh <- read.csv(file.path(V5_T, "v5_km_curves.csv"), stringsAsFactors = FALSE) %>%
  filter(split == "Non-Hip Fractures") %>%
  mutate(pheno_lbl = factor(c("P1" = "P1", "P2" = "P2", "P3" = "P3")[p_f], levels = c("P1", "P2", "P3")))

ps2a <- ggplot(km_nh, aes(x = time, y = surv, color = pheno_lbl, fill = pheno_lbl)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.13, color = NA) +
  geom_step(linewidth = 0.9) +
  geom_vline(xintercept = 90, linetype = "dotted", color = "#64748B", linewidth = 0.5) +
  scale_color_manual(values = PAL_P, name = NULL) +
  scale_fill_manual(values = PAL_P, name = NULL) +
  scale_x_continuous(breaks = c(0, 90, 180, 270, 365), limits = c(0, 365)) +
  scale_y_continuous(limits = c(62, 100), labels = function(x) paste0(x, "%")) +
  labs(x = "Days from Admission", y = "Survival Probability") +
  theme_nc() +
  theme(legend.position = "top")

s15 <- read.csv(file.path(V5_T, "SuppTable_S4a_15Arm_Sensitivity.csv"), stringsAsFactors = FALSE) %>%
  mutate(
    short_arm = c(
      "1. Primary Full Cohort",
      "2. Chronic-Only CCI",
      "3. Excl. Dialysis / ESRD",
      "4. Pre-2017 Derivation",
      "5. Post-2016 Validation (Frozen)",
      "6. Validation (Era-Standardized)",
      "7. Validation Pre-Pandemic",
      "8. Validation 2019-2020 Boundary",
      "9. Validation COVID-19 Era",
      "10. Soft-Probability Weighted",
      "11. Operated Subgroup",
      "12. Diagonal VVI K = 2",
      "13. Full-Covariance VVV K = 2",
      "14. INR-Inclusive VVI K = 4",
      "15. Non-Hip Major Fractures"
    ),
    arm_f = factor(short_arm, levels = rev(short_arm))
  )

ps2b <- ggplot(s15, aes(x = hr, y = arm_f)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#64748B", linewidth = 0.5) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.24, color = "#C0392B", linewidth = 0.65) +
  geom_point(color = "#C0392B", size = 2.2) +
  scale_x_continuous(limits = c(0.6, 5.6), breaks = c(1, 2, 3, 4, 5)) +
  labs(x = "Adjusted 90-Day Mortality HR (95% CI)", y = NULL) +
  theme_nc()

figs2 <- (ps2a | ps2b) + plot_layout(widths = c(0.95, 1.15)) + plot_annotation(tag_levels = "a")
save_fig(figs2, "FigS2_CrossFracture_and_Sensitivity", 11.4, 5.2)

# ==============================================================================
# SUPPLEMENTARY FIGURE S3: 72-h Creatinine Trajectories & Operated CCW
# ==============================================================================
tr_df <- read.csv(file.path(V5_T, "v5_nonesrd_lcmm_trajectories.csv"), stringsAsFactors = FALSE) %>%
  mutate(t_lbl = factor(recode(traj_class,
    "T1" = "T1: Preserved Baseline (92.0%)",
    "T2" = "T2: Elevated Baseline (8.0%)"
  )))

ps3a <- ggplot(tr_df, aes(x = hour, y = pred_cr, color = t_lbl, fill = t_lbl)) +
  geom_ribbon(aes(ymin = lo_cr, ymax = hi_cr), alpha = 0.16, color = NA) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = c("#2B5C8F", "#C0392B"), name = NULL) +
  scale_fill_manual(values = c("#2B5C8F", "#C0392B"), name = NULL) +
  scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60, 72)) +
  labs(x = "Hours from Admission", y = "Serum Creatinine (mg/dL)") +
  theme_nc() +
  theme(legend.position = "top")

ckm_op <- read.csv(file.path(V5_T, "v5_ccw_km_curves_operated.csv"), stringsAsFactors = FALSE) %>%
  mutate(
    p_lbl = factor(LBL_P[p_f], levels = unname(LBL_P)),
    strategy = gsub(">=Day 2", "\u2265Day 2", strategy)
  )

ps3b <- ggplot(ckm_op, aes(x = time, y = surv, color = strategy, fill = strategy)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, color = NA) +
  geom_step(linewidth = 0.85) +
  geom_vline(xintercept = 30, linetype = "dotted", color = "#64748B", linewidth = 0.45) +
  facet_wrap(~ p_lbl, nrow = 1) +
  scale_color_manual(values = c("Early Surgery (Day 0-1)" = "#2B5C8F", "Delayed Surgery (\u2265Day 2)" = "#C0392B"), name = NULL) +
  scale_fill_manual(values = c("Early Surgery (Day 0-1)" = "#2B5C8F", "Delayed Surgery (\u2265Day 2)" = "#C0392B"), name = NULL) +
  scale_x_continuous(breaks = c(0, 30, 60, 90), limits = c(0, 90)) +
  scale_y_continuous(limits = c(55, 100), labels = function(x) paste0(x, "%")) +
  labs(x = "Days from Admission", y = "Weighted Survival") +
  theme_nc() +
  theme(legend.position = "top")

figs3 <- (ps3a / ps3b) + plot_layout(heights = c(0.92, 1.05)) + plot_annotation(tag_levels = "a")
save_fig(figs3, "FigS3_Renal_Trajectories_and_Operated_CCW", 11.2, 7.2)

cat("=== ALL 8 VISUAL-FIRST V5 FIGURES GENERATED SUCCESSFULLY ===\n")
