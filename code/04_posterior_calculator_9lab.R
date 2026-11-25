# Closed-Form 9-Biomarker Pre-Admission Host Subphenotype Posterior Calculator
predict_hip_subphenotype_9lab <- function(hb, rdw, platelets, wbc, creatinine, glucose, sodium, potassium, bicarbonate, param_csv) {
  p <- read.csv(param_csv, stringsAsFactors = FALSE)
  pi_k <- c(P1 = 0.3929, P2 = 0.3607, P3 = 0.2465)
  raw <- list(hb = hb, log_rdw = log(rdw), log_plt = log(platelets), log_wbc = log(wbc),
              log_cr = log(creatinine), log_glu = log(glucose), sodium = sodium,
              potassium = potassium, bicarbonate = bicarbonate)
  n <- length(hb)
  Z <- matrix(0, nrow = n, ncol = 9, dimnames = list(NULL, p$biomarker))
  for (j in seq_len(9)) {
    v <- p$biomarker[j]; x <- raw[[v]]
    x <- ifelse(is.na(x), p$dev_median[j], pmin(pmax(x, p$winsor_0p5[j]), p$winsor_99p5[j]))
    Z[, j] <- (x - p$dev_mean[j]) / p$dev_sd[j]
  }
  ll <- matrix(0, nrow = n, ncol = 3, dimnames = list(NULL, c("P1","P2","P3")))
  for (k in 1:3) {
    mu <- p[[paste0("P", k, "_mean_z")]]; sdv <- p[[paste0("P", k, "_sd_z")]]
    ll[, k] <- log(pi_k[k]) - 0.5 * rowSums(sweep(sweep(Z, 2, mu, "-")^2, 2, sdv^2, "/") + rep(log(2 * pi * sdv^2), each = n))
  }
  mx <- apply(ll, 1, max); pr <- exp(ll - mx); pr <- pr / rowSums(pr)
  data.frame(prob_P1 = round(pr[,1], 4), prob_P2 = round(pr[,2], 4), prob_P3 = round(pr[,3], 4),
             assigned_subphenotype = colnames(pr)[max.col(pr)])
}

# ------------------------------------------------------------------------------
# Runnable Demonstration (executes when called via `Rscript code/04_posterior_calculator_9lab.R`)
# ------------------------------------------------------------------------------
if (sys.nframe() == 0) {
  args_all <- commandArgs(trailingOnly = FALSE)
  script_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
  s_dir <- if (length(script_arg) > 0) dirname(normalizePath(script_arg)) else getwd()
  param_path <- file.path(s_dir, "..", "parameters", "SuppTable_S5_Exact_9Lab_Posterior_Parameters.csv")
  if (!file.exists(param_path)) param_path <- "SuppTable_S5_Exact_9Lab_Posterior_Parameters.csv"

  demo_patients <- predict_hip_subphenotype_9lab(
    hb          = c(12.4, 12.1,  9.8),
    rdw         = c(13.6, 13.4, 16.2),
    platelets   = c(210,  245,  185),
    wbc         = c(9.2,  14.2,  9.6),
    creatinine  = c(0.80, 0.90, 1.85),
    glucose     = c(112,  168,  128),
    sodium      = c(140,  136,  138),
    potassium   = c(4.1,  4.4,  5.0),
    bicarbonate = c(26,   24,   22),
    param_csv   = param_path
  )
  rownames(demo_patients) <- c("Example_Homeostatic", "Example_StressInflam", "Example_Cardiorenal")
  cat("=== Closed-Form 9-Biomarker Pre-Admission Subphenotype Calculator Demo ===\n")
  print(demo_patients)
}
