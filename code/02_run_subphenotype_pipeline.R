#!/usr/bin/env Rscript
# ==============================================================================
# Paper 1 v5 (Nature Communications Final Submission Package) — Complete R Pipeline
#   1. Copies base LPA fit tables + adds Era-Standardized Validation Sensitivity (Arm 15)
#      quantifying the impact of BIDMC analyzer/reference-interval shifts
#   2. Uses VanderWeele & Ding (2017) common-outcome (>15% incidence) HR-to-RR
#      conversion for all 90-day mortality E-values
#   3. Evaluates STRICTLY INCIDENT post-admission AKI anchored to prior [-365d,-7d] Cr
#      AFTER excluding the 38 patients with arrival AKI (n = 1,029; zero circularity!)
#   4. Places the 1 genuine overnight ED-to-OR patient (overnight ED-to-OR Day-0 arrival) in Early Day 0-1
#      and excludes only the 3 pre-episode date anomalies (n = 2,485 All-Elig; n = 2,216 Op)
# ==============================================================================

suppressPackageStartupMessages({
  library(mclust)
  library(poLCA)
  library(lcmm)
  library(survival)
  library(cmprsk)
  library(dplyr)
  library(tidyr)
})

set.seed(20261124)

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
SCRIPT_DIR <- if (length(script_arg) > 0) dirname(normalizePath(script_arg)) else getwd()
BASE_DIR   <- normalizePath(file.path(SCRIPT_DIR, ".."))
OUT_T      <- file.path(BASE_DIR, "tables")

d_raw  <- read.csv(file.path(OUT_T, "hip65_v5_preadmit_master.csv"), stringsAsFactors = FALSE)
nh_raw <- read.csv(file.path(OUT_T, "nonhip_v5_preadmit.csv"), stringsAsFactors = FALSE)

raw9 <- c("hb","rdw","platelets","wbc","creatinine","glucose","sodium","potassium","bicarbonate")
d_raw$n_obs9  <- rowSums(!is.na(d_raw[, raw9]))
nh_raw$n_obs9 <- rowSums(!is.na(nh_raw[, raw9]))

df <- d_raw[d_raw$n_obs9 >= 7, ]
nh <- nh_raw[nh_raw$n_obs9 >= 7, ]

prep_vars <- function(x) {
  x$log_rdw <- log(x$rdw); x$log_plt <- log(x$platelets); x$log_wbc <- log(x$wbc)
  x$log_cr  <- log(x$creatinine); x$log_glu <- log(x$glucose)
  x$sex_m   <- ifelse(x$gender == "M", 1L, 0L)
  x
}
df <- prep_vars(df); nh <- prep_vars(nh)

df$hx_dm_any <- ifelse(df$hx_dm_uncomp == 1 | df$hx_dm_comp == 1, 1L, 0L)
if (!"hx_dm_any" %in% names(nh)) nh$hx_dm_any <- ifelse(nh$hx_dm_uncomp == 1 | nh$hx_dm_comp == 1, 1L, 0L)

df$anticoag_flag <- ifelse((!is.na(df$inr) & df$inr >= 1.5) | df$dx_longterm_anticoag == 1 |
                             df$prior_hosp_doac == 1 | df$prior_hosp_warfarin == 1, 1L, 0L)
df$med_srv <- ifelse(df$first_service %in% c("MED","OMED","CMED","NMED"), 1L, 0L)
df$t30 <- pmin(df$surv_days_365, 30.0); df$e30 <- ifelse(df$death_30d == 1 & df$surv_days_365 <= 30.0, 1L, 0L)
df$t90 <- pmin(df$surv_days_365, 90.0); df$e90 <- ifelse(df$death_90d == 1 & df$surv_days_365 <= 90.0, 1L, 0L)
nh$t90 <- pmin(nh$surv_days_365, 90.0); nh$e90 <- ifelse(nh$death_90d == 1 & nh$surv_days_365 <= 90.0, 1L, 0L)

dev_idx <- df$strict_temporal_split == "strict_pre2017"
val_idx <- df$strict_temporal_split == "strict_post2016"

vars9 <- c("hb","log_rdw","log_plt","log_wbc","log_cr","log_glu","sodium","potassium","bicarbonate")
Xd_raw <- as.matrix(df[dev_idx, vars9]); Xa_raw <- as.matrix(df[, vars9]); Xn_raw <- as.matrix(nh[, vars9])
Xv_raw <- as.matrix(df[val_idx, vars9])
Xd <- matrix(NA_real_, nrow(Xd_raw), 9, dimnames = list(NULL, vars9))
Xa <- matrix(NA_real_, nrow(Xa_raw), 9, dimnames = list(NULL, vars9))
Xn <- matrix(NA_real_, nrow(Xn_raw), 9, dimnames = list(NULL, vars9))
Xv_erastd <- matrix(NA_real_, nrow(Xv_raw), 9, dimnames = list(NULL, vars9))

for (v in vars9) {
  dv <- Xd_raw[, v]
  lo <- as.numeric(quantile(dv, 0.005, na.rm = TRUE)); hi <- as.numeric(quantile(dv, 0.995, na.rm = TRUE))
  md <- as.numeric(median(dv, na.rm = TRUE))
  dw <- pmin(pmax(dv, lo), hi); dw[is.na(dw)] <- md
  mu <- mean(dw); sdv <- sd(dw)
  aw <- pmin(pmax(Xa_raw[, v], lo), hi); aw[is.na(aw)] <- md
  nw <- pmin(pmax(Xn_raw[, v], lo), hi); nw[is.na(nw)] <- md
  Xd[, v] <- (dw - mu) / sdv; Xa[, v] <- (aw - mu) / sdv; Xn[, v] <- (nw - mu) / sdv

  # Era-standardized validation Z-scores (removing inter-era analyzer/reference shift!)
  vv <- Xv_raw[, v]
  vlo <- as.numeric(quantile(vv, 0.005, na.rm = TRUE)); vhi <- as.numeric(quantile(vv, 0.995, na.rm = TRUE))
  vmd <- as.numeric(median(vv, na.rm = TRUE))
  vw <- pmin(pmax(vv, vlo), vhi); vw[is.na(vw)] <- vmd
  Xv_erastd[, v] <- (vw - mean(vw)) / sd(vw)
}

set.seed(20261124)
mod_vvi3 <- Mclust(Xd, G = 3, modelNames = "VVI", verbose = FALSE)
zm3 <- t(mod_vvi3$parameters$mean)
r_p3 <- which.max(zm3[, "log_cr"] - zm3[, "hb"]); rem3 <- setdiff(1:3, r_p3)
r_p2 <- rem3[which.max(zm3[rem3, "log_wbc"] + zm3[rem3, "log_glu"])]; r_p1 <- setdiff(rem3, r_p2)[1]
raw_to_p <- integer(3); raw_to_p[r_p1] <- 1L; raw_to_p[r_p2] <- 2L; raw_to_p[r_p3] <- 3L
p_to_raw <- c(r_p1, r_p2, r_p3)

pr_all <- predict(mod_vvi3, newdata = Xa)
z_all  <- pr_all$z[, p_to_raw, drop = FALSE]
df$lpa_class <- raw_to_p[pr_all$classification]
df$p_f <- factor(paste0("P", df$lpa_class), levels = paste0("P", 1:3))
df$prob_p1 <- z_all[, 1]; df$prob_p2 <- z_all[, 2]; df$prob_p3 <- z_all[, 3]

pr_nh <- predict(mod_vvi3, newdata = Xn)
nh$p_f <- factor(paste0("P", raw_to_p[pr_nh$classification]), levels = paste0("P", 1:3))

# Check Era-Standardized Post-2016 Validation Assignment (Arm 15!)
pr_val_erastd <- predict(mod_vvi3, newdata = Xv_erastd)
val_p_erastd <- factor(paste0("P", raw_to_p[pr_val_erastd$classification]), levels = paste0("P", 1:3))
cat("Post-2016 Validation P1/P2/P3 % — Frozen Z-scores vs Era-Standardized Z-scores:\n")
print(round(100 * prop.table(table(df$p_f[val_idx])), 1))
print(round(100 * prop.table(table(val_p_erastd)), 1))

set.seed(20261124)
mod_vvi2 <- Mclust(Xd, G = 2, modelNames = "VVI", verbose = FALSE)
r_vvi2_hi <- which.max(t(mod_vvi2$parameters$mean)[, "log_cr"] - t(mod_vvi2$parameters$mean)[, "hb"])
df$vvi2_high <- ifelse(predict(mod_vvi2, newdata = Xa)$classification == r_vvi2_hi, 1L, 0L)

set.seed(20261124)
mod_vvv2 <- Mclust(Xd, G = 2, modelNames = "VVV", verbose = FALSE)
r_vvv2_hi <- which.max(t(mod_vvv2$parameters$mean)[, "log_cr"] - t(mod_vvv2$parameters$mean)[, "hb"])
df$vvv2_high <- ifelse(predict(mod_vvv2, newdata = Xa)$classification == r_vvv2_hi, 1L, 0L)

df$log_inr <- log(ifelse(is.na(df$inr), median(df$inr[dev_idx], na.rm=TRUE), df$inr))
inr_d <- df$log_inr[dev_idx]; inr_lo <- quantile(inr_d, 0.005); inr_hi <- quantile(inr_d, 0.995)
inr_dw <- pmin(pmax(inr_d, inr_lo), inr_hi); inr_mu <- mean(inr_dw); inr_sd <- sd(inr_dw)
Xd10 <- cbind(Xd, log_inr = (inr_dw - inr_mu) / inr_sd)
Xa10 <- cbind(Xa, log_inr = (pmin(pmax(df$log_inr, inr_lo), inr_hi) - inr_mu) / inr_sd)
set.seed(20261124)
mod_inr4 <- Mclust(Xd10, G = 4, modelNames = "VVI", verbose = FALSE)
zm4 <- t(mod_inr4$parameters$mean)
r4_vka <- which.max(zm4[, "log_inr"]); rem4 <- setdiff(1:4, r4_vka)
r4_cr  <- rem4[which.max(zm4[rem4, "log_cr"] - zm4[rem4, "hb"])]; rem4 <- setdiff(rem4, r4_cr)
r4_inf <- rem4[which.max(zm4[rem4, "log_wbc"] + zm4[rem4, "log_glu"])]; r4_pr <- setdiff(rem4, r4_inf)[1]
map4 <- integer(4); map4[r4_pr] <- 1L; map4[r4_inf] <- 2L; map4[r4_cr] <- 3L; map4[r4_vka] <- 4L
df$inr4_f <- factor(paste0("K4_P", map4[predict(mod_inr4, newdata = Xa10)$classification]), levels = paste0("K4_P", 1:4))

# Comorbidity LCA
lca_cols <- c("hx_mi","hx_chf","hx_pvd","hx_cevd","hx_dementia","hx_copd","hx_dm_any","hx_renal")
lca_hip <- as.data.frame(lapply(df[, lca_cols], function(x) as.integer(x) + 1L))
lca_nh  <- as.data.frame(lapply(nh[, lca_cols], function(x) as.integer(ifelse(is.na(x), 0L, x)) + 1L))
f_lca <- as.formula(paste0("cbind(", paste(lca_cols, collapse = ","), ") ~ 1"))
set.seed(20261124)
mod_lca <- poLCA(f_lca, data = lca_hip[dev_idx, ], nrep = 15, nclass = 3, verbose = FALSE)
lca_ord <- order(tapply(df$cci_score[dev_idx], mod_lca$predclass, mean))
lca_map <- integer(3); lca_map[lca_ord] <- 1:3
pred_polca <- function(mod, nd, cols) {
  ll <- matrix(0, nrow(nd), length(mod$P))
  for (k in seq_along(mod$P)) {
    lp <- rep(log(pmax(mod$P[k], 1e-12)), nrow(nd))
    for (v in cols) lp <- lp + log(pmax(mod$probs[[v]][k, ], 1e-12)[nd[[v]]])
    ll[, k] <- lp
  }
  apply(ll, 1, which.max)
}
df$lca_f <- factor(paste0("LCA", lca_map[pred_polca(mod_lca, lca_hip, lca_cols)]), levels = paste0("LCA", 1:3))
nh$lca_f <- factor(paste0("LCA", lca_map[pred_polca(mod_lca, lca_nh,  lca_cols)]), levels = paste0("LCA", 1:3))

s_dev <- df[dev_idx, ]; s_val <- df[val_idx, ]; s_full <- df
s_val$p_erastd <- val_p_erastd

# Check C-index when Era-Standardized Z-scores are used in Post-2016 Validation!
c_m1_dev <- coxph(Surv(t90, e90) ~ admit_age + sex_m + fx_neck_all + fx_pertroch + cci_score + lca_f + anticoag_flag, data = s_dev)
c_m2_dev <- coxph(Surv(t90, e90) ~ admit_age + sex_m + fx_neck_all + fx_pertroch + cci_score + lca_f + anticoag_flag + p_f, data = s_dev)
c_m2_val <- coxph(Surv(t90, e90) ~ admit_age + sex_m + fx_neck_all + fx_pertroch + cci_score + lca_f + anticoag_flag + p_f, data = s_val)
c_m2_val_era <- coxph(Surv(t90, e90) ~ admit_age + sex_m + fx_neck_all + fx_pertroch + cci_score + lca_f + anticoag_flag + p_erastd, data = s_val)
c_m2_full <- coxph(Surv(t90, e90) ~ admit_age + sex_m + fx_neck_all + fx_pertroch + cci_score + lca_f + anticoag_flag + p_f, data = s_full)
c_m2_nh   <- coxph(Surv(t90, e90) ~ admit_age + sex_m + fx_spine + fx_pelvis + fx_upper + fx_lower + cci_score + lca_f + p_f, data = nh)

# ------------------------------------------------------------------------------
# 1. 15-Arm Mortality Sensitivity Table with VanderWeele Common-Outcome (>15%) E-values!
#    VanderWeele & Ding (Ann Intern Med 2017): RR_approx = (1 - 0.5^sqrt(HR)) / (1 - 0.5^sqrt(1/HR))
# ------------------------------------------------------------------------------
e_val_common_hr <- function(hr) {
  h <- ifelse(hr >= 1, hr, 1 / hr)
  rr <- (1 - 0.5^sqrt(h)) / (1 - 0.5^sqrt(1 / h))
  rr + sqrt(rr * (rr - 1))
}

extr_hr <- function(mod, term, lbl, n_v) {
  ci <- summary(mod)$conf.int[term, ]; pv <- summary(mod)$coefficients[term, "Pr(>|z|)"]
  data.frame(arm = lbl, n = n_v, hr = round(unname(ci[1]), 2), lo = round(unname(ci[3]), 2),
             hi = round(unname(ci[4]), 2), p_value = signif(pv, 3),
             e_value_common = round(unname(e_val_common_hr(ci[1])), 2), stringsAsFactors = FALSE)
}

s_soft <- do.call(rbind, lapply(1:3, function(k) {
  tmp <- s_full; tmp$p_f <- factor(paste0("P", k), levels=paste0("P",1:3)); tmp$w_s <- tmp[[paste0("prob_p", k)]]; tmp
}))
m_soft   <- coxph(Surv(t90, e90) ~ admit_age + sex_m + fx_neck_all + fx_pertroch + cci_score + lca_f + anticoag_flag + p_f, data = s_soft, weights = w_s, cluster = subject_id)
m_chron  <- coxph(Surv(t90, e90) ~ admit_age + sex_m + fx_neck_all + fx_pertroch + cci_chronic_only + lca_f + anticoag_flag + p_f, data = s_full)
m_noesrd <- coxph(Surv(t90, e90) ~ admit_age + sex_m + fx_neck_all + fx_pertroch + cci_score + lca_f + anticoag_flag + p_f, data = s_full[s_full$esrd_dialysis == 0, ])
m_op     <- coxph(Surv(t90, e90) ~ admit_age + sex_m + fx_neck_all + fx_pertroch + cci_score + lca_f + anticoag_flag + p_f, data = s_full[s_full$had_surgery == 1, ])
m_precov <- coxph(Surv(t90, e90) ~ admit_age + sex_m + fx_neck_all + fx_pertroch + cci_score + lca_f + anticoag_flag + p_f, data = s_val[s_val$val_pandemic_split == "val_pre_covid_2017_2019", ])
m_covbnd <- coxph(Surv(t90, e90) ~ admit_age + sex_m + fx_neck_all + fx_pertroch + cci_score + lca_f + anticoag_flag + p_f, data = s_val[s_val$val_pandemic_split == "val_boundary_2019_2020", ])
m_covid  <- coxph(Surv(t90, e90) ~ admit_age + sex_m + fx_neck_all + fx_pertroch + cci_score + lca_f + anticoag_flag + p_f, data = s_val[s_val$val_pandemic_split == "val_covid_2020_2022", ])
m_vvi2   <- coxph(Surv(t90, e90) ~ admit_age + sex_m + fx_neck_all + fx_pertroch + cci_score + lca_f + anticoag_flag + vvi2_high, data = s_full)
m_vvv2   <- coxph(Surv(t90, e90) ~ admit_age + sex_m + fx_neck_all + fx_pertroch + cci_score + lca_f + anticoag_flag + vvv2_high, data = s_full)
m_inr4   <- coxph(Surv(t90, e90) ~ admit_age + sex_m + fx_neck_all + fx_pertroch + cci_score + lca_f + inr4_f, data = s_full)

sens15 <- rbind(
  extr_hr(c_m2_full,    "p_fP3",      "1. Primary Full Cohort (P3 vs P1)", nrow(s_full)),
  extr_hr(m_chron,      "p_fP3",      "2. Chronic-Only CCI (Excl. Acute MI/Stroke)", nrow(s_full)),
  extr_hr(m_noesrd,     "p_fP3",      "3. Excluding Maintenance Dialysis / ESRD", sum(s_full$esrd_dialysis == 0)),
  extr_hr(c_m2_dev,     "p_fP3",      "4. Pre-2017 Derivation Cohort", nrow(s_dev)),
  extr_hr(c_m2_val,     "p_fP3",      "5. Post-2016 Validation (Frozen Pre-2017 Z-Scores)", nrow(s_val)),
  extr_hr(c_m2_val_era, "p_erastdP3", "6. Post-2016 Validation (Era-Standardized Z-Scores)", nrow(s_val)),
  extr_hr(m_precov,     "p_fP3",      "7. Validation Strict Pre-Pandemic (2017-2019)", sum(s_val$val_pandemic_split == "val_pre_covid_2017_2019")),
  extr_hr(m_covbnd,     "p_fP3",      "8. Validation 2019-2020 Boundary Interval", sum(s_val$val_pandemic_split == "val_boundary_2019_2020")),
  extr_hr(m_covid,      "p_fP3",      "9. Validation Strict COVID-19 Era (2020-2022)", sum(s_val$val_pandemic_split == "val_covid_2020_2022")),
  extr_hr(m_soft,       "p_fP3",      "10. Posterior Soft-Probability Weighted", nrow(s_full)),
  extr_hr(m_op,         "p_fP3",      "11. Surgically Treated Subgroup", sum(s_full$had_surgery == 1)),
  extr_hr(m_vvi2,       "vvi2_high",  "12. Diagonal VVI K=2 High-Risk Profile", nrow(s_full)),
  extr_hr(m_vvv2,       "vvv2_high",  "13. Full-Covariance VVV K=2 High-Risk Profile", nrow(s_full)),
  extr_hr(m_inr4,       "inr4_fK4_P3","14. INR-Inclusive 10-Biomarker VVI K=4 Cardiorenal", nrow(s_full)),
  extr_hr(c_m2_nh,      "p_fP3",      "15. Same-Hospital Non-Hip Major Fractures", nrow(nh))
)
write.csv(sens15, file.path(OUT_T, "SuppTable_S4a_15Arm_Sensitivity.csv"), row.names = FALSE)
cat("\n=== 15-Arm Mortality Sensitivity Table (Common-Outcome HR->RR E-Values!) ===\n")
print(sens15)

# ------------------------------------------------------------------------------
# 2. De-Circularized Non-Dialysis KDIGO AKI Sensitivity Suite (Supp Table S4b)
#    Excludes all 38 prevalent arrival-AKI cases in the prior-Cr cohort (n = 1,029)!
# ------------------------------------------------------------------------------
aki_df <- df[df$esrd_dialysis == 0 & !is.na(df$aki_7d_non_esrd), ]

run_fg <- function(dat, ev_col, t_col, lbl) {
  st <- ifelse(dat[[ev_col]] == 1, 1L,
        ifelse(dat$death_30d == 1 & dat$surv_days_365 <= 7.0, 2L, 0L))
  tm <- ifelse(st == 1 & !is.na(dat[[t_col]]), pmax(dat[[t_col]], 0.05),
        ifelse(st == 2, pmax(dat$surv_days_365, 0.05), pmin(dat$los_days, 7.0)))
  tm <- pmax(tm, 0.05)
  cov_m <- model.matrix(~ admit_age + sex_m + fx_neck_all + fx_pertroch + cci_score + lca_f + anticoag_flag + p_f, data = dat)[, -1]
  fg <- crr(tm, st, cov_m, failcode = 1, cencode = 0)
  sm <- summary(fg); ci <- sm$conf.int; pv <- sm$coef[, "p-value"]
  by_p <- dat %>% group_by(p_f) %>% summarise(inc = sprintf("%d/%d (%.1f%%)", sum(.data[[ev_col]]==1), n(), 100*mean(.data[[ev_col]]==1)), .groups="drop")
  data.frame(
    definition = lbl, n_cohort = nrow(dat),
    events_total = sprintf("%d (%.1f%%)", sum(dat[[ev_col]] == 1), 100*mean(dat[[ev_col]] == 1)),
    P1_inc = by_p$inc[by_p$p_f == "P1"], P2_inc = by_p$inc[by_p$p_f == "P2"], P3_inc = by_p$inc[by_p$p_f == "P3"],
    P2_vs_P1_sHR = sprintf("%.2f (%.2f-%.2f)", ci["p_fP2", 1], ci["p_fP2", 3], ci["p_fP2", 4]), P2_p = signif(pv["p_fP2"], 3),
    P3_vs_P1_sHR = sprintf("%.2f (%.2f-%.2f)", ci["p_fP3", 1], ci["p_fP3", 3], ci["p_fP3", 4]), P3_p = signif(pv["p_fP3"], 3),
    p2_est = round(ci["p_fP2", 1], 2), p2_lo = round(ci["p_fP2", 3], 2), p2_hi = round(ci["p_fP2", 4], 2),
    p3_est = round(ci["p_fP3", 1], 2), p3_lo = round(ci["p_fP3", 3], 2), p3_hi = round(ci["p_fP3", 4], 2),
    stringsAsFactors = FALSE
  )
}

fg_prim <- run_fg(aki_df, "aki_7d_non_esrd", "time_to_aki_non_esrd",
                  "1. Primary 7-d KDIGO Stage >=1 (+0.3 in 48h or >=1.5x ED Cr)")
fg_rel  <- run_fg(aki_df, "aki_7d_rel_nonesrd", "time_to_aki_rel_nonesrd",
                  "2. Relative-Only 7-d KDIGO (>=1.5x ED Cr; omits +0.3 criterion)")

aki_inc_prior <- aki_df[!is.na(aki_df$aki_7d_inc_rel_prior_nonesrd), ]
fg_inc_rel_pr <- run_fg(aki_inc_prior, "aki_7d_inc_rel_prior_nonesrd", "time_to_inc_rel_prior_nonesrd",
                        "3. Incident Post-Admission Relative AKI (>=1.5x prior Cr; arrival AKI excl.)")
fg_inc_kd_pr  <- run_fg(aki_inc_prior, "aki_7d_inc_kdigo_prior_nonesrd", "time_to_inc_kdigo_prior_nonesrd",
                        "4. Incident Post-Admission KDIGO AKI (+0.3 or >=1.5x prior Cr; arrival AKI excl.)")

aki_sens_tab <- rbind(fg_prim, fg_rel, fg_inc_rel_pr, fg_inc_kd_pr)
write.csv(aki_sens_tab, file.path(OUT_T, "SuppTable_S4b_KDIGO_AKI_Sensitivity.csv"), row.names = FALSE)
cat("\n=== De-Circularized Non-Dialysis AKI Sensitivity Table ===\n")
print(aki_sens_tab[, 1:10])

# ------------------------------------------------------------------------------
# 3. Clean v5 Exploratory CCW (overnight ED-to-OR Day-0 arrival)
#    Day-0 All-Eligible n = 2,485 | Operated Subgroup n = 2,216
# ------------------------------------------------------------------------------
for (v in c("hb","log_cr","log_wbc","bicarbonate")) {
  df[[paste0(v, "_ps")]] <- ifelse(is.na(df[[v]]), median(df[[v]], na.rm = TRUE), df[[v]])
}

df_ccw_all <- df[df$surg_date_anomaly_excl == 0, ]
df_ccw_all$early_surg <- ifelse(df_ccw_all$had_surgery == 1 & !is.na(df_ccw_all$days_to_surgery_clean) & df_ccw_all$days_to_surgery_clean <= 1, 1L, 0L)
df_ccw_all$t_surg_eff <- ifelse(df_ccw_all$early_surg == 1 & df_ccw_all$days_to_surgery_clean == 0, 0.50,
                         ifelse(df_ccw_all$early_surg == 1 & df_ccw_all$days_to_surgery_clean == 1, 1.25, 99.0))
df_ccw_op <- df_ccw_all[df_ccw_all$had_surgery == 1, ]

cat(sprintf("\nClean v5 Target Trial Ns: All-Eligible=%d (Early=%d, Delayed/Def=%d) | Operated=%d (Early=%d, Delayed=%d)\n",
            nrow(df_ccw_all), sum(df_ccw_all$early_surg==1), sum(df_ccw_all$early_surg==0),
            nrow(df_ccw_op), sum(df_ccw_op$early_surg==1), sum(df_ccw_op$early_surg==0)))

cnt_all <- df_ccw_all %>%
  group_by(p_f, strategy = ifelse(early_surg == 1, "Early Surgery (Day 0-1)", "Delayed or Deferred (>Day 1; incl. Non-Op)")) %>%
  summarise(n = n(), n_nonop = sum(had_surgery == 0), deaths_30d = sum(e30), mort_30d_pct = round(100*mean(e30), 1),
            deaths_90d = sum(e90), mort_90d_pct = round(100*mean(e90), 1), .groups = "drop") %>%
  mutate(population = "Day0_All_Eligible_n2485")
cnt_op <- df_ccw_op %>%
  group_by(p_f, strategy = ifelse(early_surg == 1, "Early Surgery (Day 0-1)", "Delayed Surgery (>=Day 2)")) %>%
  summarise(n = n(), n_nonop = 0L, deaths_30d = sum(e30), mort_30d_pct = round(100*mean(e30), 1),
            deaths_90d = sum(e90), mort_90d_pct = round(100*mean(e90), 1), .groups = "drop") %>%
  mutate(population = "Operated_Subgroup_n2216")
write.csv(rbind(cnt_all, cnt_op), file.path(OUT_T, "SuppTable_S7_CCW_Counts.csv"), row.names = FALSE)
print(as.data.frame(rbind(cnt_all, cnt_op)))

fit_strat_ccw <- function(dat, G = 1.5, return_details = FALSE) {
  dat$sw1 <- NA_real_; dat$sw0 <- NA_real_; dat$ps <- NA_real_
  for (k in paste0("P", 1:3)) {
    ik <- which(dat$p_f == k); dk <- dat[ik, ]
    ps_m <- glm(early_surg ~ admit_age + sex_m + fx_neck_all + fx_pertroch +
                  cci_score + lca_f + anticoag_flag + med_srv +
                  hb_ps + log_cr_ps + log_wbc_ps + bicarbonate_ps,
                data = dk, family = binomial())
    psk <- pmin(pmax(as.numeric(predict(ps_m, newdata = dk, type = "response")), 0.03), 0.97)
    pm  <- mean(dk$early_surg)
    s1  <- pmin(pmax(pm / psk, quantile(pm / psk, 0.01)), quantile(pm / psk, 0.99))
    s0  <- pmin(pmax((1 - pm) / (1 - psk), quantile((1 - pm) / (1 - psk), 0.01)), quantile((1 - pm) / (1 - psk), 0.99))
    dat$ps[ik] <- psk; dat$sw1[ik] <- s1; dat$sw0[ik] <- s0
  }

  a1_1 <- data.frame(id = dat$subject_id, p_f = dat$p_f, arm_del = 0L,
                     tstart = 0, tstop = pmin(dat$t90, G), event = ifelse(dat$t90 <= G, dat$e90, 0L), w = 1.0)
  k1 <- dat$t90 > G & dat$early_surg == 1
  a1_2 <- data.frame(id = dat$subject_id[k1], p_f = dat$p_f[k1], arm_del = 0L,
                     tstart = G, tstop = dat$t90[k1], event = dat$e90[k1], w = dat$sw1[k1])

  ts0 <- ifelse(dat$early_surg == 1, pmin(dat$t90, dat$t_surg_eff), pmin(dat$t90, G))
  ev0 <- ifelse(dat$early_surg == 1, ifelse(dat$t90 < dat$t_surg_eff, dat$e90, 0L), ifelse(dat$t90 <= G, dat$e90, 0L))
  a0_1 <- data.frame(id = dat$subject_id, p_f = dat$p_f, arm_del = 1L,
                     tstart = 0, tstop = pmax(ts0, 0.01), event = ev0, w = 1.0)
  k0 <- dat$t90 > G & dat$early_surg == 0
  a0_2 <- data.frame(id = dat$subject_id[k0], p_f = dat$p_f[k0], arm_del = 1L,
                     tstart = G, tstop = dat$t90[k0], event = dat$e90[k0], w = dat$sw0[k0])

  cp <- rbind(a1_1, a1_2, a0_1, a0_2); cp <- cp[cp$tstop > cp$tstart, ]

  est_stratum <- function(sub_cp) {
    sf0 <- survfit(Surv(tstart, tstop, event) ~ 1, data = sub_cp[sub_cp$arm_del == 0, ], weights = w)
    sf1 <- survfit(Surv(tstart, tstop, event) ~ 1, data = sub_cp[sub_cp$arm_del == 1, ], weights = w)
    r_at <- function(sf, t_eval) {
      idx <- which(sf$time <= t_eval); if (length(idx) == 0) return(0.0); 1 - sf$surv[max(idx)]
    }
    rmst_cp <- function(sf, tau) {
      tt <- c(0, sf$time[sf$time <= tau], tau); ss <- c(1, sf$surv[sf$time <= tau]); ss <- c(ss, tail(ss, 1))
      sum(diff(tt) * head(ss, -1))
    }
    c(rd30 = (r_at(sf1, 30.0) - r_at(sf0, 30.0)) * 100,
      rd90 = (r_at(sf1, 90.0) - r_at(sf0, 90.0)) * 100,
      rmst90 = rmst_cp(sf1, 90.0) - rmst_cp(sf0, 90.0),
      lhr90 = unname(coef(coxph(Surv(tstart, tstop, event) ~ arm_del, data = sub_cp, weights = w))[1]))
  }

  ov <- est_stratum(cp); p1 <- est_stratum(cp[cp$p_f == "P1", ])
  p2 <- est_stratum(cp[cp$p_f == "P2", ]); p3 <- est_stratum(cp[cp$p_f == "P3", ])
  vec <- c(ov_rd30 = ov[1], ov_rd90 = ov[2], ov_rmst90 = ov[3], ov_lhr = ov[4],
           p1_rd30 = p1[1], p1_rd90 = p1[2], p1_rmst90 = p1[3], p1_lhr = p1[4],
           p2_rd30 = p2[1], p2_rd90 = p2[2], p2_rmst90 = p2[3], p2_lhr = p2[4],
           p3_rd30 = p3[1], p3_rd90 = p3[2], p3_rmst90 = p3[3], p3_lhr = p3[4],
           int_rd30_p3_p1 = unname(p3[1] - p1[1]), int_rd90_p3_p1 = unname(p3[2] - p1[2]), int_lhr_p3_p1 = unname(p3[4] - p1[4]))
  if (return_details) return(list(est = vec, cp = cp, dat = dat))
  vec
}

ccw_all_det <- fit_strat_ccw(df_ccw_all, return_details = TRUE)
d_grace <- ccw_all_det$dat[ccw_all_det$dat$t90 > 1.5, ]
d_grace$w_ccw <- ifelse(d_grace$early_surg == 1, d_grace$sw1, d_grace$sw0)

bal_vars <- c("admit_age","sex_m","fx_neck_all","fx_pertroch","cci_score","anticoag_flag","med_srv","hb_ps","log_cr_ps","log_wbc_ps","bicarbonate_ps")
bal_lbls <- c("Age","Male sex","Femoral neck fracture","Pertrochanteric fracture","Charlson Index","Pre-admission anticoagulation","Medicine admitting service","Hemoglobin","log(Creatinine)","log(WBC)","Bicarbonate")
smd_fn <- function(x, a, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  w1 <- w[a == 1]; x1 <- x[a == 1]; w0 <- w[a == 0]; x0 <- x[a == 0]
  m1 <- sum(w1 * x1) / sum(w1); m0 <- sum(w0 * x0) / sum(w0)
  v1 <- sum(w1 * (x1 - m1)^2) / sum(w1); v0 <- sum(w0 * (x0 - m0)^2) / sum(w0)
  denom <- sqrt((v1 + v0) / 2); if (denom < 1e-8) return(0.0); abs(m1 - m0) / denom
}
bal_rows <- list()
for (grp in c("P1", "P2", "P3")) {
  sg <- d_grace[d_grace$p_f == grp, ]
  for (j in seq_along(bal_vars)) {
    v <- bal_vars[j]
    bal_rows[[length(bal_rows) + 1]] <- data.frame(
      stratum = grp, covariate = bal_lbls[j],
      smd_unweighted = round(smd_fn(sg[[v]], sg$early_surg), 3),
      smd_stabilized_ipcw = round(smd_fn(sg[[v]], sg$early_surg, sg$w_ccw), 3), stringsAsFactors = FALSE)
  }
}
bal_tab <- do.call(rbind, bal_rows)
write.csv(bal_tab, file.path(OUT_T, "SuppTable_S6_CCW_IPCW_Balance_PerPhenotype.csv"), row.names = FALSE)
cat("Max Stabilized IPCW SMD: P1=", max(bal_tab$smd_stabilized_ipcw[bal_tab$stratum=="P1"]),
    "| P2=", max(bal_tab$smd_stabilized_ipcw[bal_tab$stratum=="P2"]),
    "| P3=", max(bal_tab$smd_stabilized_ipcw[bal_tab$stratum=="P3"]), "\n")

run_ccw_boot <- function(dat, B = 1000L, cohort_lbl = "Day0_All_Eligible") {
  pt <- fit_strat_ccw(dat); n_d <- nrow(dat)
  bm <- matrix(NA_real_, B, length(pt)); colnames(bm) <- names(pt)
  for (b in seq_len(B)) {
    ib <- sample.int(n_d, n_d, replace = TRUE); db <- dat[ib, ]; db$subject_id <- seq_len(n_d)
    bm[b, ] <- tryCatch(fit_strat_ccw(db), error = function(e) rep(NA_real_, length(pt)))
  }
  bm <- bm[complete.cases(bm), , drop = FALSE]
  emp_p <- function(obs, boots, null_val = 0) {
    p_lo <- (1 + sum(boots <= null_val)) / (length(boots) + 1)
    p_hi <- (1 + sum(boots >= null_val)) / (length(boots) + 1)
    min(1.0, 2 * min(p_lo, p_hi))
  }
  D_lhr <- cbind(bm[, "p2_lhr.lhr90"] - bm[, "p1_lhr.lhr90"], bm[, "p3_lhr.lhr90"] - bm[, "p1_lhr.lhr90"])
  d_lhr <- c(pt["p2_lhr.lhr90"] - pt["p1_lhr.lhr90"], pt["p3_lhr.lhr90"] - pt["p1_lhr.lhr90"])
  p_int_global_lhr <- pchisq(as.numeric(t(d_lhr) %*% solve(cov(D_lhr)) %*% d_lhr), df = 2, lower.tail = FALSE)
  p_int_p3_p1_lhr  <- emp_p(pt["int_lhr_p3_p1"], bm[, "int_lhr_p3_p1"], 0)
  p_int_p3_p1_rd30 <- emp_p(pt["int_rd30_p3_p1"], bm[, "int_rd30_p3_p1"], 0)

  prefs <- c("ov", "p1", "p2", "p3")
  s_nms <- c(paste0(cohort_lbl, ": Overall"), "P1: Homeostatic", "P2: Stress-Inflammatory", "P3: Cardiorenal-Anemic")
  do.call(rbind, lapply(1:4, function(i) {
    pr <- prefs[i]
    c_rd30 <- paste0(pr, "_rd30.rd30"); c_rd90 <- paste0(pr, "_rd90.rd90")
    c_rmst <- paste0(pr, "_rmst90.rmst90"); c_lhr <- paste0(pr, "_lhr.lhr90")
    data.frame(
      cohort = cohort_lbl, stratum = s_nms[i],
      rd_30d_pct = round(unname(pt[c_rd30]), 2), rd_30d_lo = round(unname(quantile(bm[, c_rd30], 0.025)), 2), rd_30d_hi = round(unname(quantile(bm[, c_rd30], 0.975)), 2),
      rd_90d_pct = round(unname(pt[c_rd90]), 2), rd_90d_lo = round(unname(quantile(bm[, c_rd90], 0.025)), 2), rd_90d_hi = round(unname(quantile(bm[, c_rd90], 0.975)), 2),
      rmst_90d_diff = round(unname(pt[c_rmst]), 2), rmst_90d_lo = round(unname(quantile(bm[, c_rmst], 0.025)), 2), rmst_90d_hi = round(unname(quantile(bm[, c_rmst], 0.975)), 2),
      hr_90d = round(exp(unname(pt[c_lhr])), 2), hr_lo = round(exp(unname(quantile(bm[, c_lhr], 0.025))), 2), hr_hi = round(exp(unname(quantile(bm[, c_lhr], 0.975))), 2),
      p_boot_hr = round(emp_p(pt[c_lhr], bm[, c_lhr], 0), 4),
      p_int_global_hr = round(p_int_global_lhr, 4), p_int_p3_vs_p1_hr = round(p_int_p3_p1_lhr, 4), p_int_p3_vs_p1_rd30 = round(p_int_p3_p1_rd30, 4),
      stringsAsFactors = FALSE)
  }))
}

set.seed(20261124)
ccw_all_tab <- run_ccw_boot(df_ccw_all, B = 1000L, cohort_lbl = "Day0_All_Eligible")
set.seed(20261125)
ccw_op_tab  <- run_ccw_boot(df_ccw_op,  B = 1000L, cohort_lbl = "Operated_Subgroup")
ccw_combined <- rbind(ccw_all_tab, ccw_op_tab)
write.csv(ccw_combined, file.path(OUT_T, "Table3_v5_Exploratory_CCW.csv"), row.names = FALSE)
cat("\n=== Clean v5 Exploratory CCW Results (n=2,485 All-Elig & n=2,216 Operated) ===\n")
print(ccw_combined)

save_ckm <- function(cp_obj, st0_lbl, st1_lbl, out_fn) {
  cp_obj$strategy <- ifelse(cp_obj$arm_del == 0, st0_lbl, st1_lbl)
  rows <- list()
  for (k in paste0("P", 1:3)) {
    for (st in c(st0_lbl, st1_lbl)) {
      sub <- cp_obj[cp_obj$p_f == k & cp_obj$strategy == st, ]
      sf <- survfit(Surv(tstart, tstop, event) ~ 1, data = sub, weights = w)
      rows[[length(rows) + 1]] <- data.frame(
        p_f = k, strategy = st, time = c(0, sf$time), surv = c(1, sf$surv) * 100,
        lo = c(1, sf$lower) * 100, hi = c(1, sf$upper) * 100, stringsAsFactors = FALSE)
    }
  }
  write.csv(do.call(rbind, rows), file.path(OUT_T, out_fn), row.names = FALSE)
}
save_ckm(ccw_all_det$cp, "Early Surgery (Day 0-1)", "Delayed / Deferred (>Day 1)", "v5_ccw_km_curves.csv")
ccw_op_det <- fit_strat_ccw(df_ccw_op, return_details = TRUE)
save_ckm(ccw_op_det$cp, "Early Surgery (Day 0-1)", "Delayed Surgery (>=Day 2)", "v5_ccw_km_curves_operated.csv")
cat("=== PAPER 1 v5 PIPELINE COMPLETE ===\n")
