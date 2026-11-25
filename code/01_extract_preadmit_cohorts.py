#!/usr/bin/env python3
"""
01_extract_preadmit_cohorts.py
==============================
Standalone DuckDB/SQL cohort & pre-admission biomarker extraction for:
  "Emergency Pre-Admission Host Physiological Subphenotypes and Prognostic
   Stratification in Older Adults with Acute Hip Fracture"

Requirements:
  - Credentialed access to MIMIC-IV v3.1 loaded into a DuckDB database file
    (tables under schemas `mimiciv_hosp` and `mimiciv_icu`).
  - Python >= 3.10 with `duckdb` installed (`pip install duckdb`).

Usage:
  export MIMIC_DUCKDB_PATH="/path/to/mimiciv.duckdb"
  python3 code/01_extract_preadmit_cohorts.py
"""

import argparse
import os
from pathlib import Path
import duckdb


def main():
    parser = argparse.ArgumentParser(description="Extract MIMIC-IV v3.1 pre-admission hip & non-hip fracture cohorts.")
    parser.add_argument("--db", default=os.environ.get("MIMIC_DUCKDB_PATH", "mimiciv.duckdb"),
                        help="Path to MIMIC-IV v3.1 DuckDB database file.")
    parser.add_argument("--outdir", default=str(Path(__file__).resolve().parent.parent / "tables"),
                        help="Output directory for extracted tables (excluded from Git via .gitignore).")
    args = parser.parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(args.db, read_only=True)

    # --------------------------------------------------------------------------
    # 1. Index Acute Hip Fracture Admissions (Age >= 65 y, n = 2,825 unique patients)
    #    with [admittime - 24h, admittime] first-draw ED biomarkers, exact calendar
    #    intervals [admit_yr_min, admit_yr_max], arrival-AKI-excluded prior-Cr AKI,
    #    and rule-based surgical timing audit.
    # --------------------------------------------------------------------------
    hip_out = out_dir / "hip65_v5_preadmit_master.csv"
    con.execute(f"""
    COPY (
      WITH hip_dx AS (
        SELECT hadm_id,
               MAX(icd_version) AS icd_version,
               MAX(CASE WHEN (icd_version=9 AND (icd_code LIKE '8200%' OR icd_code LIKE '8201%' OR icd_code IN ('8208','8209')))
                          OR (icd_version=10 AND icd_code LIKE 'S720%') THEN 1 ELSE 0 END) AS fx_neck_all,
               MAX(CASE WHEN (icd_version=9 AND (icd_code LIKE '82020%' OR icd_code LIKE '82021%' OR icd_code LIKE '82030%' OR icd_code LIKE '82031%'))
                          OR (icd_version=10 AND icd_code LIKE 'S721%') THEN 1 ELSE 0 END) AS fx_pertroch,
               MAX(CASE WHEN (icd_version=9 AND (icd_code LIKE '82022%' OR icd_code LIKE '82032%'))
                          OR (icd_version=10 AND icd_code LIKE 'S722%') THEN 1 ELSE 0 END) AS fx_subtroch
        FROM mimiciv_hosp.diagnoses_icd
        WHERE (icd_version = 9 AND icd_code LIKE '820%')
           OR (icd_version = 10
               AND (icd_code LIKE 'S720%' OR icd_code LIKE 'S721%' OR icd_code LIKE 'S722%')
               AND RIGHT(icd_code, 1) IN ('A', 'B', 'C'))
        GROUP BY hadm_id
      ),
      ranked_adm AS (
        SELECT a.subject_id, a.hadm_id, a.admittime, a.dischtime, a.deathtime, a.edregtime,
               COALESCE(a.discharge_location, 'UNKNOWN') AS discharge_location,
               a.hospital_expire_flag, p.gender, p.dod, p.anchor_year, p.anchor_year_group,
               p.anchor_age + (EXTRACT(YEAR FROM a.admittime) - p.anchor_year) AS admit_age,
               CAST(SPLIT_PART(p.anchor_year_group, ' - ', 1) AS INT) + (EXTRACT(YEAR FROM a.admittime) - p.anchor_year) AS admit_yr_min,
               CAST(SPLIT_PART(p.anchor_year_group, ' - ', 2) AS INT) + (EXTRACT(YEAR FROM a.admittime) - p.anchor_year) AS admit_yr_max,
               ROUND(EXTRACT(EPOCH FROM (a.dischtime - a.admittime)) / 86400.0, 3) AS los_days,
               h.icd_version, h.fx_neck_all, h.fx_pertroch, h.fx_subtroch,
               ROW_NUMBER() OVER (PARTITION BY a.subject_id ORDER BY a.admittime ASC, a.hadm_id ASC) AS rn
        FROM mimiciv_hosp.admissions a
        JOIN mimiciv_hosp.patients p ON a.subject_id = p.subject_id
        JOIN hip_dx h ON a.hadm_id = h.hadm_id
        WHERE a.dischtime >= a.admittime
          AND (p.anchor_age + (EXTRACT(YEAR FROM a.admittime) - p.anchor_year)) >= 65
      ),
      base AS (
        SELECT *,
          CASE WHEN admit_yr_max <= 2016 THEN 'strict_pre2017'
               WHEN admit_yr_min >= 2017 THEN 'strict_post2016'
               ELSE 'boundary_2015_2018' END AS strict_temporal_split,
          CASE WHEN admit_yr_max <= 2010 THEN '2008-2010'
               WHEN admit_yr_min >= 2011 AND admit_yr_max <= 2013 THEN '2011-2013'
               WHEN admit_yr_min >= 2014 AND admit_yr_max <= 2016 THEN '2014-2016'
               WHEN admit_yr_min >= 2017 AND admit_yr_max <= 2019 THEN '2017-2019'
               WHEN admit_yr_min >= 2020 THEN '2020-2022'
               ELSE 'Boundary' END AS strict_era_bin
        FROM ranked_adm WHERE rn = 1
      ),
      surg AS (
        SELECT p.hadm_id, MIN(p.chartdate) AS first_surg_date
        FROM mimiciv_hosp.procedures_icd p
        JOIN base b ON p.hadm_id = b.hadm_id
        WHERE (p.icd_version = 9 AND (
                p.icd_code LIKE '7915%' OR p.icd_code LIKE '7925%' OR p.icd_code LIKE '7935%' OR
                p.icd_code LIKE '7855%' OR p.icd_code LIKE '8151%' OR p.icd_code LIKE '8152%' OR p.icd_code LIKE '8153%'))
           OR (p.icd_version = 10 AND (
                p.icd_code LIKE '0QS6%' OR p.icd_code LIKE '0QS7%' OR
                p.icd_code LIKE '0QH6%' OR p.icd_code LIKE '0QH7%' OR
                p.icd_code LIKE '0SR9%' OR p.icd_code LIKE '0SRA%' OR p.icd_code LIKE '0SRB%' OR
                p.icd_code LIKE '0SRE%' OR p.icd_code LIKE '0SRR%' OR p.icd_code LIKE '0SRS%'))
        GROUP BY p.hadm_id
      ),
      dx_flags AS (
        SELECT b.hadm_id,
          MAX(CASE WHEN (d.icd_version=9 AND (d.icd_code LIKE '410%' OR d.icd_code LIKE '412%')) OR (d.icd_version=10 AND (d.icd_code LIKE 'I21%' OR d.icd_code LIKE 'I22%' OR d.icd_code='I252')) THEN 1 ELSE 0 END) AS hx_mi,
          MAX(CASE WHEN (d.icd_version=9 AND d.icd_code LIKE '428%') OR (d.icd_version=10 AND d.icd_code LIKE 'I50%') THEN 1 ELSE 0 END) AS hx_chf,
          MAX(CASE WHEN (d.icd_version=9 AND (d.icd_code LIKE '440%' OR d.icd_code LIKE '441%' OR d.icd_code LIKE '4439%')) OR (d.icd_version=10 AND (d.icd_code LIKE 'I70%' OR d.icd_code LIKE 'I71%' OR d.icd_code LIKE 'I739%')) THEN 1 ELSE 0 END) AS hx_pvd,
          MAX(CASE WHEN (d.icd_version=9 AND SUBSTR(d.icd_code,1,3) BETWEEN '430' AND '438') OR (d.icd_version=10 AND SUBSTR(d.icd_code,1,3) BETWEEN 'I60' AND 'I69') THEN 1 ELSE 0 END) AS hx_cevd,
          MAX(CASE WHEN (d.icd_version=9 AND (d.icd_code LIKE '290%' OR d.icd_code LIKE '3310%')) OR (d.icd_version=10 AND (d.icd_code LIKE 'F01%' OR d.icd_code LIKE 'F02%' OR d.icd_code LIKE 'F03%' OR d.icd_code LIKE 'G30%')) THEN 1 ELSE 0 END) AS hx_dementia,
          MAX(CASE WHEN (d.icd_version=9 AND SUBSTR(d.icd_code,1,3) BETWEEN '490' AND '496') OR (d.icd_version=10 AND SUBSTR(d.icd_code,1,3) BETWEEN 'J40' AND 'J47') THEN 1 ELSE 0 END) AS hx_copd,
          MAX(CASE WHEN (d.icd_version=9 AND d.icd_code LIKE '2500%') OR (d.icd_version=10 AND SUBSTR(d.icd_code,1,4) IN ('E109','E119','E139','E149')) THEN 1 ELSE 0 END) AS hx_dm_uncomp,
          MAX(CASE WHEN (d.icd_version=9 AND d.icd_code LIKE '250%' AND d.icd_code NOT LIKE '2500%') OR (d.icd_version=10 AND SUBSTR(d.icd_code,1,3) IN ('E10','E11','E13','E14') AND SUBSTR(d.icd_code,1,4) NOT IN ('E109','E119','E139','E149')) THEN 1 ELSE 0 END) AS hx_dm_comp,
          MAX(CASE WHEN (d.icd_version=9 AND SUBSTR(d.icd_code,1,3) IN ('582','583','585','586')) OR (d.icd_version=10 AND SUBSTR(d.icd_code,1,3) IN ('N18','N19')) THEN 1 ELSE 0 END) AS hx_renal,
          MAX(CASE WHEN (d.icd_version=9 AND d.icd_code IN ('5856','V4511','V4512','V560')) OR (d.icd_version=10 AND d.icd_code IN ('N186','Z992','Z4931')) THEN 1 ELSE 0 END) AS esrd_dialysis,
          MAX(CASE WHEN (d.icd_version=9 AND d.icd_code = 'V5861') OR (d.icd_version=10 AND d.icd_code LIKE 'Z7901%') THEN 1 ELSE 0 END) AS dx_longterm_anticoag,
          MAX(CASE WHEN (d.icd_version=9 AND d.icd_code LIKE '410%') OR (d.icd_version=10 AND d.icd_code LIKE 'I21%') THEN 1 ELSE 0 END) AS acute_mi_flag,
          MAX(CASE WHEN (d.icd_version=9 AND (d.icd_code LIKE '433%1' OR d.icd_code LIKE '434%1')) OR (d.icd_version=10 AND d.icd_code LIKE 'I63%') THEN 1 ELSE 0 END) AS acute_stroke_flag
        FROM base b
        LEFT JOIN mimiciv_hosp.diagnoses_icd d ON b.hadm_id = d.hadm_id
        GROUP BY 1
      ),
      srv AS (
        SELECT hadm_id, curr_service AS first_service,
               ROW_NUMBER() OVER (PARTITION BY hadm_id ORDER BY transfertime ASC) AS rn
        FROM mimiciv_hosp.services
      ),
      icu_f AS (
        SELECT hadm_id, 1 AS icu_admitted FROM mimiciv_icu.icustays GROUP BY hadm_id
      ),
      prior_rx AS (
        SELECT b.hadm_id,
          MAX(CASE WHEN LOWER(p.drug) LIKE '%apixaban%' OR LOWER(p.drug) LIKE '%eliquis%'
                     OR LOWER(p.drug) LIKE '%rivaroxaban%' OR LOWER(p.drug) LIKE '%xarelto%'
                     OR LOWER(p.drug) LIKE '%dabigatran%' OR LOWER(p.drug) LIKE '%edoxaban%' THEN 1 ELSE 0 END) AS prior_hosp_doac,
          MAX(CASE WHEN LOWER(p.drug) LIKE '%warfarin%' OR LOWER(p.drug) LIKE '%coumadin%' THEN 1 ELSE 0 END) AS prior_hosp_warfarin
        FROM base b
        LEFT JOIN mimiciv_hosp.prescriptions p
          ON b.subject_id = p.subject_id AND p.hadm_id != b.hadm_id
         AND p.starttime IS NOT NULL
         AND CAST(p.starttime AS TIMESTAMP) < b.admittime
         AND CAST(p.starttime AS TIMESTAMP) >= b.admittime - INTERVAL 365 DAY
        GROUP BY 1
      ),
      ed_labs_ranked AS (
        SELECT b.hadm_id,
          CASE
            WHEN le.itemid IN (51222, 50811) AND le.valuenum BETWEEN 3 AND 22 THEN 'hb'
            WHEN le.itemid = 51277 AND le.valuenum BETWEEN 9 AND 35 THEN 'rdw'
            WHEN le.itemid = 51265 AND le.valuenum BETWEEN 5 AND 1500 THEN 'platelets'
            WHEN le.itemid IN (51301, 51300) AND le.valuenum BETWEEN 0.2 AND 150 THEN 'wbc'
            WHEN le.itemid = 50912 AND le.valuenum BETWEEN 0.1 AND 25 THEN 'creatinine'
            WHEN le.itemid IN (50931, 50809) AND le.valuenum BETWEEN 20 AND 1200 THEN 'glucose'
            WHEN le.itemid IN (50983, 50824) AND le.valuenum BETWEEN 105 AND 175 THEN 'sodium'
            WHEN le.itemid IN (50971, 50822) AND le.valuenum BETWEEN 1.5 AND 9.5 THEN 'potassium'
            WHEN le.itemid IN (50882, 50803) AND le.valuenum BETWEEN 5 AND 50 THEN 'bicarbonate'
            WHEN le.itemid = 51237 AND le.valuenum BETWEEN 0.5 AND 20 THEN 'inr'
            WHEN le.itemid = 50862 AND le.valuenum BETWEEN 1.0 AND 6.5 THEN 'ext_albumin'
            WHEN le.itemid = 50813 AND le.valuenum BETWEEN 0.2 AND 30 THEN 'ext_lactate'
            WHEN le.itemid = 51003 AND le.valuenum >= 0 THEN 'ext_troponin_t'
            WHEN le.itemid = 50963 AND le.valuenum >= 0 THEN 'ext_ntprobnp'
            ELSE NULL
          END AS biomarker,
          le.valuenum, le.charttime, le.labevent_id
        FROM base b
        JOIN mimiciv_hosp.labevents le
          ON b.subject_id = le.subject_id AND le.valuenum IS NOT NULL
         AND le.charttime >= b.admittime - INTERVAL 24 HOUR
         AND le.charttime <= b.admittime
         AND COALESCE(le.storetime, le.charttime) <= b.admittime + INTERVAL 2 HOUR
      ),
      ed_first AS (
        SELECT hadm_id, biomarker, valuenum,
               ROW_NUMBER() OVER (PARTITION BY hadm_id, biomarker ORDER BY charttime ASC, labevent_id ASC) AS rn
        FROM ed_labs_ranked WHERE biomarker IS NOT NULL
      ),
      ed_piv AS (
        SELECT hadm_id,
          MAX(CASE WHEN biomarker='hb' THEN valuenum END) AS hb,
          MAX(CASE WHEN biomarker='rdw' THEN valuenum END) AS rdw,
          MAX(CASE WHEN biomarker='platelets' THEN valuenum END) AS platelets,
          MAX(CASE WHEN biomarker='wbc' THEN valuenum END) AS wbc,
          MAX(CASE WHEN biomarker='creatinine' THEN valuenum END) AS creatinine,
          MAX(CASE WHEN biomarker='glucose' THEN valuenum END) AS glucose,
          MAX(CASE WHEN biomarker='sodium' THEN valuenum END) AS sodium,
          MAX(CASE WHEN biomarker='potassium' THEN valuenum END) AS potassium,
          MAX(CASE WHEN biomarker='bicarbonate' THEN valuenum END) AS bicarbonate,
          MAX(CASE WHEN biomarker='inr' THEN valuenum END) AS inr,
          MAX(CASE WHEN biomarker='ext_albumin' THEN valuenum END) AS ext_albumin,
          MAX(CASE WHEN biomarker='ext_lactate' THEN valuenum END) AS ext_lactate,
          MAX(CASE WHEN biomarker='ext_troponin_t' THEN valuenum END) AS ext_troponin_t,
          MAX(CASE WHEN biomarker='ext_ntprobnp' THEN valuenum END) AS ext_ntprobnp
        FROM ed_first WHERE rn = 1 GROUP BY 1
      ),
      prior_cr AS (
        SELECT b.hadm_id, MEDIAN(le.valuenum) AS prior_cr_med_365d
        FROM base b
        JOIN mimiciv_hosp.labevents le
          ON b.subject_id = le.subject_id AND le.itemid = 50912 AND le.valuenum BETWEEN 0.1 AND 25
         AND le.charttime >= b.admittime - INTERVAL 365 DAY
         AND le.charttime <= b.admittime - INTERVAL 7 DAY
        GROUP BY 1
      ),
      cr_post7 AS (
        SELECT b.hadm_id, ep.creatinine AS ed_base_cr, pc.prior_cr_med_365d, le.valuenum AS cr_val,
               EXTRACT(EPOCH FROM (le.charttime - b.admittime))/86400.0 AS days_post_adm
        FROM base b
        JOIN ed_piv ep ON b.hadm_id = ep.hadm_id
        LEFT JOIN prior_cr pc ON b.hadm_id = pc.hadm_id
        JOIN mimiciv_hosp.labevents le
          ON b.subject_id = le.subject_id AND le.itemid = 50912 AND le.valuenum BETWEEN 0.1 AND 25
         AND le.charttime > b.admittime AND le.charttime <= b.admittime + INTERVAL 7 DAY
        WHERE ep.creatinine IS NOT NULL
      ),
      aki_all AS (
        SELECT hadm_id,
          MAX(CASE WHEN (days_post_adm <= 2.0 AND cr_val - ed_base_cr >= 0.3) OR (days_post_adm <= 7.0 AND cr_val >= 1.5 * ed_base_cr) THEN 1 ELSE 0 END) AS aki_prim,
          MIN(CASE WHEN (days_post_adm <= 2.0 AND cr_val - ed_base_cr >= 0.3) OR (days_post_adm <= 7.0 AND cr_val >= 1.5 * ed_base_cr) THEN days_post_adm END) AS t_aki_prim,
          MAX(CASE WHEN days_post_adm <= 7.0 AND cr_val >= 1.5 * ed_base_cr THEN 1 ELSE 0 END) AS aki_rel,
          MIN(CASE WHEN days_post_adm <= 7.0 AND cr_val >= 1.5 * ed_base_cr THEN days_post_adm END) AS t_aki_rel,
          MAX(CASE WHEN prior_cr_med_365d IS NOT NULL AND cr_val >= 1.5 * prior_cr_med_365d THEN 1 ELSE 0 END) AS inc_rel_prior,
          MIN(CASE WHEN prior_cr_med_365d IS NOT NULL AND cr_val >= 1.5 * prior_cr_med_365d THEN days_post_adm END) AS t_inc_rel_prior,
          MAX(CASE WHEN prior_cr_med_365d IS NOT NULL AND ((days_post_adm <= 2.0 AND cr_val - ed_base_cr >= 0.3) OR cr_val >= 1.5 * prior_cr_med_365d) THEN 1 ELSE 0 END) AS inc_kd_prior,
          MIN(CASE WHEN prior_cr_med_365d IS NOT NULL AND ((days_post_adm <= 2.0 AND cr_val - ed_base_cr >= 0.3) OR cr_val >= 1.5 * prior_cr_med_365d) THEN days_post_adm END) AS t_inc_kd_prior
        FROM cr_post7 GROUP BY 1
      )
      SELECT b.subject_id, b.hadm_id, b.admittime, b.dischtime, b.admit_age, b.gender,
             b.anchor_year_group, b.admit_yr_min, b.admit_yr_max, b.strict_temporal_split, b.strict_era_bin,
             CASE WHEN b.strict_temporal_split != 'strict_post2016' THEN b.strict_temporal_split
                  WHEN b.admit_yr_max <= 2019 THEN 'val_pre_covid_2017_2019'
                  WHEN b.admit_yr_min >= 2020 THEN 'val_covid_2020_2022'
                  ELSE 'val_boundary_2019_2020' END AS val_pandemic_split,
             b.los_days, COALESCE(ic.icu_admitted, 0) AS icu_admitted,
             COALESCE(sv.first_service, 'OTHER') AS first_service, b.discharge_location,
             b.fx_neck_all, b.fx_pertroch, b.fx_subtroch,
             df.hx_mi, df.hx_chf, df.hx_pvd, df.hx_cevd, df.hx_dementia, df.hx_copd,
             df.hx_dm_uncomp, df.hx_dm_comp, df.hx_renal, df.esrd_dialysis, df.dx_longterm_anticoag,
             (df.hx_mi + df.hx_chf + df.hx_pvd + df.hx_cevd + df.hx_dementia + df.hx_copd + df.hx_dm_uncomp + 2*df.hx_dm_comp + 2*df.hx_renal) AS cci_score,
             GREATEST(0, (df.hx_mi + df.hx_chf + df.hx_pvd + df.hx_cevd + df.hx_dementia + df.hx_copd + df.hx_dm_uncomp + 2*df.hx_dm_comp + 2*df.hx_renal) - df.acute_mi_flag - df.acute_stroke_flag) AS cci_chronic_only,
             COALESCE(prx.prior_hosp_doac, 0) AS prior_hosp_doac,
             COALESCE(prx.prior_hosp_warfarin, 0) AS prior_hosp_warfarin,
             CASE WHEN s.hadm_id IS NOT NULL THEN 1 ELSE 0 END AS had_surgery,
             s.first_surg_date,
             CASE WHEN s.first_surg_date IS NOT NULL THEN DATE_DIFF('day', CAST(b.admittime AS DATE), s.first_surg_date) END AS days_to_surgery,
             -- Rule-based surgical date validation (zero hardcoded admission IDs):
             CASE WHEN s.first_surg_date IS NOT NULL AND s.first_surg_date < CAST(COALESCE(b.edregtime, b.admittime) AS DATE) THEN 1 ELSE 0 END AS surg_date_anomaly_excl,
             CASE WHEN s.first_surg_date IS NOT NULL AND s.first_surg_date = CAST(b.edregtime AS DATE) AND s.first_surg_date < CAST(b.admittime AS DATE) THEN 0
                  WHEN s.first_surg_date IS NOT NULL THEN DATE_DIFF('day', CAST(b.admittime AS DATE), s.first_surg_date) END AS days_to_surgery_clean,
             ep.hb, ep.rdw, ep.platelets, ep.wbc, ep.creatinine, ep.glucose, ep.sodium, ep.potassium, ep.bicarbonate, ep.inr,
             ep.ext_albumin, ep.ext_lactate, ep.ext_troponin_t, ep.ext_ntprobnp,
             pc.prior_cr_med_365d,
             CASE WHEN df.esrd_dialysis = 1 THEN NULL ELSE COALESCE(ak.aki_prim, 0) END AS aki_7d_non_esrd,
             CASE WHEN df.esrd_dialysis = 1 THEN NULL ELSE ak.t_aki_prim END AS time_to_aki_non_esrd,
             CASE WHEN df.esrd_dialysis = 1 THEN NULL ELSE COALESCE(ak.aki_rel, 0) END AS aki_7d_rel_nonesrd,
             CASE WHEN df.esrd_dialysis = 1 THEN NULL ELSE ak.t_aki_rel END AS time_to_aki_rel_nonesrd,
             CASE WHEN df.esrd_dialysis = 1 OR pc.prior_cr_med_365d IS NULL OR ep.creatinine >= 1.5 * pc.prior_cr_med_365d THEN NULL
                  ELSE COALESCE(ak.inc_rel_prior, 0) END AS aki_7d_inc_rel_prior_nonesrd,
             CASE WHEN df.esrd_dialysis = 1 OR pc.prior_cr_med_365d IS NULL OR ep.creatinine >= 1.5 * pc.prior_cr_med_365d THEN NULL
                  ELSE ak.t_inc_rel_prior END AS time_to_inc_rel_prior_nonesrd,
             CASE WHEN df.esrd_dialysis = 1 OR pc.prior_cr_med_365d IS NULL OR ep.creatinine >= 1.5 * pc.prior_cr_med_365d THEN NULL
                  ELSE COALESCE(ak.inc_kd_prior, 0) END AS aki_7d_inc_kdigo_prior_nonesrd,
             CASE WHEN df.esrd_dialysis = 1 OR pc.prior_cr_med_365d IS NULL OR ep.creatinine >= 1.5 * pc.prior_cr_med_365d THEN NULL
                  ELSE ak.t_inc_kd_prior END AS time_to_inc_kdigo_prior_nonesrd,
             b.hospital_expire_flag,
             CASE WHEN b.dod IS NOT NULL AND DATE_DIFF('day', CAST(b.admittime AS DATE), b.dod) <= 30  THEN 1 ELSE 0 END AS death_30d,
             CASE WHEN b.dod IS NOT NULL AND DATE_DIFF('day', CAST(b.admittime AS DATE), b.dod) <= 90  THEN 1 ELSE 0 END AS death_90d,
             CASE WHEN b.dod IS NOT NULL AND DATE_DIFF('day', CAST(b.admittime AS DATE), b.dod) <= 365 THEN 1 ELSE 0 END AS death_1y,
             CASE WHEN b.deathtime IS NOT NULL THEN LEAST(365.0, GREATEST(0.05, EXTRACT(EPOCH FROM (b.deathtime - b.admittime))/86400.0))
                  WHEN b.dod IS NOT NULL THEN LEAST(365.0, GREATEST(0.5, CAST(DATE_DIFF('day', CAST(b.admittime AS DATE), b.dod) AS DOUBLE)))
                  ELSE 365.0 END AS surv_days_365
      FROM base b
      LEFT JOIN surg s ON b.hadm_id = s.hadm_id
      LEFT JOIN dx_flags df ON b.hadm_id = df.hadm_id
      LEFT JOIN srv sv ON b.hadm_id = sv.hadm_id AND sv.rn = 1
      LEFT JOIN icu_f ic ON b.hadm_id = ic.hadm_id
      LEFT JOIN prior_rx prx ON b.hadm_id = prx.hadm_id
      LEFT JOIN ed_piv ep ON b.hadm_id = ep.hadm_id
      LEFT JOIN prior_cr pc ON b.hadm_id = pc.hadm_id
      LEFT JOIN aki_all ak ON b.hadm_id = ak.hadm_id
    ) TO '{hip_out}' (HEADER, DELIMITER ',')
    """)

    # 2. Export Institutional Laboratory Reference-Interval Shift Audit (Supp Table S2d)
    s2d_out = out_dir / "SuppTable_S2d_BIDMC_Assay_Reference_Shifts.csv"
    con.execute(f"""
    COPY (
      WITH h AS (
        SELECT subject_id, hadm_id, CAST(admittime AS TIMESTAMP) AS admittime, strict_era_bin
        FROM read_csv_auto('{hip_out}')
        WHERE strict_era_bin != 'Boundary' AND
              (CASE WHEN hb IS NOT NULL THEN 1 ELSE 0 END + CASE WHEN rdw IS NOT NULL THEN 1 ELSE 0 END +
               CASE WHEN platelets IS NOT NULL THEN 1 ELSE 0 END + CASE WHEN wbc IS NOT NULL THEN 1 ELSE 0 END +
               CASE WHEN creatinine IS NOT NULL THEN 1 ELSE 0 END + CASE WHEN glucose IS NOT NULL THEN 1 ELSE 0 END +
               CASE WHEN sodium IS NOT NULL THEN 1 ELSE 0 END + CASE WHEN potassium IS NOT NULL THEN 1 ELSE 0 END +
               CASE WHEN bicarbonate IS NOT NULL THEN 1 ELSE 0 END) >= 7
      ),
      ed_first AS (
        SELECT h.hadm_id, h.strict_era_bin, le.itemid, le.valuenum, le.ref_range_lower, le.ref_range_upper,
               ROW_NUMBER() OVER (PARTITION BY h.hadm_id, le.itemid ORDER BY le.charttime ASC, le.labevent_id ASC) AS rn
        FROM h
        JOIN mimiciv_hosp.labevents le
          ON h.subject_id = le.subject_id
         AND le.itemid IN (51222, 51277, 51265, 51301, 50912, 50931, 50983, 50971, 50882)
         AND le.valuenum IS NOT NULL
         AND le.charttime >= h.admittime - INTERVAL 24 HOUR AND le.charttime <= h.admittime
      )
      SELECT CASE itemid
               WHEN 51222 THEN 'Hemoglobin (g/dL)' WHEN 51277 THEN 'RDW (%)'
               WHEN 51265 THEN 'Platelets (x10^3/uL)' WHEN 51301 THEN 'WBC (x10^3/uL)'
               WHEN 50912 THEN 'Creatinine (mg/dL)' WHEN 50931 THEN 'Glucose (mg/dL)'
               WHEN 50983 THEN 'Sodium (mEq/L)' WHEN 50971 THEN 'Potassium (mEq/L)'
               WHEN 50882 THEN 'Bicarbonate (mEq/L)' END AS biomarker,
             strict_era_bin AS admission_period,
             COUNT(*) AS n_measured,
             CAST(MODE(ref_range_lower) AS VARCHAR) || '-' || CAST(MODE(ref_range_upper) AS VARCHAR) AS modal_ref_interval,
             ROUND(MEDIAN(valuenum), 2) AS median_val
      FROM ed_first WHERE rn = 1
      GROUP BY itemid, strict_era_bin ORDER BY itemid, strict_era_bin
    ) TO '{s2d_out}' (HEADER, DELIMITER ',')
    """)

    # 3. Longitudinal 72-h Serum Creatinine Stream (for Supplementary Figure S3a lcmm)
    cr72_out = out_dir / "hip65_longitudinal_cr_96h.csv"
    con.execute(f"""
    COPY (
      SELECT h.subject_id, h.hadm_id,
             ROUND(EXTRACT(EPOCH FROM (le.charttime - CAST(h.admittime AS TIMESTAMP))) / 3600.0, 2) AS hours_from_admit,
             le.valuenum AS creatinine
      FROM read_csv_auto('{hip_out}') h
      JOIN mimiciv_hosp.labevents le
        ON h.subject_id = le.subject_id AND le.itemid = 50912 AND le.valuenum BETWEEN 0.1 AND 25
       AND le.charttime >= CAST(h.admittime AS TIMESTAMP) - INTERVAL 12 HOUR
       AND le.charttime <= CAST(h.admittime AS TIMESTAMP) + INTERVAL 96 HOUR
    ) TO '{cr72_out}' (HEADER, DELIMITER ',')
    """)
    # 4. Same-Hospital Non-Hip Major Fracture Generalization Cohort (Age >= 65 y)
    nh_out = out_dir / "nonhip_v5_preadmit.csv"
    con.execute(f"""
    COPY (
      WITH hip_any AS (
        SELECT DISTINCT hadm_id FROM mimiciv_hosp.diagnoses_icd
        WHERE (icd_version = 9 AND icd_code LIKE '820%')
           OR (icd_version = 10 AND SUBSTR(icd_code, 1, 4) IN ('S720', 'S721', 'S722'))
      ),
      nh_dx AS (
        SELECT hadm_id,
          MAX(CASE WHEN (icd_version=9 AND SUBSTR(icd_code,1,3) IN ('805','806')) OR (icd_version=10 AND SUBSTR(icd_code,1,3) IN ('S12','S22','S32') AND SUBSTR(icd_code,1,4)!='S324') THEN 1 ELSE 0 END) AS fx_spine,
          MAX(CASE WHEN (icd_version=9 AND icd_code LIKE '808%') OR (icd_version=10 AND SUBSTR(icd_code,1,4) IN ('S323','S324','S325','S326','S328')) THEN 1 ELSE 0 END) AS fx_pelvis,
          MAX(CASE WHEN (icd_version=9 AND SUBSTR(icd_code,1,3) IN ('812','813')) OR (icd_version=10 AND SUBSTR(icd_code,1,3) IN ('S42','S52')) THEN 1 ELSE 0 END) AS fx_upper,
          MAX(CASE WHEN (icd_version=9 AND SUBSTR(icd_code,1,3) IN ('821','822','823','824')) OR (icd_version=10 AND (SUBSTR(icd_code,1,4) IN ('S723','S724','S728','S729') OR SUBSTR(icd_code,1,3)='S82')) THEN 1 ELSE 0 END) AS fx_lower
        FROM mimiciv_hosp.diagnoses_icd
        WHERE hadm_id NOT IN (SELECT hadm_id FROM hip_any)
          AND ((icd_version=9 AND (SUBSTR(icd_code,1,3) IN ('805','806','807','808','812','813','821','822','823','824')))
            OR (icd_version=10 AND RIGHT(icd_code,1) IN ('A','B','C') AND SUBSTR(icd_code,1,3) IN ('S12','S22','S32','S42','S52','S72','S82')))
        GROUP BY 1
      ),
      nh_adm AS (
        SELECT a.subject_id, a.hadm_id, a.admittime, p.gender, p.dod, a.deathtime, p.anchor_year_group,
               p.anchor_age + (EXTRACT(YEAR FROM a.admittime) - p.anchor_year) AS admit_age,
               nx.fx_spine, nx.fx_pelvis, nx.fx_upper, nx.fx_lower,
               ROW_NUMBER() OVER (PARTITION BY a.subject_id ORDER BY a.admittime ASC, a.hadm_id ASC) AS rn
        FROM mimiciv_hosp.admissions a
        JOIN mimiciv_hosp.patients p ON a.subject_id = p.subject_id
        JOIN nh_dx nx ON a.hadm_id = nx.hadm_id
        WHERE a.dischtime >= a.admittime
          AND (p.anchor_age + (EXTRACT(YEAR FROM a.admittime) - p.anchor_year)) >= 65
      ),
      nh AS (SELECT * FROM nh_adm WHERE rn = 1),
      cm AS (
        SELECT n.hadm_id,
          MAX(CASE WHEN (d.icd_version=9 AND (d.icd_code LIKE '410%' OR d.icd_code LIKE '412%')) OR (d.icd_version=10 AND (d.icd_code LIKE 'I21%' OR d.icd_code LIKE 'I22%' OR d.icd_code='I252')) THEN 1 ELSE 0 END) AS hx_mi,
          MAX(CASE WHEN (d.icd_version=9 AND d.icd_code LIKE '428%') OR (d.icd_version=10 AND d.icd_code LIKE 'I50%') THEN 1 ELSE 0 END) AS hx_chf,
          MAX(CASE WHEN (d.icd_version=9 AND (d.icd_code LIKE '440%' OR d.icd_code LIKE '441%' OR d.icd_code LIKE '4439%')) OR (d.icd_version=10 AND (d.icd_code LIKE 'I70%' OR d.icd_code LIKE 'I71%' OR d.icd_code LIKE 'I739%')) THEN 1 ELSE 0 END) AS hx_pvd,
          MAX(CASE WHEN (d.icd_version=9 AND SUBSTR(d.icd_code,1,3) BETWEEN '430' AND '438') OR (d.icd_version=10 AND SUBSTR(d.icd_code,1,3) BETWEEN 'I60' AND 'I69') THEN 1 ELSE 0 END) AS hx_cevd,
          MAX(CASE WHEN (d.icd_version=9 AND (d.icd_code LIKE '290%' OR d.icd_code LIKE '3310%')) OR (d.icd_version=10 AND (d.icd_code LIKE 'F01%' OR d.icd_code LIKE 'F02%' OR d.icd_code LIKE 'F03%' OR d.icd_code LIKE 'G30%')) THEN 1 ELSE 0 END) AS hx_dementia,
          MAX(CASE WHEN (d.icd_version=9 AND SUBSTR(d.icd_code,1,3) BETWEEN '490' AND '496') OR (d.icd_version=10 AND SUBSTR(d.icd_code,1,3) BETWEEN 'J40' AND 'J47') THEN 1 ELSE 0 END) AS hx_copd,
          MAX(CASE WHEN (d.icd_version=9 AND d.icd_code LIKE '250%') OR (d.icd_version=10 AND SUBSTR(d.icd_code,1,3) IN ('E10','E11','E13','E14')) THEN 1 ELSE 0 END) AS hx_dm_any,
          MAX(CASE WHEN (d.icd_version=9 AND SUBSTR(d.icd_code,1,3) IN ('582','583','585','586')) OR (d.icd_version=10 AND SUBSTR(d.icd_code,1,3) IN ('N18','N19')) THEN 1 ELSE 0 END) AS hx_renal
        FROM nh n LEFT JOIN mimiciv_hosp.diagnoses_icd d ON n.hadm_id = d.hadm_id GROUP BY 1
      ),
      ed_r AS (
        SELECT n.hadm_id,
          CASE
            WHEN le.itemid IN (51222, 50811) AND le.valuenum BETWEEN 3 AND 22 THEN 'hb'
            WHEN le.itemid = 51277 AND le.valuenum BETWEEN 9 AND 35 THEN 'rdw'
            WHEN le.itemid = 51265 AND le.valuenum BETWEEN 5 AND 1500 THEN 'platelets'
            WHEN le.itemid IN (51301, 51300) AND le.valuenum BETWEEN 0.2 AND 150 THEN 'wbc'
            WHEN le.itemid = 50912 AND le.valuenum BETWEEN 0.1 AND 25 THEN 'creatinine'
            WHEN le.itemid IN (50931, 50809) AND le.valuenum BETWEEN 20 AND 1200 THEN 'glucose'
            WHEN le.itemid IN (50983, 50824) AND le.valuenum BETWEEN 105 AND 175 THEN 'sodium'
            WHEN le.itemid IN (50971, 50822) AND le.valuenum BETWEEN 1.5 AND 9.5 THEN 'potassium'
            WHEN le.itemid IN (50882, 50803) AND le.valuenum BETWEEN 5 AND 50 THEN 'bicarbonate'
          END AS biomarker,
          le.valuenum, le.charttime, le.labevent_id
        FROM nh n JOIN mimiciv_hosp.labevents le
          ON n.subject_id = le.subject_id AND le.valuenum IS NOT NULL
         AND le.charttime >= n.admittime - INTERVAL 24 HOUR AND le.charttime <= n.admittime
         AND COALESCE(le.storetime, le.charttime) <= n.admittime + INTERVAL 2 HOUR
      ),
      ed_f AS (
        SELECT hadm_id, biomarker, valuenum,
               ROW_NUMBER() OVER (PARTITION BY hadm_id, biomarker ORDER BY charttime ASC, labevent_id ASC) AS rn
        FROM ed_r WHERE biomarker IS NOT NULL
      ),
      ed_p AS (
        SELECT hadm_id,
          MAX(CASE WHEN biomarker='hb' THEN valuenum END) AS hb,
          MAX(CASE WHEN biomarker='rdw' THEN valuenum END) AS rdw,
          MAX(CASE WHEN biomarker='platelets' THEN valuenum END) AS platelets,
          MAX(CASE WHEN biomarker='wbc' THEN valuenum END) AS wbc,
          MAX(CASE WHEN biomarker='creatinine' THEN valuenum END) AS creatinine,
          MAX(CASE WHEN biomarker='glucose' THEN valuenum END) AS glucose,
          MAX(CASE WHEN biomarker='sodium' THEN valuenum END) AS sodium,
          MAX(CASE WHEN biomarker='potassium' THEN valuenum END) AS potassium,
          MAX(CASE WHEN biomarker='bicarbonate' THEN valuenum END) AS bicarbonate
        FROM ed_f WHERE rn = 1 GROUP BY 1
      )
      SELECT n.subject_id, n.hadm_id, n.admit_age, n.gender, n.anchor_year_group,
             n.fx_spine, n.fx_pelvis, n.fx_upper, n.fx_lower,
             cm.hx_mi, cm.hx_chf, cm.hx_pvd, cm.hx_cevd, cm.hx_dementia, cm.hx_copd, cm.hx_dm_any, cm.hx_renal,
             (cm.hx_mi + cm.hx_chf + cm.hx_pvd + cm.hx_cevd + cm.hx_dementia + cm.hx_copd + cm.hx_dm_any + 2*cm.hx_renal) AS cci_score,
             ep.hb, ep.rdw, ep.platelets, ep.wbc, ep.creatinine, ep.glucose, ep.sodium, ep.potassium, ep.bicarbonate,
             CASE WHEN n.dod IS NOT NULL AND DATE_DIFF('day', CAST(n.admittime AS DATE), n.dod) <= 30  THEN 1 ELSE 0 END AS death_30d,
             CASE WHEN n.dod IS NOT NULL AND DATE_DIFF('day', CAST(n.admittime AS DATE), n.dod) <= 90  THEN 1 ELSE 0 END AS death_90d,
             CASE WHEN n.dod IS NOT NULL AND DATE_DIFF('day', CAST(n.admittime AS DATE), n.dod) <= 365 THEN 1 ELSE 0 END AS death_1y,
             CASE WHEN n.deathtime IS NOT NULL THEN LEAST(365.0, GREATEST(0.05, EXTRACT(EPOCH FROM (n.deathtime - n.admittime))/86400.0))
                  WHEN n.dod IS NOT NULL THEN LEAST(365.0, GREATEST(0.5, CAST(DATE_DIFF('day', CAST(n.admittime AS DATE), n.dod) AS DOUBLE)))
                  ELSE 365.0 END AS surv_days_365
      FROM nh n LEFT JOIN cm ON n.hadm_id = cm.hadm_id LEFT JOIN ed_p ep ON n.hadm_id = ep.hadm_id
    ) TO '{nh_out}' (HEADER, DELIMITER ',')
    """)
    print(f"Extraction complete -> {out_dir}")


if __name__ == "__main__":
    main()
