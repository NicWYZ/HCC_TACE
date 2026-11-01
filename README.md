# Sequencing and Timing of TACE and Systemic Therapy in Advanced HCC

Analysis code for **"Optimizing the Sequencing and Timing of Transarterial Chemoembolization and Systemic Therapy in Advanced Hepatocellular Carcinoma"**, a retrospective cohort study using the National Cancer Database (NCDB).

- Manuscript submitted to **eBioMedicine**.
- Selected for a **poster presentation at ASCO GI 2026** (ASCO Gastrointestinal Cancers Symposium).

## Summary

**Question.** TACE and systemic therapy work well together in advanced hepatocellular carcinoma (HCC). It is not known which should come first, or how far apart they should be given.

**Cohort.** 1,403 adults with AJCC stage III–IV HCC (NCDB 2004–2022) who received both TACE and systemic therapy. Patients were grouped by treatment order:

| Group | Definition | n (%) |
|---|---|---|
| TACE-first | TACE, then systemic therapy >3 days later | 916 (65.3%) |
| Systemic-first | Systemic therapy, then TACE >3 days later | 202 (14.4%) |
| Simultaneous | Both within 3 days | 285 (20.3%) |

**Methods.** Overall survival (OS) was the primary endpoint. Baseline differences were adjusted with inverse probability of treatment weighting (IPTW), which balanced covariates better than propensity score matching (all SMD < 0.1). Weighted Kaplan–Meier curves and Cox models gave pairwise hazard ratios. Subgroup analyses looked at AJCC stage and at the interval between systemic therapy and TACE (0–3, 4–14, 15–90 days; 4–21 / 22–90 days as a sensitivity check).

**Findings.**
- Systemic-first was associated with longer OS than TACE-first (weighted HR 0.80, 95% CI 0.72–0.88) and simultaneous treatment (HR 0.81, 0.70–0.95). Results held in stage III and stage IV subgroups.
- Among systemic-first patients, a 4–14 day gap before TACE had the best OS, compared with 0–3 days (HR 0.62, 0.50–0.77) and 15–90 days (HR 0.63, 0.46–0.87).

## Repository layout

Run the notebooks roughly in this order:

| File | Purpose |
|---|---|
| `data_preparation.ipynb` | Builds the analytic cohort from the raw NCDB liver file: HCC selection, TACE/systemic definitions, survival time/event, stage, sequencing groups. Writes `analysis_ready.csv`. |
| `psm_iptw_asmd.ipynb` | Builds pairwise datasets, fits IPTW weights, and compares covariate balance (ASMD) for IPTW vs PSM (Figure S1). |
| `TACE_PSM.Rmd` | Propensity score matching in R (`MatchIt`, `cobalt`) for each pairwise comparison; writes matched datasets used by the balance comparison. |
| `three_groups.ipynb` | Main three-group analysis: Table 1, unadjusted and IPTW-weighted Kaplan–Meier curves, pairwise Cox models (Figures 2–3). |
| `stage_spacing_subgroup.ipynb` | Stage-specific subgroup analyses and KM plots (Figure 4), plus exploratory interval × sequencing interactions. |
| `spacing_14.ipynb` | Interval optimization within the systemic-first group: 0–3 / 4–14 / 15–90 days with pairwise IPTW comparisons (Figure 5). Changing the bins to 3 / 21 / 90 gives the sensitivity grouping (Figure S4). |
| `MICE_all.ipynb`, `MICE_pred.ipynb` | Exploratory: multiple imputation (MICE) of missing stage, checked against held-out values. |
| `TACE_vs_sys.ipynb` | Exploratory: early two-group (TACE-first vs systemic-first) analysis, including bootstrap and interval checks. |

## Data

The NCDB Participant User File is available only to investigators at CoC-accredited programs under a data use agreement, so **no data are included in this repository**. All `.csv` files (raw and derived) are git-ignored. The preparation notebooks expect the raw file at `../NCDB_2022/liver.csv`. Every later step reads the intermediate files written by earlier steps in the repository root.

The NCDB is a joint project of the Commission on Cancer of the American College of Surgeons and the American Cancer Society. Neither organization has verified, and neither is responsible for, the analytic or statistical methods used here or the conclusions drawn.

## Environment

- Python 3.12 with `pandas`, `numpy`, `scikit-learn`, `lifelines`, `scipy`, `matplotlib`, `seaborn`
- R with `MatchIt`, `cobalt`, `readr` (for `TACE_PSM.Rmd` only)

## Authors

Qingyao Shang\*, Nicolas Zhang\*, Yi Yang\*, Meishuo Ouyang, Yinting Liu, Iris Luo, Janica Luo, Selina Yuan, Hui Zhang, You-Wen He†, Hong Zhao†

\* Equal contribution · † Corresponding authors
