# Sequencing and Timing of TACE and Systemic Therapy in Advanced Hepatocellular Carcinoma

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)
![R 4.5](https://img.shields.io/badge/R-4.5-276DC3.svg)
![Data: NCDB](https://img.shields.io/badge/data-NCDB%20PUF-lightgrey.svg)

Analysis code for:

> **Optimizing the Sequencing and Timing of Transarterial Chemoembolization and Systemic Therapy in Advanced Hepatocellular Carcinoma**
> Qingyao Shang\*, Nicolas Zhang\*, Yi Yang\*, Meishuo Ouyang, Yinting Liu, Iris Luo, Janica Luo, Selina Yuan, Hui Zhang, You-Wen He†, Hong Zhao†
> \*Equal contribution · †Corresponding authors

**Status:** manuscript submitted to *eBioMedicine*. Selected for a **poster presentation at the 2026 ASCO Gastrointestinal Cancers Symposium (ASCO GI 2026)**.

---

## Overview

TACE and systemic therapy work together in advanced HCC (e.g. the TACTICS and LAUNCH trials). The best **order** for the two treatments, and the best **interval** between them, are not known. This retrospective cohort study uses the National Cancer Database (NCDB) to compare three strategies in AJCC stage III–IV HCC:

| Strategy | Definition |
|---|---|
| **TACE first** | TACE, then systemic therapy more than 3 days later |
| **Systemic first** | Systemic therapy, then TACE more than 3 days later |
| **Simultaneous** | Both started within 3 days |

The primary endpoint is overall survival (OS). Confounding is handled with inverse probability of treatment weighting (IPTW); propensity score matching (PSM) is a comparison method.

## Key findings

The numbers below are the ones this repository reproduces (stage III–IV, IPTW-weighted Cox models). HR < 1 favours the second group named.

| Comparison | HR (95% CI) |
|---|---|
| Systemic first vs TACE first | 0.81 (0.74–0.90) |
| Systemic first vs simultaneous | 0.83 (0.71–0.96) |
| Systemic first: 4–14 days vs 0–3 days before TACE | 0.61 (0.49–0.75) |
| Systemic first: 4–14 days vs 15–90 days before TACE | 0.56 (0.41–0.77) |

**Conclusion:** starting systemic therapy first and following with TACE within 4–14 days gave the longest survival. The finding held in both stage III and stage IV subgroups. Prospective validation is needed.

## Repository structure

```
HCC_TACE/
├── src/hcc_tace/             # Reusable analysis package
│   ├── config.py             #   paths, covariates, group coding, plot colours
│   ├── preprocessing.py      #   cohort selection and NCDB code recoding
│   ├── cohort.py             #   advanced-stage and pairwise subsets
│   ├── weighting.py          #   propensity scores, IPTW, ASMD
│   ├── survival.py           #   Cox / log-rank / KM analyses
│   ├── plotting.py           #   KM figures with risk tables, balance plots
│   └── tables.py             #   Table 1
├── notebooks/                # Analysis, run in numbered order
│   ├── 01_cohort_construction.ipynb
│   ├── 02_covariate_balance.ipynb
│   ├── 03_tace_vs_systemic.ipynb
│   ├── 04_three_group_comparison.ipynb
│   ├── 05_stage_subgroups.ipynb
│   ├── 06_interval_optimization.ipynb
│   └── exploratory/          # MICE imputation of missing stage (not in the paper)
├── R/psm_matching.R          # Propensity score matching (MatchIt)
├── data/                     # NCDB data, not tracked (see data/README.md)
├── results/                  # Generated figures and tables, not tracked
├── environment.yml           # Conda environment (Python + R)
├── requirements.txt          # Pinned Python dependencies
├── pyproject.toml            # Installs src/hcc_tace
└── Makefile                  # Runs the full pipeline
```

## Notebooks and manuscript outputs

| Step | Notebook / script | Produces |
|---|---|---|
| 1 | `01_cohort_construction` | Analytic cohort and patient flow counts (Figure 1); pairwise cohorts for matching |
| 2 | `R/psm_matching.R` | 1:4, 1:3 and 1:1 nearest-neighbour matched cohorts (caliper 0.2 SD), love plots |
| 3 | `02_covariate_balance` | ASMD before adjustment, after IPTW and after PSM (Figure S1) |
| 4 | `03_tace_vs_systemic` | TACE first vs systemic first: unadjusted, PSM and IPTW KM curves (Figure 3); stage-combination subgroups |
| 5 | `04_three_group_comparison` | Table 1; sequencing by stage (Figure S2); three-group KM curves (Figure 2); comparisons with the simultaneous group (Figure 3) |
| 6 | `05_stage_subgroups` | Sequence × stage interaction models; stage-specific IPTW KM curves (Figure 4) |
| 7 | `06_interval_optimization` | Interval groups 0–3 / 4–14 / 15–90 days in the systemic-first group (Figure 5); set the bins to `[-1, 3, 21, 90]` for the 4–21 day sensitivity analysis (Figure S4) |

The notebooks are saved with their outputs, so you can read the results without the data. All outputs are aggregate statistics.

## Methods summary

- **Data source.** NCDB liver Participant User File, 2004–2022 release. The cohort includes patients diagnosed from 2007 onward.
- **Population.** Adults with HCC (ICD-O-3 histology 8170, 8172–8175; fibrolamellar and combined HCC-cholangiocarcinoma excluded) who received both TACE and systemic therapy, with known survival and AJCC stage. Main analyses use stage III–IV.
- **Exposure.** TACE was identified from `RX_SUMM_CHEMO` together with the chemotherapy start date. Systemic therapy (targeted agents, immunotherapy) was identified from the immunotherapy start date. Groups were defined from the difference between the two start dates, using a ±3-day window for "simultaneous".
- **Covariates.** Age, sex, race, insurance status, Charlson-Deyo score, histology, AJCC stage, year of diagnosis, surgery and radiation therapy.
- **Confounding adjustment.** Logistic-regression propensity scores, with IPTW weights of 1/PS for the treated group and 1/(1−PS) for the reference group. Balance target: ASMD < 0.1. PSM used `MatchIt` nearest-neighbour matching with a 0.2 SD caliper.
- **Outcome models.** Kaplan–Meier curves with (weighted) log-rank tests and Cox proportional hazards models. Matched analyses use stratified Cox models and a stratified log-rank test (`survival::survdiff`). Two-sided α = 0.05.

## Reproducing the analysis

**1. Get the data.** Request the NCDB liver PUF through your CoC-accredited program and place it at `data/raw/liver.csv`, or set `NCDB_LIVER_CSV` to its path.

**2. Create the environment.**

```bash
conda env create -f environment.yml
```

```bash
conda activate hcc-tace
```

```bash
pip install -e .
```

To use your own Python environment instead, run `pip install -r requirements.txt` and `pip install -e .`, then install R with `MatchIt`, `cobalt`, `readr`, `ggplot2` and `survival`.

**3. Run the pipeline** from the repository root:

```bash
make all
```

This runs cohort construction, then R matching, then the analysis notebooks. Figures are written to `results/figures/` and Table 1 to `results/tables/`. Results are deterministic, and the bootstrap CIs for the exploratory spline curves are seeded.

## Software

Python 3.12 (pandas 2.2, NumPy 2.0, scikit-learn 1.6, lifelines 0.30, SciPy 1.15, matplotlib 3.10, seaborn 0.13). R 4.5 (MatchIt 4.7, cobalt 4.6, survival).

## Citation

The manuscript is under review. Until it is published, please cite this repository (see [`CITATION.cff`](CITATION.cff)):

```bibtex
@misc{shang_zhang_yang_hcc_tace,
  title  = {Optimizing the Sequencing and Timing of Transarterial Chemoembolization and
            Systemic Therapy in Advanced Hepatocellular Carcinoma},
  author = {Shang, Qingyao and Zhang, Nicolas and Yang, Yi and Ouyang, Meishuo and Liu, Yinting and
            Luo, Iris and Luo, Janica and Yuan, Selina and Zhang, Hui and He, You-Wen and Zhao, Hong},
  note   = {Manuscript submitted to eBioMedicine. Code: https://github.com/NicWYZ/HCC_TACE},
  year   = {2026}
}
```

## Data availability and disclaimer

NCDB data cannot be redistributed, so this repository contains no patient-level data. The NCDB is a joint project of the American College of Surgeons Commission on Cancer and the American Cancer Society. Neither organization has verified, and neither is responsible for, the analytic or statistical methods used here or the conclusions drawn.

## License

Code is released under the [MIT License](LICENSE).
