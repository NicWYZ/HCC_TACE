# Data

Nothing in this folder is tracked by git. NCDB data are available only to investigators at
Commission on Cancer–accredited programs under the NCDB Participant User File (PUF) data use agreement.

```
data/
├── raw/
│   └── liver.csv          # NCDB liver PUF (2004–2022), as delivered
└── processed/             # written by the pipeline
    ├── analysis_ready.csv             # 01_cohort_construction
    ├── df_{tace_sys,tace_same,same_sys}.csv      # 01, input to R/psm_matching.R
    └── matched_{tace_sys,tace_same,same_sys}.csv # R/psm_matching.R
```

To keep the PUF elsewhere, set `NCDB_LIVER_CSV=/path/to/liver.csv` before running the notebooks.
