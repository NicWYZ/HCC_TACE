"""Project paths and analysis constants shared by all notebooks."""

import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]

DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
TABLES_DIR = RESULTS_DIR / "tables"

# NCDB liver PUF. Override with the NCDB_LIVER_CSV environment variable.
RAW_NCDB_FILE = Path(os.environ.get("NCDB_LIVER_CSV", RAW_DIR / "liver.csv"))
ANALYSIS_READY = PROCESSED_DIR / "analysis_ready.csv"

# ICD-O-3 histology codes kept as HCC (8171 fibrolamellar and 8180 combined
# HCC-cholangiocarcinoma are excluded).
HCC_HISTOLOGY_CODES = [8170, 8172, 8173, 8174, 8175]

# Covariates used in the propensity score model.
COVARIATES = [
    "AGE", "SEX", "RACE",
    "INSURANCE_STATUS",
    "Charlson_Deyo_Score", "HISTOLOGY", "COMBINED_STAGE_GROUP",
    "YEAR_OF_DIAGNOSIS", "Surgery", "Radiation_Therapy",
]

# sequence_group coding written by 01_cohort_construction.
TACE_FIRST, SYSTEMIC_FIRST, SIMULTANEOUS = 0, 1, 2
SIMULTANEOUS_WINDOW_DAYS = 3

ADVANCED_STAGES = [3, 4]
TIME_GRID = np.arange(0, 73, 12)  # months

COLORS = {"tace": "#ff7f0e", "systemic": "#1f77b4", "simultaneous": "#2ca02c"}


def setup_notebook():
    """Create output folders and apply shared plotting/warning settings."""
    from lifelines.exceptions import StatisticalWarning

    for d in (PROCESSED_DIR, FIGURES_DIR, TABLES_DIR):
        d.mkdir(parents=True, exist_ok=True)
    warnings.filterwarnings("ignore", category=StatisticalWarning)
    plt.rcParams["pdf.fonttype"] = 42  # editable text in exported PDFs


def figure_path(file_name, suffix=".pdf"):
    """Resolve a bare figure name to results/figures/<name>.pdf."""
    path = Path(f"{file_name}{suffix}")
    return path if path.is_absolute() else FIGURES_DIR / path
