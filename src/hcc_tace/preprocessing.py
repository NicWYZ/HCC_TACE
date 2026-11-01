"""Cohort construction and covariate recoding for the NCDB liver PUF."""

import re

import numpy as np
import pandas as pd

from . import config


def select_hcc(df):
    """Keep HCC histologies (see config.HCC_HISTOLOGY_CODES)."""
    return df[df["HISTOLOGY"].isin(config.HCC_HISTOLOGY_CODES)].copy()


def select_tace_and_systemic(df):
    """Keep patients with both a TACE date and a systemic therapy date.

    TACE is identified from RX_SUMM_CHEMO codes 1-3 together with a
    chemotherapy start date; systemic therapy from the immunotherapy start date.
    """
    keep = (
        df["DX_CHEMO_STARTED_DAYS"].notna()
        & df["DX_IMMUNO_STARTED_DAYS"].notna()
        & df["RX_SUMM_CHEMO"].isin([1, 2, 3])
    )
    return df[keep].copy()


def select_known_survival(df):
    """Drop patients with unknown vital status or follow-up time."""
    keep = (
        df["PUF_VITAL_STATUS"].isin([0, 1])
        & df["DX_LASTCONTACT_DEATH_MONTHS"].notna()
        & (df["DX_LASTCONTACT_DEATH_MONTHS"] != 9999.0)
    )
    return df[keep].copy()


def add_survival_outcome(df):
    """time = months from diagnosis to death/last contact; event = 1 if dead."""
    df = df.copy()
    df["time"] = df["DX_LASTCONTACT_DEATH_MONTHS"]
    df["event"] = (df["PUF_VITAL_STATUS"] == 0).astype(int)
    return df


def add_stage(df):
    """COMBINED_STAGE_GROUP = analytic AJCC stage I-IV; other values dropped."""
    df = df.copy()
    stage = df["ANALYTIC_STAGE_GROUP"]
    df["COMBINED_STAGE_GROUP"] = stage.where(stage.isin([1, 2, 3, 4]), np.nan)
    df = df.dropna(subset=["COMBINED_STAGE_GROUP"])
    df["COMBINED_STAGE_GROUP"] = df["COMBINED_STAGE_GROUP"].astype(int)
    return df


def add_sequence_group(df, window=config.SIMULTANEOUS_WINDOW_DAYS):
    """Classify treatment order.

    0 = TACE first (TACE more than `window` days before systemic therapy)
    1 = systemic first (systemic therapy more than `window` days before TACE)
    2 = simultaneous (within `window` days)
    """
    df = df.copy()
    df["treatment_spacing_tace_ref"] = df["DX_CHEMO_STARTED_DAYS"] - df["DX_IMMUNO_STARTED_DAYS"]
    df["sequence_group"] = df["treatment_spacing_tace_ref"].apply(
        lambda x: config.TACE_FIRST if x < -window
        else (config.SYSTEMIC_FIRST if x > window else config.SIMULTANEOUS)
    )
    return df


def _yes_no(mask):
    return mask.astype(int).map({0: "No", 1: "Yes"})


def add_surgery_flags(df):
    """Surgery of the primary site, and its timing relative to TACE/systemic."""
    df = df.copy()
    surg = df["DX_DEFSURG_STARTED_DAYS"]
    tace = df["DX_CHEMO_STARTED_DAYS"]
    sys_ = df["DX_IMMUNO_STARTED_DAYS"]
    df["Surgery"] = _yes_no(~df["RX_SUMM_SURG_PRIM_SITE"].isin([0, 99]))
    df["Prior_Surgery"] = _yes_no((surg < tace) & (surg < sys_))
    df["Subsequent_Surgery"] = _yes_no((surg > tace) & (surg > sys_))
    df["Surgery_Between"] = _yes_no(
        ((surg >= tace) & (surg <= sys_)) | ((surg >= sys_) & (surg <= tace))
    )
    return df


# --- covariate recoding -----------------------------------------------------

def collapse_year(year):
    if year >= 2020:
        return "2020+"
    elif year >= 2017:
        return "2017–2019"
    else:
        return "≤2016"


def recode_hispanic_origin(x):
    if x == 0:
        return "Non-Hispanic"
    elif x in [1, 2, 3, 4, 5, 6, 8]:
        return "Hispanic"
    elif x == 7:
        return "Likely Hispanic"
    else:
        return np.nan


def bin_tumor_size(x):
    if pd.isna(x):
        return np.nan
    elif x == 0:
        return "No tumor"
    elif x <= 30:
        return "≤3 cm"
    elif x <= 50:
        return "3–5 cm"
    elif x <= 100:
        return "5–10 cm"
    elif x <= 200:
        return "10–20 cm"
    else:
        return ">20 cm"


RACE_MAP = {
    "1": "White", "2": "Black", "3": "Other", "7": "Other", "27": "Other",
    "97": "Other", "98": "Other", "99": np.nan,
    **{code: "Asian" for code in ["4", "5", "6", "8", "10", "11", "13", "14", "15", "16", "17", "96"]},
}

INSURANCE_MAP = {
    0: "No insurance", 1: "Private/Managed Care", 2: "Medicaid",
    3: "Medicare", 4: "Other Government", 9: "Unknown",
}

MEDIAN_INCOME_MAP = {
    1: "<$46,277", 2: "$46,277-$57,856", 3: "$57,857-$74,062", 4: ">$74,063",
}

NO_HIGH_SCHOOL_MAP = {1: ">=15.3%", 2: "9.1%-15.2%", 3: "5.0%-9.0%", 4: "<5.0%"}

FACILITY_TYPE_MAP = {
    1: "Community Cancer Program",
    2: "Comprehensive Community Cancer Program",
    3: "Academic/Research Program",
    4: "Integrated Network Cancer Program",
}

HISTOLOGY_MAP = {
    8170: "HCC, NOS",
    8171: "Fibrolamellar HCC",
    8172: "Scirrhous HCC",
    8173: "Spindle cell HCC",
    8174: "Clear cell HCC",
    8175: "Pleomorphic HCC",
    8180: "Mixed HCC-Cholangiocarcinoma",
}

REGION_MAP = {
    1: "New England", 2: "Middle Atlantic", 3: "South Atlantic",
    4: "East North Central", 5: "East South Central", 6: "West North Central",
    7: "West South Central", 8: "Mountain", 9: "Pacific",
}


def _first_number(x):
    m = re.search(r"\d+", x)
    return m.group() if m else np.nan


def recode_covariates(df):
    """Map NCDB codes to labelled analysis covariates."""
    df = df.copy()

    df["COMBINED_CLIN_T"] = (
        df["AJCC_TNM_CLIN_T"].combine_first(df["TNM_CLIN_T"]).astype(str).apply(_first_number)
    )
    df["YEAR_OF_DIAGNOSIS"] = df["YEAR_OF_DIAGNOSIS"].apply(collapse_year)
    df["SEX"] = df["SEX"].replace({1: "Male", 2: "Female"})
    df["RACE"] = df["RACE"].astype(str).map(RACE_MAP).fillna("Other")
    df["INSURANCE_STATUS"] = df["INSURANCE_STATUS"].map(INSURANCE_MAP)
    df["MEDIAN_INCOME"] = df["MED_INC_QUAR_2020"].map(MEDIAN_INCOME_MAP)
    df["Charlson_Deyo_Score"] = df["CDCC_TOTAL_BEST"].map({0: "0", 1: "1", 2: "2", 3: "3"})
    df["NoHighSchooldegree2016_2020"] = df["NO_HSD_QUAR_2020"].map(NO_HIGH_SCHOOL_MAP)
    df["FACILITY_TYPE"] = df["FACILITY_TYPE_CD"].map(FACILITY_TYPE_MAP)
    df["HISTOLOGY"] = df["HISTOLOGY"].map(HISTOLOGY_MAP)
    df["FACILITY_REGION"] = df["FACILITY_LOCATION_CD"].map(REGION_MAP)
    df["HISPANIC"] = df["SPANISH_HISPANIC_ORIGIN"].apply(recode_hispanic_origin)

    size = pd.to_numeric(df["TUMOR_SIZE_SUMMARY_2016"], errors="coerce")
    df["TUMOR_SIZE_SUMMARY_2016"] = size.mask(size > 300)
    df["TUMOR_SIZE"] = df["TUMOR_SIZE_SUMMARY_2016"].apply(bin_tumor_size)

    df["Radiation_Therapy"] = df["PHASE_I_RT_MODALITY"].apply(
        lambda x: "No" if x == 0 else ("Unknown" if x == 99 else "Yes")
    )
    return df
