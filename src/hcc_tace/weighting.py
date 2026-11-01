"""Propensity scores, IPTW weights and covariate balance."""

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from . import config


def perform_iptw(df, covariates=None, treatment_col="sequence_group"):
    """Fit a logistic propensity model and add `ps` and `iptw_weight` columns.

    The treatment column must be coded 0/1. Weights are 1/ps for treated (1)
    and 1/(1 - ps) for the reference group (0). Modifies and returns `df`.
    """
    covariates = config.COVARIATES if covariates is None else covariates
    X = pd.get_dummies(df[covariates], drop_first=True)
    y = df[treatment_col]
    ps_model = LogisticRegression(max_iter=2000)
    ps_model.fit(X, y)
    df["ps"] = ps_model.predict_proba(X)[:, 1]
    df["iptw_weight"] = df.apply(
        lambda row: 1 / row["ps"] if row[treatment_col] == 1 else 1 / (1 - row["ps"]),
        axis=1,
    )
    return df


def _weighted_mean_std(x, w):
    mean = np.average(x, weights=w)
    var = np.average((x - mean) ** 2, weights=w)
    return mean, np.sqrt(var)


def compute_asmd(df, covariate, treatment_col, weight_col=None):
    """Absolute standardized mean difference between the two treatment groups.

    Continuous covariates use the pooled SD. Categorical covariates are
    dummy-encoded and the ASMD is averaged across levels.
    """
    groups = df[treatment_col].unique()
    g1 = df[df[treatment_col] == groups[0]]
    g2 = df[df[treatment_col] == groups[1]]
    if weight_col is None:
        w1, w2 = np.ones(len(g1)), np.ones(len(g2))
    else:
        w1, w2 = g1[weight_col].values, g2[weight_col].values

    if df[covariate].dtype in [np.float64, np.int64]:
        m1, s1 = _weighted_mean_std(g1[covariate], w1)
        m2, s2 = _weighted_mean_std(g2[covariate], w2)
        return abs(m1 - m2) / np.sqrt((s1 ** 2 + s2 ** 2) / 2)

    dummies = pd.get_dummies(df[covariate], drop_first=False)
    asmd_list = []
    for level in dummies.columns:
        m1 = np.average(dummies.loc[g1.index, level], weights=w1)
        m2 = np.average(dummies.loc[g2.index, level], weights=w2)
        pooled = (m1 * (1 - m1) + m2 * (1 - m2)) / 2
        asmd_list.append(abs(m1 - m2) / np.sqrt(pooled) if pooled > 0 else 0)
    return np.mean(asmd_list)


def asmd_table(df, covariates=None, treatment_col="sequence_group", weight_col="iptw_weight"):
    """ASMD for every covariate before and after weighting."""
    covariates = config.COVARIATES if covariates is None else covariates
    unweighted, weighted = [], []
    for covariate in covariates:
        try:
            a1 = compute_asmd(df, covariate, treatment_col, weight_col=None)
            a2 = compute_asmd(df, covariate, treatment_col, weight_col=weight_col)
        except Exception as e:
            print(f"Skipping {covariate}: {e}")
            continue
        unweighted.append(a1)
        weighted.append(a2)
    out = pd.DataFrame({"Covariate": covariates, "Unweighted": unweighted, "IPTW": weighted})
    return out.sort_values(by="Unweighted", ascending=True)
