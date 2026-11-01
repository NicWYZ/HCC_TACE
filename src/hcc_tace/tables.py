"""Baseline characteristics (Table 1)."""

import re

import numpy as np
import pandas as pd
import scipy.stats as stats


def _fmt_p(p):
    return f"{p:.4f}" if p >= 0.0001 else "<0.0001"


def _sorted_levels(values, order=None):
    values = values.dropna().unique()
    if order is not None:
        return sorted(values, key=lambda x: order.index(x) if x in order else float("inf"))
    try:
        return sorted(values, key=lambda x: int(re.findall(r"\d+", str(x))[0])
                      if re.search(r"\d+", str(x)) else str(x))
    except TypeError:
        return sorted(values)


def baseline_table(df, covariates, group_col, group_labels, level_order=None, covariate_order=None):
    """Counts (%) or mean (SD) per group with chi-square / ANOVA p-values.

    Columns with object dtype or <= 10 distinct values are treated as
    categorical. Returns a DataFrame indexed by (Covariate, Level).
    """
    level_order = level_order or {}
    group = df[group_col].map(group_labels)
    records = []

    for col in covariates:
        p_val = ""
        if df[col].dtype == "object" or df[col].nunique() <= 10:
            contingency = pd.crosstab(df[col], group)
            if contingency.shape[0] > 1 and contingency.shape[1] > 1:
                _, p, _, _ = stats.chi2_contingency(contingency)
                p_val = _fmt_p(p)
            for i, lvl in enumerate(_sorted_levels(df[col], level_order.get(col))):
                record = {}
                for grp in group_labels.values():
                    mask = group == grp
                    total = mask.sum()
                    count = (df.loc[mask, col] == lvl).sum()
                    pct = 100 * count / total if total > 0 else 0.0
                    record[grp] = f"{count} ({pct:.1f}%)"
                record["p-value"] = p_val if i == 0 else ""
                records.append(((col, lvl), record))
        elif np.issubdtype(df[col].dtype, np.number):
            samples = [df.loc[group == grp, col].dropna() for grp in group_labels.values()]
            if all(len(s) > 1 for s in samples):
                _, p = stats.f_oneway(*samples)
                p_val = _fmt_p(p)
            record = {}
            for grp in group_labels.values():
                values = df.loc[group == grp, col]
                record[grp] = f"{values.mean():.2f} ({values.std():.2f})"
            record["p-value"] = p_val
            records.append(((col, "mean (SD)"), record))

    number = {grp: f"{(group == grp).sum()}" for grp in group_labels.values()}
    number["p-value"] = ""
    records.insert(0, (("Number", ""), number))

    index, data = zip(*records)
    table = pd.DataFrame(list(data), index=pd.MultiIndex.from_tuples(index, names=["Covariate", "Level"]),
                         columns=list(group_labels.values()) + ["p-value"])
    if covariate_order:
        table = table.loc[[idx for cov in covariate_order for idx in table.index if idx[0] == cov]]
    return table
