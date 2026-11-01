"""Cox models, log-rank tests and KM summaries for treatment comparisons."""

import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test

from . import config
from .plotting import DEFAULT, plot_km_three_groups, plot_km_two_groups

LANDMARKS = [36, 60]  # months


def risk_counts(df, time_points, group_value, group_col="sequence_group", weight_col=None):
    """Number at risk at each time point (IPTW-weighted if `weight_col` is given)."""
    group = df[df[group_col] == group_value]
    if weight_col is None:
        return [group[group["time"] >= t].shape[0] for t in time_points]
    return [int(round(group[group["time"] >= t][weight_col].sum())) for t in time_points]


def survival_with_ci(kmf, time):
    """KM survival at `time` with its 95% CI (None if not estimable)."""
    surv = kmf.predict(time)
    closest_idx = kmf.confidence_interval_.index.get_indexer([time], method="nearest")[0]
    ci = kmf.confidence_interval_.iloc[closest_idx]
    lower, upper = ci.iloc[0], ci.iloc[1]
    if (upper == 1.0 and lower == 1.0) or (upper == 0.0 and lower == 0.0):
        return surv, None, None
    return surv, lower, upper


def stratified_logrank_p(df, group_col, strata):
    """Stratified log-rank p-value via R's survival::survdiff (for matched sets)."""
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "survdiff_input.csv"
        df[["time", "event", group_col, strata]].to_csv(path, index=False)
        r_cmd = (
            f"library(survival); df <- read.csv('{path.as_posix()}'); "
            f"res <- survdiff(Surv(time,event) ~ {group_col} + strata({strata}), data=df); "
            "cat(res$chisq, '\\n', pchisq(res$chisq, df=1, lower.tail=FALSE), '\\n')"
        )
        proc = subprocess.run(["Rscript", "-e", r_cmd], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Rscript failed: {proc.stderr.strip()}")
    out_lines = [ln.strip() for ln in proc.stdout.strip().splitlines() if ln.strip()]
    if len(out_lines) < 2:
        raise RuntimeError(f"Unexpected R output: {proc.stdout!r}")
    return round(float(out_lines[1]), 4)


def survival_analysis(df, time_grid, groups, colors, group_col="sequence_group", strata=None,
                      weight_col=None, title=None, style=DEFAULT, p_format=None,
                      print_summary=True, print_landmarks=True, file_name=None):
    """Two-group comparison: Cox HR, log-rank p, and a KM figure.

    `group_col` must be coded 0 (reference, groups[0]) and 1 (groups[1]).
    Use `weight_col='iptw_weight'` for the IPTW analysis, or `strata='subclass'`
    for a PSM-matched dataset (stratified Cox and stratified log-rank).
    `p_format` is an optional format string for the p-value shown on the plot
    (default: rounded to 4 decimals).
    """
    cols = ["time", "event", group_col]
    cph = CoxPHFitter()
    if weight_col:
        cph.fit(df[cols + [weight_col]], duration_col="time", event_col="event", weights_col=weight_col)
    elif strata:
        cph.fit(df[cols + [strata]], duration_col="time", event_col="event", strata=[strata])
    else:
        cph.fit(df[cols], duration_col="time", event_col="event")
    if print_summary:
        cph.print_summary()
    row = cph.summary.loc[group_col]
    hr, ci_low, ci_high = row["exp(coef)"], row["exp(coef) lower 95%"], row["exp(coef) upper 95%"]

    mask_0, mask_1 = df[group_col] == 0, df[group_col] == 1
    if strata:
        p_val = stratified_logrank_p(df, group_col, strata)
    else:
        result = logrank_test(
            df.loc[mask_0, "time"], df.loc[mask_1, "time"],
            event_observed_A=df.loc[mask_0, "event"], event_observed_B=df.loc[mask_1, "event"],
            weights_A=df.loc[mask_0, weight_col] if weight_col else None,
            weights_B=df.loc[mask_1, weight_col] if weight_col else None,
        )
        p_val = format(result.p_value, p_format) if p_format else result.p_value.round(4)

    kmfs = []
    for mask, label in [(mask_0, groups[0]), (mask_1, groups[1])]:
        kmf = KaplanMeierFitter()
        kmf.fit(durations=df.loc[mask, "time"], event_observed=df.loc[mask, "event"],
                weights=df.loc[mask, weight_col] if weight_col else None, label=label)
        kmfs.append(kmf)
    kmf_0, kmf_1 = kmfs
    risk_0 = risk_counts(df, time_grid, 0, group_col, weight_col)
    risk_1 = risk_counts(df, time_grid, 1, group_col, weight_col)

    if print_landmarks:
        for t in LANDMARKS:
            print(f"\nTime: {t} months")
            for label, kmf in [(groups[0], kmf_0), (groups[1], kmf_1)]:
                surv, low, high = survival_with_ci(kmf, t)
                if low is None or high is None:
                    print(f"{label}: {surv:.1%} (CI not estimable)")
                else:
                    print(f"{label}: {surv:.1%} ({low:.1%} – {high:.1%})")

    surv_text = {label: {t: f"{kmf.predict(t):.1%}" for t in LANDMARKS}
                 for label, kmf in [(groups[0], kmf_0), (groups[1], kmf_1)]}
    plot_km_two_groups(groups, colors, kmf_0, kmf_1, hr, ci_low, ci_high, p_val, surv_text,
                       risk_0, risk_1, time_grid, weighted=bool(weight_col), title=title,
                       style=style, file_name=file_name)
    return {"hr": hr, "ci_low": ci_low, "ci_high": ci_high, "p": p_val}


REFERENCE_LABELS = {
    "TACE": ["TACE First", "Systemic First", "Simultaneous"],
    "Systemic": ["Systemic First", "TACE First", "Simultaneous"],
    "Simultaneous": ["Simultaneous", "TACE First", "Systemic First"],
}


def reference_labels(ref_group, simultaneous_label="Simultaneous"):
    """Group labels in 0/1/2 order for a sequence coding with `ref_group` as 0."""
    key = "Simultaneous" if ref_group == simultaneous_label else ref_group
    if key not in REFERENCE_LABELS:
        raise ValueError(f"Unknown reference group: {ref_group}")
    return [simultaneous_label if lbl == "Simultaneous" else lbl for lbl in REFERENCE_LABELS[key]]


def survival_analysis_3g(df, time_grid, group_col, labels=None, ref_group=None, title=None,
                         simultaneous_label="Simultaneous", dummy_prefix="seq", style=DEFAULT,
                         print_landmarks=True, file_name=None):
    """Unadjusted three-group Cox model and KM figure.

    `group_col` is coded 0/1/2 with 0 as the reference. Pass `labels` directly,
    or `ref_group` ('TACE', 'Systemic' or the simultaneous label) to derive them.
    """
    if labels is None:
        labels = reference_labels(ref_group, simultaneous_label)
    if ref_group:
        plot_title = f"{title}: {ref_group} as reference" if title else f"{ref_group} as reference"
    else:
        plot_title = f"{title}"

    dummies = pd.get_dummies(df[group_col], prefix=dummy_prefix, drop_first=True)
    cph = CoxPHFitter()
    cph.fit(pd.concat([df[["time", "event"]], dummies], axis=1), duration_col="time", event_col="event")
    cph.print_summary()

    kmf_dict, risk_dict = {}, {}
    for val, lbl in zip([0, 1, 2], labels):
        mask = df[group_col] == val
        kmf = KaplanMeierFitter(label=lbl)
        kmf.fit(df.loc[mask, "time"], event_observed=df.loc[mask, "event"])
        kmf_dict[lbl] = kmf
        risk_dict[lbl] = risk_counts(df, time_grid, val, group_col)
        if print_landmarks:
            for t in LANDMARKS:
                lower, upper = None, None
                ci = kmf.confidence_interval_
                if not ci.empty:
                    nearest = ci.index.get_indexer([t], method="nearest")[0]
                    lower, upper = ci.iloc[nearest, 0], ci.iloc[nearest, 1]
                surv = kmf.predict(t)
                if lower is None or upper is None:
                    print(f"{lbl} at {t} mo: {surv:.1%} (CI not estimable)")
                else:
                    print(f"{lbl} at {t} mo: {surv:.1%} ({lower:.1%}-{upper:.1%})")

    plot_km_three_groups(kmf_dict, risk_dict, time_grid, plot_title, style=style, file_name=file_name)
    return cph


def multigroup_logrank_p(df, group_col, values=(0, 1, 2)):
    """Overall log-rank p-value across the levels of `group_col`."""
    parts = [df[df[group_col] == v] for v in values]
    result = multivariate_logrank_test(
        np.concatenate([p["time"].values for p in parts]),
        np.concatenate([np.full(len(p), i) for i, p in enumerate(parts)]),
        np.concatenate([p["event"].values for p in parts]),
    )
    return result.p_value


def weighted_cox_hr(df, group_col="sequence_group", weight_col="iptw_weight"):
    """IPTW-weighted Cox HR, 95% CI and p-value for a 0/1 group column."""
    cph = CoxPHFitter()
    cph.fit(df[["time", "event", group_col, weight_col]], duration_col="time", event_col="event",
            weights_col=weight_col)
    row = cph.summary.loc[group_col]
    return {"HR": row["exp(coef)"], "CI Low": row["exp(coef) lower 95%"],
            "CI High": row["exp(coef) upper 95%"], "p-value": row["p"].round(4), "n": df.shape[0]}


def stage_sequence_interaction(df, categories):
    """Weighted Cox models testing a sequence x stage interaction.

    Fits the interaction with stage as numeric, then with stage as categorical
    dummies (first entry of `categories` is the reference stage).
    """
    df = df.copy()
    cph_1 = CoxPHFitter()
    cph_1.fit(df, duration_col="time", event_col="event", weights_col="iptw_weight",
              formula="sequence_group * COMBINED_STAGE_GROUP")
    cph_1.print_summary()

    df["COMBINED_STAGE_GROUP"] = pd.Categorical(df["COMBINED_STAGE_GROUP"].astype(int),
                                                categories=categories, ordered=False)
    df = pd.get_dummies(df, columns=["COMBINED_STAGE_GROUP"], drop_first=True)
    stage_cols = [f"COMBINED_STAGE_GROUP_{cat}" for cat in categories[1:]]
    for col in stage_cols:
        df[f"sequence_group:{col}"] = df["sequence_group"] * df[col]
    formula = "sequence_group + " + " + ".join(stage_cols + [f"sequence_group:{c}" for c in stage_cols])
    cph_2 = CoxPHFitter()
    cph_2.fit(df, duration_col="time", event_col="event", weights_col="iptw_weight", formula=formula)
    cph_2.print_summary()


def best_logrank_cutpoint(df, col="treatment_spacing", min_frac=0.25):
    """Cut point on `col` that maximises the two-group log-rank statistic.

    Groups smaller than `min_frac` of the sample are not considered.
    """
    best_cut, best_stat = None, 0
    for cut in range(1, int(df[col].max())):
        low, high = df[col] <= cut, df[col] > cut
        if low.sum() < df.shape[0] * min_frac or high.sum() < df.shape[0] * min_frac:
            continue
        result = logrank_test(df.loc[low, "time"], df.loc[high, "time"],
                              event_observed_A=df.loc[low, "event"], event_observed_B=df.loc[high, "event"])
        if result.test_statistic > best_stat:
            best_cut, best_stat = cut, result.test_statistic
    return best_cut, best_stat
