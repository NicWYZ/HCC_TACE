"""Kaplan-Meier figures with number-at-risk tables, and IPTW diagnostics."""

from dataclasses import dataclass, replace

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from lifelines import KaplanMeierFitter
from matplotlib.lines import Line2D

from . import config
from .weighting import asmd_table


@dataclass(frozen=True)
class KMStyle:
    """Font sizes and options for the KM figures.

    `KMStyle()` is the style of the main pairwise / three-group figures;
    `LARGE` is the larger-font style used for the stage and interval panels.
    """
    label_size: int = 20          # axis labels (two-group plot)
    title_size: int = 20
    legend_size: int = 20
    tick_size: int = 12
    table_size: int = 20          # number-at-risk table
    annotate: bool = True         # print 3- and 5-year OS on the curves
    annot_offsets: tuple = (-0.095, 0.075)  # y offsets for group 0 / group 1 labels
    # three-group plot
    label_size_3g: int = 18
    title_size_3g: int = 18
    tick_size_3g: int = 12
    y_tick_size_3g: int = None
    table_size_3g: int = 18
    months_label_3g: tuple = (16, 0.3)  # (font size, y position) of the lower "Months" label
    legend_3g: bool = True


DEFAULT = KMStyle()
LARGE = KMStyle(
    label_size=24, title_size=24, legend_size=24, tick_size=20, table_size=24, annotate=False,
    label_size_3g=20, title_size_3g=20, tick_size_3g=20, y_tick_size_3g=20, table_size_3g=19,
    months_label_3g=(20, 0.275),
)


def _save(file_name):
    if file_name:
        plt.savefig(config.figure_path(file_name), format="pdf", bbox_inches="tight")


def plot_km_two_groups(groups, colors, kmf_0, kmf_1, hr, ci_low, ci_high, p_val, surv_text,
                       risk_0, risk_1, time_grid, weighted=False, title=None, style=DEFAULT,
                       file_name=None):
    """Two-group KM curves with HR, log-rank p and a number-at-risk table."""
    fig = plt.figure(figsize=(9, 7))
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 0.75], hspace=0.05)
    ax_km = plt.subplot(gs[0])
    kmf_0.plot(ax=ax_km, ci_show=True, color=colors[0], linewidth=3)
    kmf_1.plot(ax=ax_km, ci_show=True, color=colors[1], linewidth=3)

    if style.annotate:
        off_0, off_1 = style.annot_offsets
        for x in [36, 60]:
            ax_km.text(x, kmf_1.predict(x) + off_1, surv_text[groups[1]][x], color=colors[1],
                       fontsize=20, fontweight="bold", ha="center")
            ax_km.text(x, kmf_0.predict(x) + off_0, surv_text[groups[0]][x], color=colors[0],
                       fontsize=20, fontweight="bold", ha="center")
    else:
        ax_km.get_legend().remove()

    custom_lines = [Line2D([0], [0], color=colors[0], lw=3), Line2D([0], [0], color=colors[1], lw=3)]
    ax_km.legend(custom_lines, groups, loc="lower left", frameon=False, fontsize=style.legend_size)
    hr_text = f"HR = {hr:.2f}\n(95% CI = {ci_low:.2f} – {ci_high:.2f})\n\nlogrank p={p_val}"
    ax_km.text(0.67, 0.75, hr_text, transform=ax_km.transAxes, fontsize=14, va="bottom", ha="left")
    ax_km.set_title(title or ("IPTW-Weighted" if weighted else "Unweighted"),
                    fontsize=style.title_size, fontweight="bold")
    ax_km.set_ylabel("Overall Survival Rate (%)", fontsize=style.label_size, fontweight="bold")
    ax_km.set_xlabel("Months", fontsize=style.label_size, fontweight="bold")
    ax_km.set_xlim(0, max(time_grid))
    ax_km.set_xticks(time_grid)
    ax_km.tick_params(labelsize=style.tick_size)
    ax_km.spines["top"].set_visible(False)
    ax_km.spines["right"].set_visible(False)

    ax_table = plt.subplot(gs[1])
    ax_table.axis("off")
    table_data = [[""] + [str(n) for n in risk_0[:len(time_grid)]],
                  [""] + [str(n) for n in risk_1[:len(time_grid)]]]
    risk_table = ax_table.table(cellText=table_data, cellLoc="center", loc="center", edges="open",
                                bbox=[-0.245, -0.4, 1.325, 0.6])
    risk_table.auto_set_font_size(False)
    risk_table.set_fontsize(style.table_size)
    risk_table.scale(1.2, 1.3)
    for row_idx, color in enumerate(colors[:2]):
        for col_idx in range(len(time_grid) + 1):
            cell = risk_table[row_idx, col_idx]
            cell.set_text_props(color=color)
            cell.set_linewidth(0)
    ax_table.text(0, 0.4, "Weighted No. at Risk" if weighted else "Number at Risk", fontsize=12,
                  transform=ax_table.transAxes, ha="center")
    ax_table.plot([-0.075, -0.05], [0.05, 0.05], color=colors[0], lw=4,
                  transform=ax_table.transAxes, clip_on=False)
    ax_table.plot([-0.075, -0.05], [-0.25, -0.25], color=colors[1], lw=4,
                  transform=ax_table.transAxes, clip_on=False)
    plt.subplots_adjust(top=0.93, bottom=0.07)
    _save(file_name)
    plt.show()


def plot_km_three_groups(kmf_dict, risk_dict, time_grid, title, style=DEFAULT, file_name=None):
    """Three-group KM curves (no CI bands) with a number-at-risk table."""
    labels = list(kmf_dict.keys())
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    fig = plt.figure(figsize=(9, 7))
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1.5], hspace=0.2)
    ax_km = fig.add_subplot(gs[0])
    ax_table = fig.add_subplot(gs[1], sharex=ax_km)
    ax_table.axis("off")
    for lbl, col in zip(labels, colors):
        kmf_dict[lbl].plot(ax=ax_km, ci_show=False, color=col, linewidth=3)

    if style.legend_3g:
        ax_km.legend(loc="lower left", frameon=False, fontsize=16)
    else:
        ax_km.get_legend().remove()
    ax_km.set_xlabel("Months", fontsize=style.label_size_3g, fontweight="bold")
    ax_km.set_xticks(time_grid)
    ax_km.tick_params(axis="x", labelbottom=True, bottom=True, labelsize=style.tick_size_3g)
    if style.y_tick_size_3g:
        ax_km.tick_params(axis="y", labelsize=style.y_tick_size_3g)
    ax_km.minorticks_off()
    ax_km.set_title(title, fontsize=style.title_size_3g, fontweight="bold")
    ax_km.set_ylabel("Overall Survival Rate (%)", fontsize=style.label_size_3g, fontweight="bold")
    ax_km.spines["top"].set_visible(False)
    ax_km.spines["right"].set_visible(False)
    ax_km.set_xlim(0, max(time_grid))

    for i, lbl in enumerate(labels):
        y_frac = 0.75 - i * 0.15
        delta = (time_grid[1] - time_grid[0]) * 0.4
        x0 = time_grid[0] - delta
        x1 = x0 + delta * 0.5
        ax_table.plot([x0, x1], [y_frac, y_frac], color=colors[i], linewidth=3,
                      transform=ax_table.get_xaxis_transform(), solid_capstyle="butt", clip_on=False)
        for t, n in zip(time_grid, risk_dict[lbl][:len(time_grid)]):
            ax_table.text(t, y_frac, str(n), color=colors[i], fontsize=style.table_size_3g,
                          va="center", ha="center", transform=ax_table.get_xaxis_transform(),
                          clip_on=False)
    ax_table.text(-0.1, 0.9, "Number at Risk (n)", fontsize=12, ha="left", transform=ax_table.transAxes)
    ax_table.xaxis.set_visible(False)
    ax_table.yaxis.set_visible(False)
    plt.subplots_adjust(top=0.92, bottom=0.07)
    size, y = style.months_label_3g
    plt.text(0.5, y, "Months", fontsize=size, fontweight="bold", ha="center", transform=fig.transFigure)
    _save(file_name)
    plt.show()


def plot_km_by_category(df, category_col, title, llp, file_name=None):
    """Unadjusted KM curves for each level of a categorical column."""
    kmf = KaplanMeierFitter()
    plt.figure(figsize=(10, 6))
    muted_colors = ["#0008ff", "#ff0000", "#000000", "#008E28"]
    for i, group in enumerate(df[category_col].unique()):
        mask = df[category_col] == group
        kmf.fit(durations=df.loc[mask, "time"], event_observed=df.loc[mask, "event"],
                label=f"{group} (n={mask.sum()})")
        kmf.plot(ci_show=False, linewidth=2, color=muted_colors[i % len(muted_colors)])
    plt.text(0.15, 0.1, f"Log-likelihood p={llp:.4f}", fontsize=14, transform=plt.gca().transAxes,
             ha="center", va="center", bbox=dict(facecolor="white", alpha=0.5, edgecolor="black"))
    plt.title(title)
    plt.ylabel("Overall Survival Probability")
    plt.xlabel("Time (months)")
    plt.xlim(0, 72)
    plt.xticks(np.arange(0, 73, 12))
    plt.grid(False)
    plt.legend()
    plt.tight_layout()
    _save(file_name)
    plt.show()


def plot_iptw_diagnostics(df, groups, covariates=None, treatment_col="sequence_group", file_name=None):
    """Propensity score overlap before/after IPTW and a covariate balance (ASMD) plot.

    Returns the ASMD table.
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 5), sharey=True)
    for ax, weights, title in [
        (axes[0], None, "Propensity Score Distribution Pre-weighting"),
        (axes[1], "iptw_weight", "Propensity Score Distribution After IPTW"),
    ]:
        sns.histplot(data=df[df[treatment_col] == 1], x="ps", weights=weights, bins=25, stat="density",
                     color="skyblue", label=groups[1], element="step", ax=ax)
        sns.histplot(data=df[df[treatment_col] == 0], x="ps", weights=weights, bins=25, stat="density",
                     color="salmon", label=groups[0], element="step", ax=ax)
        ax.set_title(title)
        ax.set_xlabel("Propensity Score")
    axes[0].set_ylabel("Density")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, title="Treatment Group")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    _save(f"{file_name}_ps" if file_name else None)
    plt.show()

    asmd_df = asmd_table(df, covariates, treatment_col)
    plot_balance(asmd_df, "Covariate Balance: Unweighted vs IPTW",
                 file_name=f"{file_name}_iptw" if file_name else None)
    return asmd_df


def plot_balance(asmd_df, title, file_name=None):
    """Dot plot of ASMD per covariate for each adjustment method."""
    melted = asmd_df.melt(id_vars="Covariate", var_name="Method", value_name="ASMD")
    plt.figure(figsize=(10, 6))
    sns.pointplot(data=melted, x="ASMD", y="Covariate", hue="Method", dodge=True, linestyle="none",
                  markers=["o", "s", "D"], palette="Set2")
    plt.axvline(0.1, linestyle="--", color="orange", label="ASMD = 0.1 threshold")
    plt.xlabel("Absolute Standardized Mean Difference")
    plt.title(title)
    plt.legend(loc="upper right")
    plt.tight_layout()
    _save(file_name)
    plt.show()
