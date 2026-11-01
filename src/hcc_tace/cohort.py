"""Subsetting helpers for the analytic cohort."""

from . import config


def advanced_stage(df):
    """Restrict to AJCC stage III-IV."""
    return df[df["COMBINED_STAGE_GROUP"].isin(config.ADVANCED_STAGES)]


def pairwise_cohort(df, reference, comparator, col="sequence_group"):
    """Two-group subset recoded so `reference` -> 0 and `comparator` -> 1.

    Example: pairwise_cohort(df, config.SIMULTANEOUS, config.SYSTEMIC_FIRST)
    gives simultaneous (0) vs systemic first (1).
    """
    out = df[df[col].isin([reference, comparator])].copy()
    out[col] = out[col].map({reference: 0, comparator: 1})
    return out


def pairwise_cohorts(df):
    """The three pairwise comparisons used throughout the paper."""
    return {
        "tace_sys": pairwise_cohort(df, config.TACE_FIRST, config.SYSTEMIC_FIRST),
        "tace_same": pairwise_cohort(df, config.TACE_FIRST, config.SIMULTANEOUS),
        "same_sys": pairwise_cohort(df, config.SIMULTANEOUS, config.SYSTEMIC_FIRST),
    }
