"""Build the ranked gene list that feeds GSEA prerank.

GSEA is sensitive to how genes are ranked. For DESeq2 the Wald ``stat`` is the
recommended metric: it is signed (direction) and magnitude-aware, and behaves
well with ties. We also offer ``sign(log2FC) * -log10(pvalue)`` and raw
``log2FoldChange``.

When several rows map to the same key (e.g. after chicken->human ortholog
collapse), we keep the row with the largest absolute metric so the strongest
signal for that gene survives.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config


def compute_metric(df: pd.DataFrame, metric: str = config.RANK_METRIC) -> pd.Series:
    """Return the raw (unranked) ranking metric as a Series aligned to df."""
    if metric == "stat":
        if "stat" not in df.columns:
            raise KeyError("'stat' column required for the 'stat' ranking metric")
        return df["stat"].astype(float)
    if metric == "log2fc":
        return df["log2FoldChange"].astype(float)
    if metric == "signed_logp":
        p = df["pvalue"].astype(float).clip(lower=config.PVALUE_FLOOR)
        return np.sign(df["log2FoldChange"].astype(float)) * -np.log10(p)
    raise ValueError(f"Unknown ranking metric: {metric!r}")


def build_rank(
    df: pd.DataFrame,
    metric: str = config.RANK_METRIC,
    key_col: str = "gene_id",
) -> pd.Series:
    """Return a descending-sorted Series indexed by ``key_col``.

    Rows with NA metric or NA key are dropped. Duplicate keys are collapsed by
    keeping the value with the largest magnitude.
    """
    work = df[[key_col]].copy()
    work["_metric"] = compute_metric(df, metric).values
    work = work.dropna(subset=[key_col, "_metric"])
    work[key_col] = work[key_col].astype(str)

    # Collapse duplicate keys by max |metric|.
    work["_abs"] = work["_metric"].abs()
    work = (
        work.sort_values("_abs", ascending=False)
        .drop_duplicates(subset=[key_col], keep="first")
    )

    ranked = (
        work.set_index(key_col)["_metric"].sort_values(ascending=False)
    )
    ranked.index.name = key_col
    return ranked
