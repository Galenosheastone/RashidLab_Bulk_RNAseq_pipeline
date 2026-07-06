"""Differentially-expressed-gene selection and the enrichment universe.

Two decisions are baked in here because they are the two most common ORA
mistakes:

1. The background/universe is the set of *tested* genes (non-NA padj), not the
   whole genome. Using the genome inflates enrichment.
2. Up- and down-regulated genes are analysed separately by default. A merged
   list mixes opposing biology and dilutes signal; directional ORA reads far
   more cleanly.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from . import config


@dataclass
class DEGSets:
    up: list[str]
    down: list[str]
    all_sig: list[str]
    universe: list[str]
    id_col: str
    padj_threshold: float
    lfc_threshold: float

    @property
    def counts(self) -> dict[str, int]:
        return {
            "up": len(self.up),
            "down": len(self.down),
            "all_sig": len(self.all_sig),
            "universe": len(self.universe),
        }


def select_degs(
    df: pd.DataFrame,
    padj_threshold: float = config.PADJ_THRESHOLD,
    lfc_threshold: float = config.LFC_THRESHOLD,
    id_col: str = "gene_id",
) -> DEGSets:
    """Split a DESeq2 table into up / down / significant / universe gene lists.

    The universe is every gene with a non-NA padj (i.e. actually tested).
    Genes are de-duplicated while preserving order.
    """
    if id_col not in df.columns:
        raise KeyError(f"id_col '{id_col}' not in dataframe")

    tested = df[df["padj"].notna()].copy()
    universe = _unique(tested[id_col])

    sig = tested[
        (tested["padj"] < padj_threshold)
        & (tested["log2FoldChange"].abs() > lfc_threshold)
    ]
    up = _unique(sig.loc[sig["log2FoldChange"] > 0, id_col])
    down = _unique(sig.loc[sig["log2FoldChange"] < 0, id_col])
    all_sig = _unique(sig[id_col])

    return DEGSets(
        up=up,
        down=down,
        all_sig=all_sig,
        universe=universe,
        id_col=id_col,
        padj_threshold=padj_threshold,
        lfc_threshold=lfc_threshold,
    )


def _unique(series: pd.Series) -> list[str]:
    """Drop NA, cast to str, de-duplicate preserving order."""
    seen: set[str] = set()
    out: list[str] = []
    for v in series.dropna().astype(str):
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out
