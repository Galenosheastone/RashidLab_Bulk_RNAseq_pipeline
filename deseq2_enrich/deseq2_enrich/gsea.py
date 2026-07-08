"""Pre-ranked GSEA via gseapy.

The ranked list (human symbols after ortholog mapping, or a user-supplied
ranking) is scored against gene-set collections. We keep the full gseapy
``Prerank`` result object so the running-enrichment-score curves and
leading-edge subsets can be plotted without recomputation.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import gseapy as gp

from . import config


@dataclass
class GSEAResult:
    table: pd.DataFrame          # tidy results
    raw: object                  # gseapy Prerank object (for running plots)
    ranking: pd.Series           # the ranked list actually used


def run_prerank(
    ranking: pd.Series,
    gene_sets,
    min_size: int = config.GSEA_MIN_SIZE,
    max_size: int = config.GSEA_MAX_SIZE,
    permutations: int = config.GSEA_PERMUTATIONS,
    seed: int = config.GSEA_SEED,
    threads: int = 1,
) -> GSEAResult:
    """Run gseapy prerank and return tidy + raw results.

    Parameters
    ----------
    ranking : pd.Series
        Descending-sorted, indexed by gene symbol. Duplicate indices must
        already be collapsed (see ``rank.build_rank``).
    gene_sets : dict | str
        ``{term: [genes]}`` dict, a GMT path, or an Enrichr library name.
    """
    rnk = ranking.copy()
    rnk.index = rnk.index.astype(str)
    rnk = rnk[~rnk.index.duplicated(keep="first")]
    rnk = rnk.sort_values(ascending=False)

    pre = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        min_size=min_size,
        max_size=max_size,
        permutation_num=permutations,
        seed=seed,
        threads=threads,
        outdir=None,           # in-memory only; nothing written to disk
        no_plot=True,
        verbose=False,
    )
    tidy = _tidy(pre.res2d)
    return GSEAResult(table=tidy, raw=pre, ranking=rnk)


def _tidy(res2d: pd.DataFrame) -> pd.DataFrame:
    """Normalise gseapy's res2d into consistent, sortable columns."""
    df = res2d.copy()
    # gseapy column names differ slightly across versions; normalise.
    colmap = {
        "Term": "term",
        "ES": "ES",
        "NES": "NES",
        "NOM p-val": "pval",
        "FDR q-val": "fdr",
        "FWER p-val": "fwer",
        "Gene %": "gene_pct",
        "Tag %": "tag_pct",
        "Lead_genes": "lead_genes",
        "Genes": "genes",
    }
    df = df.rename(columns={k: v for k, v in colmap.items() if k in df.columns})
    for numcol in ("ES", "NES", "pval", "fdr", "fwer"):
        if numcol in df.columns:
            df[numcol] = pd.to_numeric(df[numcol], errors="coerce")
    if "term" in df.columns:
        # split "TAG | term" provenance if present
        split = df["term"].astype(str).str.split(r"\s*\|\s*", n=1, expand=True)
        if split.shape[1] == 2:
            df["collection"] = split[0]
            df["term_short"] = split[1]
        else:
            df["collection"] = "custom"
            df["term_short"] = df["term"]
    df["direction"] = np.where(df.get("NES", 0) >= 0, "up", "down")
    sort_key = "fdr" if "fdr" in df.columns else "pval"
    return df.sort_values(sort_key).reset_index(drop=True)
