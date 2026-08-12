"""Pre-ranked GSEA via gseapy.

The ranked list (human symbols after ortholog mapping, or a user-supplied
ranking) is scored against gene-set collections. We keep the gseapy
``Prerank`` result object for term-level plots, but trim large running-curve
arrays for lower-priority terms so broad GO libraries do not dominate memory.
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
    running_terms: tuple[str, ...] = ()


def run_prerank(
    ranking: pd.Series,
    gene_sets,
    min_size: int = config.GSEA_MIN_SIZE,
    max_size: int = config.GSEA_MAX_SIZE,
    permutations: int = config.GSEA_PERMUTATIONS,
    seed: int = config.GSEA_SEED,
    threads: int = 1,
    max_running_terms: int | None = config.GSEA_MAX_RUNNING_PLOT_TERMS,
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
    if len(rnk) < 2:
        raise ValueError(
            "GSEA requires at least 2 ranked genes after ID/ortholog mapping; "
            f"got {len(rnk)}. Check the selected ID column or provide gene symbols."
        )

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
    running_terms = _compact_raw(pre, tidy, max_running_terms, len(rnk))
    return GSEAResult(table=tidy, raw=pre, ranking=rnk, running_terms=running_terms)


def _compact_raw(pre, tidy: pd.DataFrame,
                 max_running_terms: int | None,
                 ranking_size: int) -> tuple[str, ...]:
    """Drop heavy running-plot arrays except for the strongest terms."""
    results = getattr(pre, "results", {}) or {}
    if tidy.empty or "term" not in tidy.columns or not isinstance(results, dict):
        return tuple()

    available = [
        term for term in tidy["term"].astype(str).tolist()
        if term in results and _has_running_payload(results[term])
    ]
    if max_running_terms is None or max_running_terms < 1:
        keep = set()
    elif len(available) <= max_running_terms:
        keep = set(available)
    else:
        ranked = tidy[tidy["term"].astype(str).isin(available)].copy()
        fdr_values = (
            ranked["fdr"] if "fdr" in ranked.columns
            else pd.Series(np.nan, index=ranked.index)
        )
        nes_values = (
            ranked["NES"] if "NES" in ranked.columns
            else pd.Series(0.0, index=ranked.index)
        )
        ranked["_fdr_sort"] = pd.to_numeric(
            fdr_values, errors="coerce"
        ).fillna(np.inf)
        ranked["_abs_nes"] = pd.to_numeric(
            nes_values, errors="coerce"
        ).abs().fillna(0.0)
        keep = set(
            ranked.sort_values(
                ["_fdr_sort", "_abs_nes"], ascending=[True, False]
            )["term"].astype(str).head(max_running_terms)
        )

    for term, payload in results.items():
        if term in keep or not isinstance(payload, dict):
            continue
        payload["RES"] = _PrunedRunningCurve(ranking_size)
        payload["hits"] = []
        payload.pop("RES_null", None)

    for attr in ("gene_sets", "gmt", "_gmt", "ranking", "rnk"):
        if hasattr(pre, attr):
            try:
                setattr(pre, attr, None)
            except Exception:
                pass

    return tuple(term for term in available if term in keep)


def _has_running_payload(payload) -> bool:
    return isinstance(payload, dict) and "RES" in payload and "hits" in payload


class _PrunedRunningCurve:
    """Tiny stand-in so older UI code does not fail on pruned terms."""

    def __init__(self, length: int):
        self.length = int(length)

    def __len__(self) -> int:
        return self.length

    def __array__(self, dtype=None):
        return np.zeros(self.length, dtype=dtype or float)


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
