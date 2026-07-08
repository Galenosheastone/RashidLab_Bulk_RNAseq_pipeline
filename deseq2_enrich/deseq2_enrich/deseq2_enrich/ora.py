"""Over-representation analysis (ORA) via g:Profiler.

g:Profiler natively supports *Gallus gallus*, so chicken GO/KEGG/Reactome/
WikiPathways terms are tested directly with no ortholog loss. Crucially we pass
the tested-gene set as a **custom background** (``domain_scope='custom'``) so
enrichment is measured against genes that could have been detected, not the
whole genome.

Each direction (up / down / all) is run separately and tagged.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd

from . import config

try:
    from gprofiler import GProfiler
except Exception:  # pragma: no cover
    GProfiler = None


def _client() -> "GProfiler":
    if GProfiler is None:
        raise RuntimeError(
            "gprofiler-official is not installed. `pip install gprofiler-official`."
        )
    return GProfiler(return_dataframe=True, user_agent="deseq2_enrich")


@lru_cache(maxsize=16)
def _profile_cached(
    query_key: tuple[str, ...],
    background_key: tuple[str, ...],
    sources_key: tuple[str, ...],
    organism: str,
) -> pd.DataFrame:
    gp = _client()
    kwargs = dict(
        organism=organism,
        query=list(query_key),
        sources=list(sources_key),
        user_threshold=config.SIG_ALPHA,
        significance_threshold_method="g_SCS",
        no_evidences=False,
    )
    if background_key:
        kwargs["background"] = list(background_key)
        kwargs["domain_scope"] = "custom"
    return gp.profile(**kwargs)


def run_ora(
    query_genes: list[str],
    background_genes: list[str] | None,
    sources: list[str] | None = None,
    organism: str = config.ORGANISM,
    direction_label: str = "all",
) -> pd.DataFrame:
    """Run g:Profiler ORA for one gene list.

    Returns a tidy dataframe (one row per enriched term) with a ``direction``
    column and a derived ``neg_log10_p`` and ``gene_ratio``.
    """
    if not query_genes:
        return _empty_ora()
    sources = sources or config.ORA_DEFAULT_SOURCES
    res = _profile_cached(
        tuple(query_genes),
        tuple(background_genes) if background_genes else tuple(),
        tuple(sources),
        organism,
    )
    if res is None or len(res) == 0:
        return _empty_ora()

    res = res.copy()
    # Standardise the columns we rely on downstream. Keep this as a copy
    # instead of a bare DataFrame.rename so already-normalised inputs work too.
    canonical = {
        "native": "term_id",
        "name": "term_name",
        "p_value": "p_value",
        "term_size": "term_size",
        "query_size": "query_size",
        "intersection_size": "intersection_size",
        "source": "source",
    }
    for observed, standard in canonical.items():
        if standard in res.columns:
            continue
        if observed in res.columns:
            res[standard] = res[observed]
        else:
            res[standard] = np.nan
    res["neg_log10_p"] = -np.log10(res["p_value"].clip(lower=config.PVALUE_FLOOR))
    # gene_ratio = intersection / query ; recall = intersection / term_size
    res["gene_ratio"] = res["intersection_size"] / res["query_size"].replace(0, np.nan)
    res["recall"] = res["intersection_size"] / res["term_size"].replace(0, np.nan)
    res["direction"] = direction_label
    keep = [
        "source", "term_id", "term_name", "p_value", "neg_log10_p",
        "term_size", "query_size", "intersection_size",
        "gene_ratio", "recall", "direction",
    ]
    if "intersections" in res.columns:
        # list of query genes hitting the term (present when no_evidences=False)
        res["genes"] = res["intersections"].apply(_flatten_intersection)
        keep.append("genes")
    return res[[c for c in keep if c in res.columns]].sort_values("p_value")


def run_ora_directional(
    deg_sets,
    sources: list[str] | None = None,
    organism: str = config.ORGANISM,
    directions: tuple[str, ...] = ("up", "down"),
) -> pd.DataFrame:
    """Run ORA for each requested direction against the shared universe."""
    frames = []
    for direction in directions:
        query = getattr(deg_sets, "all_sig" if direction == "all" else direction)
        frames.append(
            run_ora(query, deg_sets.universe, sources, organism, direction)
        )
    frames = [frame for frame in frames if frame is not None and not frame.empty]
    if not frames:
        return _empty_ora()
    return pd.concat(frames, ignore_index=True)


def _flatten_intersection(x):
    """Return flat gene labels from g:Profiler's intersections column.

    g:Profiler's intersections column is a flat list of query gene labels that
    hit the term. Nested evidence codes live in a separate column, not here.
    """
    if x is None:
        return []
    if isinstance(x, (list, tuple)):
        return [str(item) for item in x if item is not None]
    return []


def _empty_ora() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "source", "term_id", "term_name", "p_value", "neg_log10_p",
            "term_size", "query_size", "intersection_size",
            "gene_ratio", "recall", "direction", "genes",
        ]
    )
