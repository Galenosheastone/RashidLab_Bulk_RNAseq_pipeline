"""Chicken -> human ortholog mapping via g:Profiler's ``orth`` endpoint.

GSEA is run against human MSigDB collections, so the ranked chicken list is
relabelled to human symbols first. g:Profiler's orthology service is used
because it is the same backend as the ORA step (consistent identifiers) and
requires no local biomart download -- important on a 1 GB Streamlit box.

One-to-many orthologs are handled at the ranking level: when several chicken
genes map to the same human symbol, ``rank.build_rank`` keeps the strongest
signal (max |metric|). Here we simply return the tidy mapping.
"""
from __future__ import annotations

from functools import lru_cache

import pandas as pd

from . import config

try:
    from gprofiler import GProfiler
except Exception:  # pragma: no cover - import guarded for offline test envs
    GProfiler = None


def _client() -> "GProfiler":
    if GProfiler is None:
        raise RuntimeError(
            "gprofiler-official is not installed. `pip install gprofiler-official`."
        )
    return GProfiler(return_dataframe=True, user_agent="deseq2_enrich")


@lru_cache(maxsize=8)
def _orth_cached(genes_key: tuple[str, ...], source: str, target: str) -> pd.DataFrame:
    gp = _client()
    res = gp.orth(
        organism=source,
        query=list(genes_key),
        target=target,
    )
    return res


def map_to_human(
    genes: list[str],
    source: str = config.ORGANISM,
    target: str = config.ORTHOLOG_TARGET,
) -> pd.DataFrame:
    """Map chicken gene IDs to human orthologs.

    Returns a tidy dataframe with at least ``incoming`` (query id) and
    ``ortholog_name`` (human symbol) columns. Rows with no ortholog
    (g:Profiler returns 'N/A') are dropped.
    """
    genes_key = tuple(dict.fromkeys(str(g) for g in genes))  # de-dup, keep order
    if not genes_key:
        return pd.DataFrame(columns=["incoming", "ortholog_name"])

    res = _orth_cached(genes_key, source, target)
    if res is None or len(res) == 0:
        return pd.DataFrame(columns=["incoming", "ortholog_name"])

    # g:Profiler column names: 'incoming', 'ortholog_ensg', 'name'/'ortholog_name'
    name_col = "ortholog_name" if "ortholog_name" in res.columns else "name"
    keep = res[res[name_col].notna() & (res[name_col] != "N/A")].copy()
    keep = keep.rename(columns={name_col: "ortholog_name"})
    return keep[["incoming", "ortholog_name"]].drop_duplicates()


def attach_human_symbol(
    df: pd.DataFrame,
    id_col: str = "gene_id",
    source: str = config.ORGANISM,
    target: str = config.ORTHOLOG_TARGET,
    fallback_cols: tuple[str, ...] = ("gene_name",),
) -> pd.DataFrame:
    """Return a copy of ``df`` with a ``human_symbol`` column added.

    Genes without an ortholog get NA and are naturally excluded from GSEA.
    A one-chicken-to-many-human relationship expands to multiple rows; the
    downstream ranking collapse then keeps the strongest per human symbol.
    If the primary ID column cannot be mapped well enough to run GSEA, symbol
    columns such as ``gene_name`` are tried as a fallback. This catches DESeq2
    exports where ``gene_id`` contains NCBI numeric IDs that g:Profiler's
    orthology endpoint does not resolve.
    """
    if id_col not in df.columns:
        raise KeyError(f"{id_col!r} column required for ortholog mapping")

    candidate_cols = [id_col]
    candidate_cols.extend(
        c for c in fallback_cols if c in df.columns and c not in candidate_cols
    )

    out = None
    for col in candidate_cols:
        mapped = _attach_from_column(df, col, source, target)
        if mapped is not None and int(mapped["human_symbol"].notna().sum()) > 1:
            out = mapped
            break

    if out is not None:
        return out

    out = df.copy()
    out["human_symbol"] = pd.NA
    out["human_symbol_source"] = pd.NA
    return out


def _attach_from_column(
    df: pd.DataFrame,
    id_col: str,
    source: str,
    target: str,
) -> pd.DataFrame | None:
    mapping = map_to_human(df[id_col].dropna().astype(str).tolist(), source, target)
    if mapping.empty:
        return None

    merged = df.merge(mapping, left_on=id_col, right_on="incoming", how="left")
    merged = merged.rename(columns={"ortholog_name": "human_symbol"})
    if "incoming" in merged.columns:
        merged = merged.drop(columns=["incoming"])
    merged["human_symbol_source"] = id_col
    return merged
