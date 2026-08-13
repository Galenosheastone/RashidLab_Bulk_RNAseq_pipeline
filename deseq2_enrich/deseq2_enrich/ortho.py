"""Chicken -> human ortholog mapping via g:Profiler's ``orth`` endpoint.

GSEA is run against human MSigDB collections, so the ranked chicken list is
relabelled to human symbols first. g:Profiler's orthology service is used
because it is the same backend as the ORA step (consistent identifiers) and
requires no local biomart download -- important on a 1 GB Streamlit box.

Mapping is keyed on the most reliable identifier available (see
``config.ORTHO_ID_PRIORITY``): stable Ensembl IDs first, symbols only as a
fallback. Chicken symbols are frequently LOC*/unannotated and map worse.

One-to-many orthologs are handled at the ranking level: when several chicken
genes map to the same human symbol, ``rank.build_rank`` keeps the strongest
signal (max |metric|). Here we simply return the tidy mapping, plus an
``OrthologReport`` describing how much of the list actually survived -- the
caller needs that to decide whether the GSEA is worth reporting at all.
"""
from __future__ import annotations

from dataclasses import dataclass
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


@lru_cache(maxsize=config.ORTHOLOG_CACHE_SIZE)
def _orth_cached(genes_key: tuple[str, ...], source: str, target: str) -> pd.DataFrame:
    gp = _client()
    res = gp.orth(
        organism=source,
        query=list(genes_key),
        target=target,
    )
    return res


def clear_cache() -> None:
    _orth_cached.cache_clear()


def cache_info():
    return _orth_cached.cache_info()


@dataclass
class OrthologReport:
    """What the ortholog step actually did, so it can be shown, not guessed.

    The pre-ranked GSEA null treats the ranked list as the universe, so a
    shrunken or biased list quietly distorts NES/FDR. Reporting the mapping
    rate alongside the results is the antidote.
    """

    id_col_used: str | None
    n_query: int              # distinct source IDs sent for mapping
    n_query_mapped: int       # source IDs that gained >=1 kept human ortholog
    n_human_symbols: int      # distinct human symbols after (optional) 1:1
    n_one_to_many: int        # source IDs -> >1 human symbol (before strict filter)
    n_many_to_one: int        # human symbols <- >1 source ID (before strict filter)
    strict_one_to_one: bool

    @property
    def mapping_rate(self) -> float:
        return self.n_query_mapped / self.n_query if self.n_query else 0.0

    def as_text(self) -> str:
        return (
            f"Ortholog mapping via '{self.id_col_used}': "
            f"{self.n_query_mapped}/{self.n_query} source IDs mapped "
            f"({self.mapping_rate:.0%}) -> {self.n_human_symbols} human symbols. "
            f"one->many: {self.n_one_to_many}, many->one: {self.n_many_to_one}"
            + (" [strict 1:1]" if self.strict_one_to_one else "")
        )


def map_to_human(
    genes: list[str],
    source: str = config.ORGANISM,
    target: str = config.ORTHOLOG_TARGET,
    strict_one_to_one: bool = False,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Map chicken gene IDs to human orthologs.

    Returns ``(mapping, multiplicity)`` where ``mapping`` is a tidy frame with
    ``incoming`` (query id) and ``ortholog_name`` (human symbol) columns, and
    ``multiplicity`` counts the ambiguous pairs *before* any strict filtering
    so the caller can report the cost of the trade-off.

    Rows with no ortholog (g:Profiler returns 'N/A') are dropped. With
    ``strict_one_to_one`` only unambiguous pairs survive, which stops one
    chicken measurement from appearing as several correlated human paralog
    ranks.
    """
    genes_key = tuple(dict.fromkeys(str(g) for g in genes))  # de-dup, keep order
    empty = pd.DataFrame(columns=["incoming", "ortholog_name"])
    no_mult = {"one_to_many": 0, "many_to_one": 0}
    if not genes_key:
        return empty, no_mult

    res = _orth_cached(genes_key, source, target)
    if res is None or len(res) == 0:
        return empty, no_mult

    # g:Profiler column names: 'incoming', 'ortholog_ensg', 'name'/'ortholog_name'
    name_col = "ortholog_name" if "ortholog_name" in res.columns else "name"
    keep = res[res[name_col].notna() & (res[name_col].astype(str) != "N/A")].copy()
    keep = keep.rename(columns={name_col: "ortholog_name"})
    keep = keep[["incoming", "ortholog_name"]].drop_duplicates()

    per_src = keep.groupby("incoming")["ortholog_name"].nunique()
    per_human = keep.groupby("ortholog_name")["incoming"].nunique()
    mult = {
        "one_to_many": int((per_src > 1).sum()),
        "many_to_one": int((per_human > 1).sum()),
    }

    if strict_one_to_one:
        good_src = per_src[per_src == 1].index
        good_human = per_human[per_human == 1].index
        keep = keep[
            keep["incoming"].isin(good_src) & keep["ortholog_name"].isin(good_human)
        ]
    return keep, mult


def attach_human_symbol(
    df: pd.DataFrame,
    id_col: str = "gene_id",
    source: str = config.ORGANISM,
    target: str = config.ORTHOLOG_TARGET,
    id_priority: tuple[str, ...] | None = None,
    strict_one_to_one: bool = False,
) -> tuple[pd.DataFrame, OrthologReport]:
    """Add a ``human_symbol`` column, trying identifier columns best-first.

    Columns are tried in ``config.ORTHO_ID_PRIORITY`` order -- stable Ensembl
    IDs before symbols -- because a real GRCg7b export often carries a symbol
    in ``gene_id`` and the reliable ENSGALG identifier in a separate column.
    Whichever column yields the most human symbols wins; the search
    short-circuits as soon as a column clears the guardrail floors, so the
    common case costs a single network round-trip.

    A caller-supplied ``id_col`` that is not part of the standard priority is
    honoured first, since that is an explicit choice rather than a default.

    Genes without an ortholog get NA and are naturally excluded from GSEA. A
    one-chicken-to-many-human relationship expands to multiple rows; the
    downstream ranking collapse then keeps the strongest per human symbol.

    Returns ``(df_with_human_symbol, OrthologReport)``.
    """
    if id_priority is None:
        priority = list(config.ORTHO_ID_PRIORITY)
        if id_col not in priority:
            priority.insert(0, id_col)
        id_priority = tuple(priority)

    candidate_cols = [c for c in id_priority if c in df.columns]
    if not candidate_cols:
        raise KeyError(f"None of {id_priority} present for ortholog mapping")

    best: tuple[pd.DataFrame, OrthologReport] | None = None
    for col in candidate_cols:
        merged, report = _attach_from_column(
            df, col, source, target, strict_one_to_one
        )
        if best is None or report.n_human_symbols > best[1].n_human_symbols:
            best = (merged, report)
        # Priority is best-first; once a column clears the floor, stop paying
        # for further orthology round-trips.
        if (
            report.mapping_rate >= config.ORTHO_MIN_MAPPING_RATE
            and report.n_human_symbols >= config.ORTHO_MIN_MAPPED_GENES
        ):
            break

    return best


def _attach_from_column(
    df: pd.DataFrame,
    id_col: str,
    source: str,
    target: str,
    strict_one_to_one: bool = False,
) -> tuple[pd.DataFrame, OrthologReport]:
    src_ids = df[id_col].dropna().astype(str).tolist()
    mapping, mult = map_to_human(src_ids, source, target, strict_one_to_one)

    merged = df.merge(
        mapping, left_on=id_col, right_on="incoming", how="left"
    ).rename(columns={"ortholog_name": "human_symbol"})
    if "incoming" in merged.columns:
        merged = merged.drop(columns=["incoming"])
    merged["human_symbol_source"] = id_col

    report = OrthologReport(
        id_col_used=id_col,
        n_query=int(df[id_col].dropna().astype(str).nunique()),
        n_query_mapped=int(mapping["incoming"].nunique()) if not mapping.empty else 0,
        n_human_symbols=int(merged["human_symbol"].dropna().nunique()),
        n_one_to_many=mult["one_to_many"],
        n_many_to_one=mult["many_to_one"],
        strict_one_to_one=strict_one_to_one,
    )
    return merged, report
