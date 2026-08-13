"""Pre-ranked GSEA via gseapy.

The ranked list (human symbols after ortholog mapping, or a user-supplied
ranking) is scored against gene-set collections. We keep the gseapy
``Prerank`` result object for term-level plots, but trim large running-curve
arrays for lower-priority terms so broad GO libraries do not dominate memory.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field

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
    # 'auto' mode runs two independent preranks (native chicken + human
    # ortholog) whose rankings live in different ID spaces and have different
    # lengths. They are kept separate here rather than pooled, because a
    # running curve is only meaningful against the ranking it was computed
    # from -- and because FDR must not be pooled across two nulls.
    routes: dict[str, "GSEAResult"] = field(default_factory=dict)


def result_for_term(result: GSEAResult, term: str) -> GSEAResult:
    """Return the per-route result that actually scored ``term``.

    Single-route runs return themselves, so callers need no special case.
    """
    for route in (result.routes or {}).values():
        raw_results = getattr(getattr(route, "raw", None), "results", {}) or {}
        if term in raw_results:
            return route
    return result


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
    # gseapy matches genes by exact string, so a single case mismatch between
    # the ranking and the gene sets silently zeroes the overlap rather than
    # erroring. Normalise BOTH sides to upper case. Only dict gene sets can be
    # transformed -- a GMT path or an Enrichr library name is a string gseapy
    # resolves itself and must be passed through untouched. Human HGNC symbols,
    # g:Profiler ortholog names and ENSGALG IDs are already upper case, so this
    # is a no-op for them and a rescue for anything that is not.
    rnk = ranking.copy()
    rnk.index = rnk.index.astype(str).str.upper()
    rnk = rnk[~rnk.index.duplicated(keep="first")]
    rnk = rnk.sort_values(ascending=False)
    if isinstance(gene_sets, dict):
        # Re-intern after upper-casing: the native chicken GMT is ~20k sets over
        # ~16k distinct genes with heavy overlap, and without interning the new
        # upper-cased strings the dict balloons on a memory-limited host.
        gene_sets = {
            term: list({sys.intern(str(g).upper()) for g in members})
            for term, members in gene_sets.items()
        }
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
    tidy = _tidy(pre.res2d, n_perm=permutations)
    running_terms = _compact_raw(pre, tidy, max_running_terms, len(rnk))
    return GSEAResult(table=tidy, raw=pre, ranking=rnk, running_terms=running_terms)


# Private gseapy attributes holding the bulky inputs we drop after scoring to
# stay inside a ~1 GB host. Centralised so a gseapy upgrade that renames them
# is a one-line change -- and so test_run_prerank_real_smoke, which exercises
# the real object, fails loudly in CI rather than silently in production.
#
# Verified against gseapy 1.1.12 and 1.3.1 (2026-08-13): the full suite passes
# on both, and a real 20k-set chicken prerank returns identical scores
# (3756 sets, top NES 2.842598) with the running plot and leading-edge heatmap
# rendering from pre.results afterwards. 1.3.1 is what the deployed
# requirements.txt resolves to.
_GSEAPY_HEAVY_ATTRS = ("gene_sets", "gmt", "_gmt", "ranking", "rnk")


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

    for attr in _GSEAPY_HEAVY_ATTRS:
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


def _tidy(res2d: pd.DataFrame, n_perm: int = config.GSEA_PERMUTATIONS) -> pd.DataFrame:
    """Normalise gseapy's res2d into consistent, sortable columns.

    ``n_perm`` sets the display floor for p/FDR: a permutation test cannot
    resolve a p-value below ``1/(n_perm + 1)``, so gseapy printing 0.0 is an
    artifact of the resolution, not evidence of p == 0.
    """
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
    # Floor the displayed significance at the permutation resolution so a
    # reported "0" cannot be mistaken for an exact zero.
    floor = 1.0 / (int(n_perm) + 1) if n_perm else 0.0
    if floor:
        for col in ("pval", "fdr", "fwer"):
            if col in df.columns:
                df[col] = df[col].mask(df[col] < floor, floor)

    # np.where(NES >= 0, "up", "down") labelled a NaN NES as "down", because
    # any comparison with NaN is False. gseapy returns NaN NES for sets whose
    # null collapsed, so those were being reported as down-regulated. Select
    # explicitly and give the undecidable case its own label.
    nes = pd.to_numeric(df.get("NES"), errors="coerce") if "NES" in df.columns \
        else pd.Series(np.nan, index=df.index)
    df["direction"] = np.select([nes > 0, nes < 0], ["up", "down"], default="ns")
    sort_key = "fdr" if "fdr" in df.columns else "pval"
    return df.sort_values(sort_key).reset_index(drop=True)
