"""Pre-ranked GSEA via gseapy.

The ranked list (human symbols after ortholog mapping, or a user-supplied
ranking) is scored against gene-set collections. We keep the gseapy
``Prerank`` result object for term-level plots, but trim large running-curve
arrays for lower-priority terms so broad GO libraries do not dominate memory.
"""
from __future__ import annotations

import gc
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


def _normalise_sets(gene_sets: dict) -> dict:
    """Upper-case and re-intern one block of gene sets.

    Done per block rather than over the whole dict: the transient copy of a
    ~4,400-set collection is itself a meaningful share of the memory peak.
    """
    return {
        term: list({sys.intern(str(g).upper()) for g in members})
        for term, members in gene_sets.items()
    }


def _global_bh_fdr(pvals: pd.Series) -> pd.Series:
    """Benjamini-Hochberg across the pooled nominal p-values.

    Per-chunk q-values are computed only within their block, so they are both
    noisier and dependent on an arbitrary partition. Pooling the nominal p and
    adjusting once is chunk-stable and is the same choice fgsea's ``padj``
    makes.

    NOTE (documented trade-off): this is *not* identical to gseapy's
    single-run q-value, which is NES-based and derived from the pooled null
    distribution rather than from p alone. It is monotonic in nominal p and
    defensible as "BH-adjusted permutation p-values", but a run reported this
    way must say so. Reproducing gseapy's exact q under chunking would require
    pooling each set's null-NES distribution across chunks -- out of scope.
    """
    from statsmodels.stats.multitest import multipletests

    p = pd.to_numeric(pvals, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    ok = p.notna()
    if ok.any():
        out.loc[ok] = multipletests(p[ok].to_numpy(), method="fdr_bh")[1]
    return out


def run_prerank_chunked(
    ranking: pd.Series,
    gene_sets: dict,
    min_size: int = config.GSEA_MIN_SIZE,
    max_size: int = config.GSEA_MAX_SIZE,
    permutations: int = config.GSEA_PERMUTATIONS,
    seed: int = config.GSEA_SEED,
    threads: int = 1,
    max_running_terms: int | None = config.GSEA_MAX_RUNNING_PLOT_TERMS,
    chunk_size: int = config.GSEA_CHUNK_SIZE,
) -> GSEAResult:
    """Score a large collection in blocks so peak memory stays bounded.

    Two phases:

    A. Score ``gene_sets`` in blocks of ``chunk_size``, keeping only the tidy
       table from each block and freeing the gseapy object immediately. Peak is
       bounded by one block instead of the whole collection. ES/NES/nominal p
       are unaffected -- gseapy seeds its permutations per gene set, so a set
       scores identically regardless of what it was scored alongside.
    B. Re-score just the top ``max_running_terms`` sets in a single small
       prerank to obtain a real gseapy object with running curves, so the
       running-sum plots and ``result_for_term`` keep working unchanged.

    FDR is recomputed globally over the pooled nominal p (see
    ``_global_bh_fdr``); per-block q-values are discarded.
    """
    rnk = ranking.copy()
    rnk.index = rnk.index.astype(str).str.upper()
    rnk = rnk[~rnk.index.duplicated(keep="first")]
    rnk = rnk.sort_values(ascending=False)
    if len(rnk) < 2:
        raise ValueError(
            "GSEA requires at least 2 ranked genes after ID/ortholog mapping; "
            f"got {len(rnk)}."
        )

    terms = list(gene_sets)
    frames: list[pd.DataFrame] = []
    for start in range(0, len(terms), max(1, chunk_size)):
        block = {t: gene_sets[t] for t in terms[start:start + chunk_size]}
        pre = gp.prerank(
            rnk=rnk,
            gene_sets=_normalise_sets(block),
            min_size=min_size,
            max_size=max_size,
            permutation_num=permutations,
            seed=seed,
            threads=threads,
            outdir=None,
            no_plot=True,
            verbose=False,
        )
        res2d = getattr(pre, "res2d", None)
        if res2d is not None and len(res2d):
            frames.append(res2d.copy())
        # Drop the block's running curves and the gseapy object before the
        # next block allocates: this is what keeps the peak to one block.
        del pre, res2d, block
        gc.collect()

    if not frames:
        raise ValueError(
            "GSEA scored no gene sets: every set fell outside "
            f"[{min_size}, {max_size}] after intersecting with the ranked list."
        )

    pooled = pd.concat(frames, ignore_index=True)
    del frames
    gc.collect()

    table = _tidy(pooled, n_perm=permutations, apply_floor=False)

    # Floor the nominal p BEFORE adjusting, not after. gseapy reports p == 0
    # whenever no permutation beat the observed ES, and BH maps 0 -> 0; applying
    # the display floor afterwards would then print q = 1/(nperm+1) for the top
    # terms, which reads as "survives correction across thousands of sets" when
    # the run simply lacks the resolution to say anything of the sort. Measured
    # on a 120-set run at 100 permutations: adjust-then-floor gives q = 0.0099
    # for the top term, floor-then-adjust gives q = 0.50 -- a 50x difference,
    # entirely an artifact of the ordering. max(p, 1/(nperm+1)) is the standard
    # add-one permutation estimator, so this is also the conventional choice.
    table = _apply_display_floor(table, permutations)
    table["fdr"] = _global_bh_fdr(table["pval"]).to_numpy()
    table["fdr_method"] = "BH pooled (chunked run)"
    # Idempotent for pval; bounds fdr into the same display range.
    table = _apply_display_floor(table, permutations)
    sort_key = "fdr" if "fdr" in table.columns else "pval"
    table = table.sort_values(sort_key).reset_index(drop=True)

    # --- Phase B: a small real prerank purely for the running curves -------
    raw, running = None, ()
    if max_running_terms and max_running_terms > 0:
        ranked = table.assign(
            _fdr=pd.to_numeric(table["fdr"], errors="coerce").fillna(np.inf),
            _abs=pd.to_numeric(table.get("NES"), errors="coerce").abs().fillna(0.0),
        ).sort_values(["_fdr", "_abs"], ascending=[True, False])
        top = [t for t in ranked["term"].astype(str).tolist() if t in gene_sets]
        top = top[:max_running_terms]
        if top:
            sub = run_prerank(
                rnk, {t: gene_sets[t] for t in top},
                min_size=min_size, max_size=max_size,
                permutations=permutations, seed=seed, threads=threads,
                max_running_terms=max_running_terms,
            )
            raw, running = sub.raw, sub.running_terms

    return GSEAResult(table=table, raw=raw, ranking=rnk, running_terms=running)


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


def _apply_display_floor(df: pd.DataFrame, n_perm: int) -> pd.DataFrame:
    """Floor p/FDR/FWER at the permutation resolution, ``1/(n_perm + 1)``.

    A permutation test cannot resolve below that, so a reported 0.0 is an
    artifact of the resolution rather than evidence of p == 0.
    """
    floor = 1.0 / (int(n_perm) + 1) if n_perm else 0.0
    if floor:
        for col in ("pval", "fdr", "fwer"):
            if col in df.columns:
                df[col] = df[col].mask(df[col] < floor, floor)
    return df


def _tidy(res2d: pd.DataFrame, n_perm: int = config.GSEA_PERMUTATIONS,
          apply_floor: bool = True) -> pd.DataFrame:
    """Normalise gseapy's res2d into consistent, sortable columns.

    ``apply_floor=False`` leaves raw p-values in place so a caller can pool
    results across chunks and recompute FDR globally before flooring; the
    floor is a display transform and must be the last step.
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
    if apply_floor:
        df = _apply_display_floor(df, n_perm)

    # np.where(NES >= 0, "up", "down") labelled a NaN NES as "down", because
    # any comparison with NaN is False. gseapy returns NaN NES for sets whose
    # null collapsed, so those were being reported as down-regulated. Select
    # explicitly and give the undecidable case its own label.
    nes = pd.to_numeric(df.get("NES"), errors="coerce") if "NES" in df.columns \
        else pd.Series(np.nan, index=df.index)
    df["direction"] = np.select([nes > 0, nes < 0], ["up", "down"], default="ns")
    sort_key = "fdr" if "fdr" in df.columns else "pval"
    return df.sort_values(sort_key).reset_index(drop=True)
