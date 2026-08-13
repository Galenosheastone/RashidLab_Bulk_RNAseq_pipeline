from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from deseq2_enrich import gsea


def test_run_prerank_rejects_empty_ranking_before_gseapy_assertion():
    with pytest.raises(ValueError, match="at least 2 ranked genes"):
        gsea.run_prerank(
            pd.Series(dtype=float),
            {"TERM": ["A", "B", "C"]},
            permutations=10,
        )


def test_run_prerank_trims_running_arrays_for_lower_priority_terms():
    fake_pre = SimpleNamespace(
        res2d=pd.DataFrame({
            "Term": ["TERM_A", "TERM_B"],
            "ES": [0.6, 0.5],
            "NES": [2.1, 2.5],
            "NOM p-val": [0.01, 0.02],
            "FDR q-val": [0.01, 0.20],
            "FWER p-val": [0.01, 0.20],
            "Lead_genes": ["A;B", "C;D"],
        }),
        results={
            "TERM_A": {
                "nes": 2.1,
                "fdr": 0.01,
                "RES": [0.0, 0.4, 0.6],
                "hits": [1, 2],
                "lead_genes": "A;B",
            },
            "TERM_B": {
                "nes": 2.5,
                "fdr": 0.20,
                "RES": [0.0, 0.3, 0.5],
                "hits": [0, 2],
                "lead_genes": "C;D",
            },
        },
        gene_sets={"large": ["A", "B", "C"]},
    )

    ranking = pd.Series([2.0, 1.0, -1.0], index=["A", "B", "C"])
    with patch.object(gsea.gp, "prerank", return_value=fake_pre):
        result = gsea.run_prerank(
            ranking,
            {"TERM_A": ["A", "B"], "TERM_B": ["C", "D"]},
            permutations=10,
            max_running_terms=1,
        )

    assert result.running_terms == ("TERM_A",)
    assert "RES" in result.raw.results["TERM_A"]
    assert "hits" in result.raw.results["TERM_A"]
    assert len(result.raw.results["TERM_B"]["RES"]) == len(ranking)
    assert result.raw.results["TERM_B"]["hits"] == []
    assert result.raw.gene_sets is None


# --------------------------------------------------------------------------
# Phase 2: case normalisation, p-value floor, real prerank smoke
# --------------------------------------------------------------------------
def test_case_mismatch_does_not_zero_the_overlap():
    """Lower-case ranking vs upper-case gene sets must still score."""
    rng = np.random.default_rng(7)
    genes = [f"g{i:03d}" for i in range(60)]          # lower case ranking
    scores = pd.Series(rng.normal(size=60), index=genes).sort_values(ascending=False)
    gene_sets = {"TOP": [g.upper() for g in genes[:20]]}  # upper case members

    res = gsea.run_prerank(scores, gene_sets, min_size=5, max_size=50,
                           permutations=50, seed=1, threads=1)
    assert not res.table.empty, "case mismatch silently zeroed the overlap"


def test_pval_floored_at_permutation_resolution():
    """A permutation test cannot resolve below 1/(n+1); 0.0 must not survive."""
    res2d = pd.DataFrame({
        "Term": ["A | x", "B | y"],
        "NES": [2.0, -1.5],
        "NOM p-val": [0.0, 0.5],
        "FDR q-val": [0.0, 0.4],
        "FWER p-val": [0.0, 0.6],
    })
    tidy = gsea._tidy(res2d, n_perm=100)
    floor = 1.0 / 101
    assert tidy["pval"].min() == pytest.approx(floor)
    assert tidy["fdr"].min() == pytest.approx(floor)
    assert tidy["fwer"].min() == pytest.approx(floor)
    # Values above the floor are untouched.
    assert tidy.loc[tidy["term"] == "B | y", "pval"].iloc[0] == pytest.approx(0.5)


def test_run_prerank_real_smoke():
    """A genuine prerank: guards the gseapy internals the plots depend on.

    This is the gate on the ``gseapy<2`` cap. It must pass on the version the
    deployment actually installs before that cap is touched. Verified on
    1.1.12 and 1.3.1; run it under the target version before any bump.
    """
    rng = np.random.default_rng(0)
    genes = [f"G{i:03d}" for i in range(60)]
    scores = pd.Series(rng.normal(size=60), index=genes).sort_values(ascending=False)
    # a term enriched at the top, a term enriched at the bottom
    gene_sets = {"TOP": list(scores.index[:20]), "BOTTOM": list(scores.index[-20:])}

    res = gsea.run_prerank(scores, gene_sets, min_size=5, max_size=50,
                           permutations=100, seed=1, threads=1)
    assert not res.table.empty
    assert {"term", "NES", "fdr"}.issubset(res.table.columns)
    # internals the app/plots depend on still exist after _compact_raw
    assert isinstance(res.raw.results, dict)
    # the heavy inputs really were dropped
    for attr in gsea._GSEAPY_HEAVY_ATTRS:
        assert getattr(res.raw, attr, None) is None
    # TOP should score positive NES, BOTTOM negative
    nes = dict(zip(res.table["term"].astype(str), res.table["NES"]))
    top = next(v for k, v in nes.items() if k.endswith("TOP"))
    bot = next(v for k, v in nes.items() if k.endswith("BOTTOM"))
    assert top > 0 > bot


def test_direction_handles_na_and_zero_nes():
    """NaN NES must not be labelled 'down' (np.where compared NaN as False)."""
    res2d = pd.DataFrame({
        "Term": ["A | pos", "B | neg", "C | nan", "D | zero"],
        "NES": [1.8, -1.4, np.nan, 0.0],
        "NOM p-val": [0.01, 0.02, 0.03, 0.04],
        "FDR q-val": [0.05, 0.06, 0.07, 0.08],
    })
    tidy = gsea._tidy(res2d, n_perm=1000).set_index("term")
    assert tidy.loc["A | pos", "direction"] == "up"
    assert tidy.loc["B | neg", "direction"] == "down"
    # Neither of these is a down-regulated set.
    assert tidy.loc["C | nan", "direction"] == "ns"
    assert tidy.loc["D | zero", "direction"] == "ns"


# --------------------------------------------------------------------------
# Chunked prerank: bounded peak memory, identical per-set statistics
# --------------------------------------------------------------------------
def _synthetic_sets(n_sets=120, n_genes=400, size=25, seed=3):
    rng = np.random.default_rng(seed)
    genes = [f"G{i:04d}" for i in range(n_genes)]
    return genes, {
        f"SET | s{i:03d}": list(rng.choice(genes, size=size, replace=False))
        for i in range(n_sets)
    }


def _synthetic_ranking(genes, seed=4):
    rng = np.random.default_rng(seed)
    return pd.Series(rng.normal(size=len(genes)),
                     index=genes).sort_values(ascending=False)


def test_chunked_matches_whole_run_per_set_statistics():
    """ES / NES / nominal p must not depend on how sets are partitioned.

    gseapy seeds permutations per gene set, so a set's own statistics are
    independent of what it was scored alongside. If this ever fails, the
    chunking is perturbing the RNG and the results cannot be trusted.
    """
    genes, sets = _synthetic_sets()
    rnk = _synthetic_ranking(genes)
    kw = dict(min_size=5, max_size=100, permutations=100, seed=42, threads=1)

    whole = gsea.run_prerank(rnk, sets, **kw).table.set_index("term")
    chunked = gsea.run_prerank_chunked(rnk, sets, chunk_size=25, **kw)\
                  .table.set_index("term")

    common = whole.index.intersection(chunked.index)
    assert len(common) > 20, "too few overlapping terms to be a real check"
    assert len(common) == len(whole) == len(chunked)
    for col in ("ES", "NES", "pval"):
        a = pd.to_numeric(whole.loc[common, col], errors="coerce")
        b = pd.to_numeric(chunked.loc[common, col], errors="coerce")
        assert np.allclose(a, b, rtol=0, atol=1e-6, equal_nan=True), \
            f"{col} drifted between whole and chunked runs"


def test_chunked_fdr_is_global_not_per_chunk():
    """Pooled BH over all nominal p, not a per-block q-value."""
    genes, sets = _synthetic_sets()
    rnk = _synthetic_ranking(genes)
    res = gsea.run_prerank_chunked(
        rnk, sets, min_size=5, max_size=100, permutations=100,
        seed=42, threads=1, chunk_size=25,
    )
    from statsmodels.stats.multitest import multipletests

    tbl = res.table
    # The published pval is already floored, and BH is applied to that floored
    # vector -- adjusting raw p and flooring afterwards would report q at the
    # permutation resolution for the top terms regardless of multiplicity.
    floor = 1.0 / (100 + 1)
    expected = multipletests(
        pd.to_numeric(tbl["pval"], errors="coerce").to_numpy(), method="fdr_bh"
    )[1]
    expected = np.maximum(expected, floor)
    assert np.allclose(tbl["fdr"].to_numpy(), expected, atol=1e-12)
    # Guard the ordering explicitly: no term may claim a q below what the
    # multiplicity burden allows for the best possible p.
    assert tbl["fdr"].min() >= floor
    assert (tbl["fdr_method"] == "BH pooled (chunked run)").all()
    # A per-chunk q-value could not exceed the pooled count; sanity-check that
    # the adjustment really saw every term.
    assert len(tbl) > 25, "global FDR must span more than one chunk"


def test_chunked_keeps_running_curves_for_plots():
    """Phase B must yield a real gseapy object so running plots still work."""
    genes, sets = _synthetic_sets(n_sets=60)
    rnk = _synthetic_ranking(genes)
    res = gsea.run_prerank_chunked(
        rnk, sets, min_size=5, max_size=100, permutations=100,
        seed=42, threads=1, chunk_size=20, max_running_terms=10,
    )
    assert res.raw is not None and isinstance(res.raw.results, dict)
    assert 0 < len(res.running_terms) <= 10
    # result_for_term contract still holds for a single-route result.
    assert gsea.result_for_term(res, res.running_terms[0]) is res
    # and the retained terms really carry running payloads
    payload = res.raw.results[res.running_terms[0]]
    assert "RES" in payload and "hits" in payload
