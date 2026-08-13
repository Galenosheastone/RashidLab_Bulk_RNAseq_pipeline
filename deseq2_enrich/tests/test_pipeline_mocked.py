from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from deseq2_enrich import config
from deseq2_enrich.pipeline import run_contrast


@pytest.fixture(autouse=True)
def _toy_scale_guardrails():
    """The toy fixtures are ~20 genes; scale the sparse-mapping floor to match.

    The production floor (200 mapped genes) exists to stop a real run from
    reporting GSEA on a remnant list -- it is not a statement about fixtures.
    """
    with patch.object(config, "ORTHO_MIN_MAPPED_GENES", 5):
        yield


def _fake_ora_result(*args, **kwargs):
    return pd.DataFrame({
        "source": ["GO:BP"],
        "native": ["GO:0001"],
        "name": ["fake term"],
        "p_value": [1e-5],
        "term_size": [50],
        "query_size": [10],
        "intersection_size": [5],
        "intersections": [["FAKE1", "FAKE2", "FAKE3", "FAKE4", "FAKE5"]],
    })


def _fake_orth_result(*args, **kwargs):
    genes = list(args[0] if args else kwargs.get("query", []))
    return pd.DataFrame({
        "incoming": genes,
        "ortholog_name": [f"HS_{g}" for g in genes],
    })


def test_pipeline_end_to_end_mocked(tmp_path, toy_deseq2):
    csv_path = tmp_path / "toy.tsv"
    toy_deseq2.to_csv(csv_path, sep="\t", index=False)

    with patch("deseq2_enrich.ora._profile_cached", side_effect=_fake_ora_result), \
         patch("deseq2_enrich.ortho._orth_cached", side_effect=_fake_orth_result), \
         patch(
             "deseq2_enrich.genesets.fetch_library",
             return_value={"HALLMARK_FAKE": ["HS_G000", "HS_G001", "HS_G002"] * 10},
         ), \
         patch("deseq2_enrich.gsea.gp.prerank") as mock_prerank:
        fake_pre = MagicMock()
        fake_pre.res2d = pd.DataFrame({
            "Term": ["HALLMARK | fake"],
            "ES": [0.5],
            "NES": [1.8],
            "NOM p-val": [0.01],
            "FDR q-val": [0.02],
            "FWER p-val": [0.01],
            "Lead_genes": ["HS_G000;HS_G001"],
        })
        fake_pre.results = {
            "HALLMARK | fake": {
                "nes": 1.8,
                "fdr": 0.02,
                "RES": [0.0, 0.3, 0.5, 0.2],
                "hits": [1, 2],
                "lead_genes": "HS_G000;HS_G001",
            }
        }
        mock_prerank.return_value = fake_pre

        res = run_contrast(
            str(csv_path),
            contrast_name="toy",
            gsea_libraries=["MSigDB_Hallmark_2020"],
            gsea_permutations=10,
        )

    assert res.errors == {}
    assert len(res.ora) > 0
    assert res.gsea is not None
    assert len(res.gsea.table) == 1
    assert res.gsea_metadata["libraries"] == ["MSigDB_Hallmark_2020"]
    assert res.gsea_metadata["ranking_size"] == len(res.ranking)
    # The mapping rate must always be reported, not just on failure.
    assert res.gsea_metadata["mapping_rate"] == 1.0
    assert "Ortholog mapping via" in res.gsea_metadata["ortholog_report"]


def _sparse_orth_result(genes_key, *args, **kwargs):
    """Only the first query gene has a human ortholog."""
    genes = list(genes_key)
    return pd.DataFrame({
        "incoming": genes[:1],
        "ortholog_name": ["HS_ONLY_ONE"],
    })


def test_gsea_raises_when_ortholog_mapping_is_sparse(tmp_path, toy_deseq2):
    """A near-empty mapping must fail loudly, not produce a meaningless GSEA."""
    csv_path = tmp_path / "toy.tsv"
    toy_deseq2.to_csv(csv_path, sep="\t", index=False)

    with patch("deseq2_enrich.ortho._orth_cached", side_effect=_sparse_orth_result), \
         patch(
             "deseq2_enrich.genesets.fetch_library",
             return_value={"HALLMARK_FAKE": ["HS_ONLY_ONE"] * 30},
         ), \
         patch("deseq2_enrich.gsea.gp.prerank") as mock_prerank:
        res = run_contrast(
            str(csv_path),
            contrast_name="toy",
            gsea_libraries=["MSigDB_Hallmark_2020"],
            gsea_permutations=10,
            do_ora=False,
        )

    assert res.gsea is None
    assert mock_prerank.call_count == 0, "prerank must not run on a sparse mapping"
    assert "ValueError" in res.errors["gsea"]
    assert "too sparse" in res.errors["gsea"]
    # The rate is recorded even though the run aborted, so the UI can show it.
    assert res.gsea_metadata["mapping_rate"] < config.ORTHO_MIN_MAPPING_RATE


# --------------------------------------------------------------------------
# Native-chicken routing (Task 1.1)
# --------------------------------------------------------------------------
_FAKE_CHICKEN_GMT = {
    "GO | fake chicken process": [f"gene_{i}" for i in range(12)],
    "REAC | fake chicken pathway": [f"gene_{i}" for i in range(5, 18)],
}


def _fake_prerank_factory():
    """A prerank stand-in that echoes back whichever terms it was given."""
    def _fake(rnk, gene_sets, **kwargs):
        terms = list(gene_sets)
        pre = MagicMock()
        pre.res2d = pd.DataFrame({
            "Term": terms,
            "ES": [0.4] * len(terms),
            "NES": [1.5] * len(terms),
            "NOM p-val": [0.01] * len(terms),
            "FDR q-val": [0.02] * len(terms),
            "FWER p-val": [0.01] * len(terms),
            "Lead_genes": ["gene_1;gene_2"] * len(terms),
        })
        pre.results = {
            t: {"nes": 1.5, "fdr": 0.02, "RES": [0.0, 0.3], "hits": [1],
                "lead_genes": "gene_1;gene_2"}
            for t in terms
        }
        return pre
    return _fake


def test_native_mode_skips_ortholog_mapping(tmp_path, toy_deseq2):
    """Native chicken must never touch the orthology endpoint."""
    csv_path = tmp_path / "toy.tsv"
    toy_deseq2.to_csv(csv_path, sep="\t", index=False)

    with patch("deseq2_enrich.genesets.fetch_chicken_gmt",
               return_value=dict(_FAKE_CHICKEN_GMT)) as mock_gmt, \
         patch("deseq2_enrich.ortho.attach_human_symbol") as mock_ortho, \
         patch("deseq2_enrich.gsea.gp.prerank", side_effect=_fake_prerank_factory()):
        res = run_contrast(
            str(csv_path), contrast_name="toy", do_ora=False,
            gsea_libraries=["GO_Biological_Process_2026"],
            gsea_mode="native_chicken", gsea_permutations=10,
        )

    assert res.errors == {}, res.errors
    mock_ortho.assert_not_called()
    assert mock_gmt.call_count == 1
    assert res.gsea_metadata["gsea_mode"] == "native_chicken"
    assert res.gsea_metadata["gsea_routes"] == ["native"]
    # Ranked on chicken symbols, not human orthologs.
    assert res.gsea_metadata["native_key_col"] in ("gene_name", "gene_id")
    assert set(res.gsea.table["gsea_route"]) == {"native"}


def test_auto_mode_runs_both_routes_with_separate_fdr(tmp_path, toy_deseq2):
    """GO goes native, Hallmark goes through orthologs, FDR stays per route."""
    csv_path = tmp_path / "toy.tsv"
    toy_deseq2.to_csv(csv_path, sep="\t", index=False)

    with patch("deseq2_enrich.genesets.fetch_chicken_gmt",
               return_value=dict(_FAKE_CHICKEN_GMT)), \
         patch("deseq2_enrich.ortho._orth_cached", side_effect=_fake_orth_result), \
         patch("deseq2_enrich.genesets.fetch_library",
               return_value={"HALLMARK_FAKE": ["HS_G000", "HS_G001"] * 8}), \
         patch("deseq2_enrich.gsea.gp.prerank", side_effect=_fake_prerank_factory()):
        res = run_contrast(
            str(csv_path), contrast_name="toy", do_ora=False,
            gsea_libraries=["GO_Biological_Process_2026", "MSigDB_Hallmark_2020"],
            gsea_mode="auto", gsea_permutations=10,
        )

    assert res.errors == {}, res.errors
    assert res.gsea_metadata["gsea_routes"] == ["native", "ortholog"]
    assert res.gsea_metadata["fdr_per_route"] is True
    routes = set(res.gsea.table["gsea_route"])
    assert routes == {"native", "ortholog"}
    # Both routes are retained separately so running plots use the right ranking.
    assert set(res.gsea.routes) == {"native", "ortholog"}
    # The ortholog route still reports its mapping rate.
    assert "ortholog_report" in res.gsea_metadata


def test_native_gmt_failure_is_reported_not_swallowed(tmp_path, toy_deseq2):
    """A dead GMT URL must surface as an error, never as an empty success."""
    csv_path = tmp_path / "toy.tsv"
    toy_deseq2.to_csv(csv_path, sep="\t", index=False)

    with patch("deseq2_enrich.genesets.fetch_chicken_gmt",
               side_effect=RuntimeError("Could not download the native chicken gene sets")):
        res = run_contrast(
            str(csv_path), contrast_name="toy", do_ora=False,
            gsea_libraries=["Reactome_2022"],
            gsea_mode="native_chicken", gsea_permutations=10,
        )

    assert res.gsea is None
    assert "Could not download" in res.errors["gsea"]
