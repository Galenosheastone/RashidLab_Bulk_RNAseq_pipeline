from unittest.mock import MagicMock, patch

import pandas as pd

from deseq2_enrich.pipeline import run_contrast


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
