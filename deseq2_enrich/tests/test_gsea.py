from types import SimpleNamespace
from unittest.mock import patch

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
