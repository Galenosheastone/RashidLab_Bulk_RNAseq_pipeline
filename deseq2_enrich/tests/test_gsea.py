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
