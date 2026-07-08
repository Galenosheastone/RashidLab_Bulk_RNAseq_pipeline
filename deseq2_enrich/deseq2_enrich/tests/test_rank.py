import pandas as pd

from deseq2_enrich import rank


def test_collapse_by_max_abs_metric():
    df = pd.DataFrame({
        "human_symbol": ["FOO", "FOO", "FOO", "BAR"],
        "stat": [2.0, -5.0, 3.0, 1.0],
        "log2FoldChange": [1.0, -1.0, 1.5, 0.5],
        "pvalue": [0.1, 0.01, 0.05, 0.5],
    })
    r = rank.build_rank(df, metric="stat", key_col="human_symbol")
    assert r["FOO"] == -5.0
    assert r["BAR"] == 1.0
    assert list(r.index) == ["BAR", "FOO"]


def test_signed_logp_handles_p_zero():
    df = pd.DataFrame({
        "human_symbol": ["A", "B"],
        "log2FoldChange": [2.0, -2.0],
        "pvalue": [0.0, 1e-10],
        "stat": [10.0, -10.0],
    })
    r = rank.build_rank(df, metric="signed_logp", key_col="human_symbol")
    assert r["A"] > 0 and abs(r["A"]) > 100
    assert r["B"] < 0
