import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def toy_deseq2():
    """Small DESeq2-like table with NA padj rows and a duplicate up gene."""
    rng = np.random.default_rng(0)
    n = 20
    df = pd.DataFrame({
        "gene_id": [f"G{i:03d}" for i in range(n)],
        "gene_name": [f"gene_{i}" for i in range(n)],
        "baseMean": rng.uniform(10, 5000, n),
        "log2FoldChange": np.concatenate([
            [2.5, 3.0, 1.5, 1.2],
            [-2.0, -3.5, -1.4],
            rng.uniform(-0.8, 0.8, n - 12),
            [0.1, 0.2, 0.3, 0.4, 0.5],
        ]),
        "lfcSE": 0.3,
        "stat": np.concatenate([
            [10, 12, 8, 7],
            [-9, -13, -6],
            rng.normal(0, 1, n - 12),
            [0.5, 0.6, 0.7, 0.8, 0.9],
        ]),
        "pvalue": np.concatenate([
            [1e-20, 1e-25, 1e-15, 1e-10],
            [1e-18, 1e-30, 1e-12],
            rng.uniform(0.2, 0.9, n - 12),
            [np.nan] * 5,
        ]),
        "padj": np.concatenate([
            [1e-19, 1e-24, 1e-14, 1e-9],
            [1e-17, 1e-29, 1e-11],
            rng.uniform(0.3, 0.95, n - 12),
            [np.nan] * 5,
        ]),
    })
    df.loc[1, "gene_id"] = df.loc[0, "gene_id"]
    return df
