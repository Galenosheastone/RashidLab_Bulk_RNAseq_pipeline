from unittest.mock import patch

import pandas as pd

from deseq2_enrich import ortho


def test_attach_human_symbol_falls_back_to_gene_name_when_gene_id_unmapped():
    df = pd.DataFrame({
        "gene_id": ["100858979", "107051366"],
        "gene_name": ["COL10A1", "FRMD8"],
        "stat": [3.5, -4.4],
    })

    def fake_orth_cached(genes_key, source, target):
        genes = list(genes_key)
        if all(g.isdigit() for g in genes):
            return pd.DataFrame(columns=["incoming", "ortholog_name"])
        return pd.DataFrame({
            "incoming": genes,
            "ortholog_name": genes,
        })

    with patch("deseq2_enrich.ortho._orth_cached", side_effect=fake_orth_cached):
        mapped = ortho.attach_human_symbol(df, id_col="gene_id")

    assert mapped["human_symbol"].tolist() == ["COL10A1", "FRMD8"]
    assert mapped["human_symbol_source"].tolist() == ["gene_name", "gene_name"]
