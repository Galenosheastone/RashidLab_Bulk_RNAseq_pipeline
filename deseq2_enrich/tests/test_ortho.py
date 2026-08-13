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
        mapped, report = ortho.attach_human_symbol(df, id_col="gene_id")

    assert mapped["human_symbol"].tolist() == ["COL10A1", "FRMD8"]
    assert mapped["human_symbol_source"].tolist() == ["gene_name", "gene_name"]
    assert report.id_col_used == "gene_name"
    assert report.n_query_mapped == 2
    assert report.mapping_rate == 1.0


def test_attach_human_symbol_prefers_ensembl_over_symbol():
    """A GRCg7b export with symbols in gene_id must still key on ENSGALG IDs."""
    df = pd.DataFrame({
        "gene_id": ["MT-CO1", "MYOD1"],
        "ensembl_id": ["ENSGALG00000032142", "ENSGALG00000012345"],
        "gene_name": ["MT-CO1", "MYOD1"],
        "stat": [3.5, -4.4],
    })

    def fake_orth_cached(genes_key, source, target):
        genes = list(genes_key)
        # Only the Ensembl IDs resolve; symbols come back empty.
        if not all(g.startswith("ENSGALG") for g in genes):
            return pd.DataFrame(columns=["incoming", "ortholog_name"])
        return pd.DataFrame({
            "incoming": genes,
            "ortholog_name": ["MT-CO1", "MYOD1"],
        })

    with patch("deseq2_enrich.ortho._orth_cached", side_effect=fake_orth_cached):
        mapped, report = ortho.attach_human_symbol(df, id_col="gene_id")

    assert report.id_col_used == "ensembl_id"
    assert mapped["human_symbol"].tolist() == ["MT-CO1", "MYOD1"]


def test_map_to_human_reports_multiplicity_and_strict_filter():
    """strict_one_to_one must drop both one->many and many->one pairs."""
    pairs = pd.DataFrame({
        "incoming": ["ENSGALG1", "ENSGALG1", "ENSGALG2", "ENSGALG3", "ENSGALG4"],
        "ortholog_name": ["PAX1", "PAX9", "MERGED", "MERGED", "CLEAN"],
    })

    with patch("deseq2_enrich.ortho._orth_cached", return_value=pairs):
        mapping, mult = ortho.map_to_human(
            ["ENSGALG1", "ENSGALG2", "ENSGALG3", "ENSGALG4"]
        )
        strict, _ = ortho.map_to_human(
            ["ENSGALG1", "ENSGALG2", "ENSGALG3", "ENSGALG4"],
            strict_one_to_one=True,
        )

    assert mult == {"one_to_many": 1, "many_to_one": 1}
    assert len(mapping) == 5
    # ENSGALG1 fans out, ENSGALG2/3 collapse onto MERGED; only ENSGALG4 is clean.
    assert strict["incoming"].tolist() == ["ENSGALG4"]
    assert strict["ortholog_name"].tolist() == ["CLEAN"]
