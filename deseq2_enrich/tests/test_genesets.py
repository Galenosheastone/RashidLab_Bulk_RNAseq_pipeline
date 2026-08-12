from unittest.mock import patch

from deseq2_enrich import config, genesets


def test_go_gsea_libraries_are_grouped_and_labeled():
    assert "GO_Biological_Process_2026" in config.GSEA_LIBRARIES
    assert "GO_Molecular_Function_2026" in config.GSEA_LIBRARIES
    assert "GO_Cellular_Component_2026" in config.GSEA_LIBRARIES
    assert config.GSEA_LIBRARY_GROUPS["Gene Ontology"] == [
        "GO_Biological_Process_2026",
        "GO_Molecular_Function_2026",
        "GO_Cellular_Component_2026",
    ]


def test_combine_libraries_prefixes_go_terms():
    with patch(
        "deseq2_enrich.genesets.fetch_library",
        return_value={"positive regulation": ["A", "B", "C"]},
    ):
        combined = genesets.combine_libraries(["GO_Biological_Process_2026"])

    assert combined == {"GO:BP | positive regulation": ["A", "B", "C"]}
