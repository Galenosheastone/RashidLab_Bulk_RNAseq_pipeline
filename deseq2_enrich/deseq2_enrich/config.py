"""Central configuration and defaults for the DESeq2 enrichment pipeline.

Keeping all tunable defaults in one place makes the CLI and the Streamlit app
behave identically and keeps the "reproduce my existing DEG column" promise
honest: the defaults below use padj < 0.05 with no fold-change cutoff.
"""
from __future__ import annotations

# --- Organism -------------------------------------------------------------
# g:Profiler organism code for chicken. Change here if you ever reuse the app
# for another species (e.g. 'hsapiens', 'mmusculus').
ORGANISM = "ggallus"
ORTHOLOG_TARGET = "hsapiens"  # GSEA is run against human MSigDB via orthologs

# --- DEG selection defaults ----------------------------------------------
PADJ_THRESHOLD = 0.05
LFC_THRESHOLD = 0.0

# --- GSEA ranking ---------------------------------------------------------
# 'stat'        : DESeq2 Wald statistic (recommended; signed, magnitude-aware)
# 'signed_logp' : sign(log2FC) * -log10(pvalue)   (p floored to avoid inf)
# 'log2fc'      : raw log2 fold change
RANK_METRIC = "stat"
RANK_METRICS = ("stat", "signed_logp", "log2fc")
PVALUE_FLOOR = 1e-300  # floor for -log10(p) so p==0 does not produce inf

# --- ORA (g:Profiler) native-chicken sources -----------------------------
# Keys are g:Profiler source codes; values are human-readable labels.
ORA_SOURCES = {
    "GO:BP": "GO Biological Process",
    "GO:MF": "GO Molecular Function",
    "GO:CC": "GO Cellular Component",
    "KEGG": "KEGG",
    "REAC": "Reactome",
    "WP": "WikiPathways",
}
ORA_DEFAULT_SOURCES = ["GO:BP", "KEGG", "REAC", "WP"]

# --- GSEA gene-set libraries (Enrichr-hosted MSigDB, human symbols) -------
# Fetched at runtime so nothing licensed is committed to a public repo.
GSEA_LIBRARIES = {
    "MSigDB_Hallmark_2020": "MSigDB Hallmark",
    "Reactome_2022": "Reactome (2022)",
    "WikiPathway_2023_Human": "WikiPathways (Human 2023)",
    "MSigDB_Oncogenic_Signatures": "MSigDB Oncogenic",
    "KEGG_2021_Human": "KEGG (Human 2021)",
}
GSEA_DEFAULT_LIBRARIES = ["MSigDB_Hallmark_2020", "Reactome_2022"]

# --- GSEA parameters ------------------------------------------------------
GSEA_MIN_SIZE = 15
GSEA_MAX_SIZE = 500
GSEA_PERMUTATIONS = 1000
GSEA_SEED = 42

# --- Significance display -------------------------------------------------
SIG_ALPHA = 0.05  # adj-p / FDR line used in plots and "top term" tables

# --- Palette (colour-blind-safe-ish, consistent across all figures) -------
COLOR_UP = "#C0392B"     # up-regulated
COLOR_DOWN = "#2471A3"   # down-regulated
COLOR_NS = "#B0B0B0"     # not significant
COLOR_ACCENT = "#1E8449"
CONTINUOUS_SCALE = "RdBu_r"  # for NES / signed statistics

# --- Required / optional columns of a DESeq2 export -----------------------
# The loader tries to normalise to these canonical names.
CANONICAL_COLUMNS = {
    "gene_id": ["gene_id", "ensembl", "ensembl_id", "gene", "geneid", "id", "row"],
    "entrez_id": ["entrez_id", "entrez", "entrezid", "ncbi_id"],
    "gene_name": ["gene_name", "symbol", "gene_symbol", "name", "external_gene_name"],
    "gene_biotype": ["gene_biotype", "biotype"],
    "baseMean": ["basemean", "base_mean"],
    "log2FoldChange": ["log2foldchange", "log2fc", "logfc", "log2_fold_change"],
    "lfcSE": ["lfcse", "lfc_se", "se"],
    "stat": ["stat", "wald", "wald_stat", "statistic"],
    "pvalue": ["pvalue", "p_value", "pval", "p"],
    "padj": ["padj", "p_adj", "fdr", "qvalue", "adj_pvalue", "padjust"],
}
REQUIRED_CANONICAL = ["gene_id", "log2FoldChange", "pvalue", "padj"]
