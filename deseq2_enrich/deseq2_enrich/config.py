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

# --- Orthology ------------------------------------------------------------
# Which identifier column to hand the orthology endpoint, best first. Stable
# Ensembl IDs (ENSGALG...) map far more completely and less ambiguously than
# chicken symbols, which are frequently LOC*/unannotated or clash with human
# symbols by coincidence. The mapper tries these in order and keeps whichever
# yields the most human symbols.
ORTHO_ID_PRIORITY = ["ensembl_id", "gene_id", "gene_name"]

# Guardrails: below either floor the ranked list is too small/biased for the
# pre-ranked GSEA null to mean anything, so the pipeline raises instead of
# silently reporting enrichment computed on a handful of genes.
ORTHO_MIN_MAPPING_RATE = 0.30   # fraction of queried IDs that gain a human ortholog
ORTHO_MIN_MAPPED_GENES = 200    # absolute floor on mapped source IDs

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
    "GO_Biological_Process_2026": "GO:BP Biological Process (2026)",
    "GO_Molecular_Function_2026": "GO:MF Molecular Function (2026)",
    "GO_Cellular_Component_2026": "GO:CC Cellular Component (2026)",
    "GO_Slim_Chicken": "GO slim (chicken, fast overview)",
}
# Both of these route to native chicken gene sets under the default 'auto'
# mode, so a default run performs no ortholog mapping at all -- no dropout
# bias, and the sparse-mapping guardrail never fires on the default path.
# Hallmark and Oncogenic Signatures are human-only and remain selectable;
# ticking either one adds a second, ortholog-based prerank with its own FDR.
GSEA_DEFAULT_LIBRARIES = ["Reactome_2022", "WikiPathway_2023_Human"]

# Native GO is ~18.7k chicken terms and dominates prerank runtime: measured on
# a 16.3k-gene table, GO alone scores 5924 sets in ~6 min at 1000 permutations
# versus ~40 s for Reactome + WikiPathways. It is therefore opt-in rather than
# a default, and the UI states the cost before you tick it.
GSEA_SLOW_LIBRARIES = {
    "GO_Biological_Process_2026",
    "GO_Molecular_Function_2026",
    "GO_Cellular_Component_2026",
}
GSEA_LIBRARY_GROUPS = {
    "Curated pathways": [
        "MSigDB_Hallmark_2020",
        "Reactome_2022",
        "KEGG_2021_Human",
        "WikiPathway_2023_Human",
        "MSigDB_Oncogenic_Signatures",
    ],
    "Gene Ontology": [
        "GO_Slim_Chicken",
        "GO_Biological_Process_2026",
        "GO_Molecular_Function_2026",
        "GO_Cellular_Component_2026",
    ],
}

# --- Native-chicken GSEA (g:Profiler GMT) --------------------------------
# GO/Reactome/WikiPathways all have native chicken annotations, so prerank can
# be run on them directly -- no ortholog step and none of its dropout bias.
# Only genuinely human-only collections (Hallmark, Oncogenic) still need it.
#
# URLs verified live 2026-08-13: the '_full_' filenames resolve (HTTP 206 on a
# range request, last modified 2026-03-20). The shorter
# 'gprofiler_ggallus.<keying>.gmt' names return 404 and must not be used.
GPROFILER_GMT_URLS = {
    "name": "https://biit.cs.ut.ee/gprofiler/static/gprofiler_full_ggallus.name.gmt",
    "ensg": "https://biit.cs.ut.ee/gprofiler/static/gprofiler_full_ggallus.ENSG.gmt",
}
# Symbol keying, NOT Ensembl. The GMT ships GRCg7b IDs (ENSGALG000100...) while
# real exports are still commonly galgal6 (ENSGALG000000...); measured overlap
# between the two series is 0/5000. Symbols overlap 84% on the same data.
CHICKEN_GMT_KEYING = "name"
CHICKEN_GMT_KEYINGS = ("name", "ensg")

# Sources kept from the chicken GMT. HP (human phenotype) and MIRNA are
# dropped: HP is meaningless for chicken and MIRNA sets are almost all too
# small to score. KEGG is absent from g:Profiler's downloadable GMTs
# (licensing) -- it remains available in the ORA tab, which queries the API.
NATIVE_GSEA_SOURCES = {
    "GO": "Gene Ontology (chicken)",
    "REAC": "Reactome (chicken)",
    "WP": "WikiPathways (chicken)",
}
NATIVE_GSEA_DEFAULT_SOURCES = ["REAC", "WP"]

# --- GO ontology helpers --------------------------------------------------
# The g:Profiler GMT does not encode the GO namespace in the term id, so
# "GO:BP" previously scored all of GO -- BP, MF and CC pooled. go-basic.obo
# supplies the real namespace per term. goslim_generic.obo is the GO
# Consortium's curated ~140-term overview subset; 82 of those survive in
# chicken at the default size bounds, which is what makes GO tractable in the
# deployed app at all (0.7 s / 391 MB at 100 permutations, versus 5,924 sets
# for full GO). Both URLs verified live 2026-08-13; follow redirects.
GO_SLIM_URL = "http://current.geneontology.org/ontology/subsets/goslim_generic.obo"
GO_BASIC_URL = "http://current.geneontology.org/ontology/go-basic.obo"

# Native GO sub-collections -> the GO namespace they select ("slim" is the
# curated subset rather than a namespace).
GO_SUBSET_NAMESPACE = {
    "GO:BP": "biological_process",
    "GO:MF": "molecular_function",
    "GO:CC": "cellular_component",
    "GO:SLIM": "slim",
}

# Enrichr library -> native chicken source/sub-collection, used by 'auto'
# routing. The GO entries resolve to real namespaces via go-basic.obo (see
# GO_SUBSET_NAMESPACE), so GO:BP now means biological_process only rather than
# all three branches pooled.
LIBRARY_NATIVE_SOURCE = {
    "Reactome_2022": "REAC",
    "WikiPathway_2023_Human": "WP",
    "GO_Biological_Process_2026": "GO:BP",
    "GO_Molecular_Function_2026": "GO:MF",
    "GO_Cellular_Component_2026": "GO:CC",
    "GO_Slim_Chicken": "GO:SLIM",
}

# Selectable in the app but not an Enrichr library, so it cannot be routed
# through the human-ortholog path.
NATIVE_ONLY_LIBRARIES = ("GO_Slim_Chicken",)
# Human-only concepts with no chicken equivalent; these keep the ortholog route.
HUMAN_ONLY_LIBRARIES = ("MSigDB_Hallmark_2020", "MSigDB_Oncogenic_Signatures")

GSEA_MODES = ("auto", "native_chicken", "ortholog")
GSEA_DEFAULT_MODE = "auto"
CHICKEN_GMT_CACHE_SIZE = 1  # a parsed GMT is large; hold one at a time

# --- GSEA parameters ------------------------------------------------------
# gseapy filters each set to min_size <= |set intersect ranked_list| <= max_size
# AFTER intersecting with the ranked list, then silently drops the rest. 15 wipes
# out otherwise-valid small pathways once the list has shrunk, so chicken runs
# start at 10. Both bounds are surfaced in the UI rather than left as silent
# constants -- how many sets actually got scored is part of reading the result.
GSEA_MIN_SIZE = 10
GSEA_MAX_SIZE = 500
GSEA_PERMUTATIONS = 1000
GSEA_SEED = 42
GSEA_LIBRARY_CACHE_SIZE = 4
ORTHOLOG_CACHE_SIZE = 2
ORA_CACHE_SIZE = 8
GSEA_MAX_RUNNING_PLOT_TERMS = 250

# --- Chunked prerank (memory ceiling) -------------------------------------
# gseapy materialises a running-enrichment curve of length |ranked list| for
# EVERY gene set inside a single gp.prerank call. At ~4,400 GO:BP sets over a
# ~16k-gene ranking that is ~0.5 GB of curves alone, and the measured peak is
# far higher once permutation structures and the transient upper-cased copy of
# the gene-set dict are included -- well past Community Cloud's ~1 GB.
#
# The peak happens INSIDE gp.prerank, so _compact_raw (which prunes after the
# call returns) cannot prevent it. Scoring in blocks bounds the peak to one
# block. Verified empirically: ES, NES and nominal p are bit-identical whole
# vs chunked because gseapy seeds its permutations per gene set. Only the
# q-value differs, since gseapy computes it across whatever is in the run --
# which is why FDR is recomputed globally afterwards.
# 100 chosen from measurement, not intuition. Full GO:BP (4,372 scored sets,
# 16.3k-gene ranking), peak RSS by chunk size at 1000 permutations:
#   whole run 3,723 MB | 500 -> 2,380 | 100 -> 1,822 | 25 -> 1,693
# Returns diminish below ~100 while runtime keeps climbing, because what
# remains is gseapy's per-set transient, not accumulation across sets.
GSEA_CHUNK_SIZE = 100

# Peak scales with PERMUTATION count, which is the real ceiling. Chunked
# GO:BP peak RSS: 100 perms -> 879 MB | 250 -> 1,057 | 500 -> 1,425 |
# 1000 -> 1,822. So GO:BP fits a ~1 GB host only at ~100 permutations; above
# that it belongs in a local CLI run. The app warns past this.
GSEA_INAPP_PERMUTATION_CEILING = 100
# Only chunk above this many sets; smaller libraries (Reactome, WikiPathways,
# custom modules) stay on the original single-call path unchanged.
GSEA_CHUNK_THRESHOLD = 800

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
    # A *separate* Ensembl column (common in tximport/biomaRt exports where
    # gene_id already holds a symbol). Deliberately does not share aliases with
    # gene_id so the two canonicals can never claim the same observed column.
    "ensembl_id": ["ensembl_gene_id", "ensembl_gene", "ensemblid", "ensgene", "ensg"],
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
