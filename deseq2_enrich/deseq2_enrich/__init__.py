"""DESeq2 -> ORA/GSEA enrichment pipeline for chicken (and other) RNA-seq.

Public modules:
    config     defaults and source definitions
    io         load/validate a DESeq2 results table (+ coverage report)
    degs       DEG selection and the tested-gene universe
    rank       ranked list for GSEA
    ortho      chicken -> human ortholog mapping (g:Profiler)
    ora        over-representation analysis (g:Profiler, native chicken)
    genesets   fetch MSigDB/Reactome libraries + parse custom GMTs
    gsea       pre-ranked GSEA (gseapy)
    plots      volcano/MA/dotplot/GSEA/heatmap/network/upset figures
    pipeline   end-to-end orchestration used by the CLI and the app
"""
from __future__ import annotations

__version__ = "0.1.0"

from . import config, io, degs, rank, ortho, ora, genesets, gsea, plots  # noqa: F401
