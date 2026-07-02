"""Fetching and parsing gene-set collections for GSEA.

MSigDB / Reactome / WikiPathways gene sets are fetched at runtime from the
Enrichr-hosted libraries via gseapy. This keeps nothing licensed in the repo
(safe for a public app) and avoids shipping large GMT files to a memory-limited
host. Users can also upload their own ``.gmt`` (e.g. curated cGAS-STING,
necroptosis, osteoclast modules) which flows through the identical GSEA path.
"""
from __future__ import annotations

from functools import lru_cache

import gseapy as gp

from . import config


@lru_cache(maxsize=16)
def fetch_library(name: str, organism: str = "human") -> dict:
    """Return an Enrichr-hosted library as ``{term: [genes]}``.

    Cached so repeated GSEA runs in a session hit the network once.
    """
    lib = gp.get_library(name=name, organism=organism)
    return lib


def load_gmt(path_or_lines) -> dict:
    """Parse a GMT file (path or list of lines) into ``{term: [genes]}``.

    GMT format: ``term<TAB>description<TAB>gene1<TAB>gene2...``
    """
    sets: dict[str, list[str]] = {}
    if isinstance(path_or_lines, str):
        with open(path_or_lines, "r", encoding="utf-8", errors="replace") as fh:
            lines = fh.readlines()
    else:
        lines = path_or_lines
    for raw in lines:
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8", errors="replace")
        parts = raw.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        term = parts[0].strip()
        genes = [g.strip() for g in parts[2:] if g.strip()]
        if term and genes:
            sets[term] = genes
    return sets


def combine_libraries(names: list[str], organism: str = "human") -> dict:
    """Fetch and merge multiple Enrichr libraries into one gene-set dict.

    Term names are prefixed with a short library tag so provenance is visible
    in the results table and collisions across libraries are avoided.
    """
    combined: dict[str, list[str]] = {}
    tags = {
        "MSigDB_Hallmark_2020": "HALLMARK",
        "Reactome_2022": "REAC",
        "WikiPathway_2023_Human": "WP",
        "KEGG_2021_Human": "KEGG",
        "MSigDB_Oncogenic_Signatures": "ONCO",
        "MSigDB_Curated_Canonical_Pathways": "C2CP",
    }
    for name in names:
        tag = tags.get(name, name.split("_")[0].upper())
        lib = fetch_library(name, organism)
        for term, genes in lib.items():
            combined[f"{tag} | {term}"] = genes
    return combined
