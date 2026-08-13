"""Fetching and parsing gene-set collections for GSEA.

MSigDB / Reactome / WikiPathways gene sets are fetched at runtime from the
Enrichr-hosted libraries via gseapy. This keeps nothing licensed in the repo
(safe for a public app) and avoids shipping large GMT files to a memory-limited
host. Users can also upload their own ``.gmt`` (e.g. curated cGAS-STING,
necroptosis, osteoclast modules) which flows through the identical GSEA path.
"""
from __future__ import annotations

import sys
import urllib.request
from functools import lru_cache

import gseapy as gp

from . import config
from .netutil import retry_call

# g:Profiler term-id prefixes -> source tag.
_GP_SOURCE_PREFIX = {"GO:": "GO", "KEGG:": "KEGG", "REAC:": "REAC", "WP:": "WP"}
# Aggregator "root" terms that contain most of the genome and mean nothing.
_GMT_ROOT_TERMS = frozenset({"WP:000000", "REAC:0000000", "KEGG:00000", "GO:0003674",
                             "GO:0005575", "GO:0008150"})


@lru_cache(maxsize=config.GSEA_LIBRARY_CACHE_SIZE)
def fetch_library(name: str, organism: str = "human") -> dict:
    """Return an Enrichr-hosted library as ``{term: [genes]}``.

    Cached so repeated GSEA runs in a session hit the network once.
    """
    organism = {"hsapiens": "human", "mmusculus": "mouse"}.get(organism, organism)
    # One retry: an Enrichr blip should not lose a whole GSEA stage.
    return retry_call(gp.get_library, name=name, organism=organism)


def _read_obo(url: str) -> list[str]:
    """Download an OBO file, following redirects, as decoded lines.

    ``current.geneontology.org`` 302s to the release bucket, so redirects must
    be followed -- fetching without them yields a 167-byte stub that parses to
    zero terms rather than failing loudly.
    """
    def _download() -> bytes:
        req = urllib.request.Request(url, headers={"User-Agent": "deseq2_enrich"})
        with urllib.request.urlopen(req, timeout=120) as fh:  # urlopen follows 3xx
            return fh.read()

    return retry_call(_download).decode("utf-8", "replace").splitlines()


@lru_cache(maxsize=1)
def fetch_go_slim_ids() -> frozenset[str]:
    """GO ids in the Consortium's generic slim (~140 broad terms, ~123 KB).

    The slim is what makes GO usable inside the deployed app: it cuts chicken
    GO from ~5,900 scoreable sets to ~82, i.e. seconds instead of minutes and
    well inside the memory ceiling. It is an overview, not a replacement for
    full GO -- report fine-grained terms from a local full-GO run.
    """
    ids = {
        line.strip()[4:]
        for line in _read_obo(config.GO_SLIM_URL)
        if line.startswith("id: GO:")
    }
    if not ids:
        raise RuntimeError(
            f"GO slim from {config.GO_SLIM_URL} parsed to zero terms; the file "
            "format may have changed or a redirect was not followed."
        )
    return frozenset(ids)


@lru_cache(maxsize=1)
def fetch_go_namespaces() -> dict[str, str]:
    """Map GO id -> 'biological_process' / 'molecular_function' / etc.

    The g:Profiler GMT carries no namespace, so without this a request for
    GO:BP silently scores all three branches together. ~31 MB, parsed to a
    ~48k-entry dict and cached for the session.
    """
    ns: dict[str, str] = {}
    current: str | None = None
    for line in _read_obo(config.GO_BASIC_URL):
        line = line.rstrip("\n")
        if line.startswith("id: GO:"):
            current = line[4:]
        elif line.startswith("namespace:") and current:
            ns[current] = line.split(": ", 1)[1].strip()
            current = None
    if not ns:
        raise RuntimeError(
            f"GO ontology from {config.GO_BASIC_URL} parsed to zero namespaces; "
            "the file format may have changed or a redirect was not followed."
        )
    return ns


def _go_term_selector(go_subsets: tuple[str, ...]):
    """Build a predicate over GO term ids for the requested sub-collections."""
    if not go_subsets:
        return lambda _tid: True
    wanted_ns = {
        config.GO_SUBSET_NAMESPACE[s]
        for s in go_subsets
        if s in config.GO_SUBSET_NAMESPACE
    }
    slim_ids = fetch_go_slim_ids() if "slim" in wanted_ns else frozenset()
    namespaces = fetch_go_namespaces() if (wanted_ns - {"slim"}) else {}
    real_ns = wanted_ns - {"slim"}

    def _keep(tid: str) -> bool:
        if tid in slim_ids:
            return True
        return bool(real_ns) and namespaces.get(tid) in real_ns

    return _keep


@lru_cache(maxsize=config.CHICKEN_GMT_CACHE_SIZE)
def fetch_chicken_gmt(
    keying: str = config.CHICKEN_GMT_KEYING,
    sources: tuple[str, ...] = ("GO", "REAC", "WP"),
    go_subsets: tuple[str, ...] = (),
) -> dict:
    """Native chicken gene sets as ``{"TAG | description": [gene_ids]}``.

    Removes the ortholog step entirely for the collections that have chicken
    annotations, which is the single biggest validity win for a chicken-only
    study: no genes are dropped for lacking a human counterpart, so the ranked
    list stays the real universe.

    Parsed streaming and filtered by source during the read -- the full file is
    ~14 MB (symbols) / ~37 MB (Ensembl) and this runs on a ~1 GB host. Gene
    strings are interned because sets overlap heavily; that cuts the retained
    dict substantially.

    Raises ``RuntimeError`` on any fetch/parse failure so the caller can fall
    back to an uploaded GMT rather than reporting an empty run as success.
    """
    if keying not in config.GPROFILER_GMT_URLS:
        raise ValueError(
            f"Unknown chicken GMT keying {keying!r}; "
            f"expected one of {tuple(config.GPROFILER_GMT_URLS)}"
        )
    url = config.GPROFILER_GMT_URLS[keying]
    wanted = set(sources or ())
    # Resolved before the GMT loop so a bad ontology fetch fails fast rather
    # than after parsing 14 MB.
    keep_go = _go_term_selector(tuple(go_subsets or ()))

    def _download() -> bytes:
        with urllib.request.urlopen(url, timeout=60) as fh:
            return fh.read()

    try:
        payload = retry_call(_download)
    except Exception as exc:  # network, DNS, TLS, HTTP error
        raise RuntimeError(
            f"Could not download the native chicken gene sets from {url} "
            f"({type(exc).__name__}: {exc}). Upload a GMT file instead, or "
            "switch the GSEA mode to ortholog."
        ) from exc

    sets: dict[str, list[str]] = {}
    for raw in payload.decode("utf-8", "replace").splitlines():
        parts = raw.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        term_id = parts[0].strip()
        if term_id in _GMT_ROOT_TERMS:
            continue
        tag = next(
            (t for pre, t in _GP_SOURCE_PREFIX.items() if term_id.startswith(pre)),
            None,
        )
        if tag is None or (wanted and tag not in wanted):
            continue
        if tag == "GO" and not keep_go(term_id):
            continue
        genes = [sys.intern(g.strip()) for g in parts[2:] if g.strip()]
        if not genes:
            continue
        desc = parts[1].strip() or term_id
        sets[f"{tag} | {desc}"] = genes

    if not sets:
        raise RuntimeError(
            f"Native chicken gene sets from {url} parsed to zero usable terms "
            f"for sources {sorted(wanted)}"
            + (f" / GO subsets {sorted(go_subsets)}" if go_subsets else "")
            + ". The file format may have changed."
        )
    return sets


def clear_cache() -> None:
    fetch_library.cache_clear()
    fetch_chicken_gmt.cache_clear()
    fetch_go_slim_ids.cache_clear()
    fetch_go_namespaces.cache_clear()


def cache_info():
    return fetch_library.cache_info()


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
        "GO_Biological_Process_2026": "GO:BP",
        "GO_Molecular_Function_2026": "GO:MF",
        "GO_Cellular_Component_2026": "GO:CC",
    }
    for name in names:
        tag = tags.get(name, name.split("_")[0].upper())
        lib = fetch_library(name, organism)
        for term, genes in lib.items():
            combined[f"{tag} | {term}"] = genes
    return combined
