"""Loading, validating and profiling a DESeq2 results table.

Real-world DESeq2 exports vary: the gene ID may live in the index, column names
differ (`log2FoldChange` vs `logFC`), and `padj` legitimately contains NA from
independent filtering. This module normalises all of that to a canonical schema
and produces a coverage report so the user never has to *guess* how many of
their genes actually carry usable IDs or survive filtering.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

from . import config


@dataclass
class LoadReport:
    """Structured summary of what was loaded and how usable it is."""
    n_rows: int = 0
    n_unique_gene_id: int = 0
    column_mapping: dict = field(default_factory=dict)
    missing_required: list = field(default_factory=list)
    coverage: dict = field(default_factory=dict)         # id -> non-null count
    biotype_counts: dict = field(default_factory=dict)
    n_na_padj: int = 0
    n_tested: int = 0                                    # genes with non-NA padj
    warnings: list = field(default_factory=list)

    def as_text(self) -> str:
        lines = [f"Loaded {self.n_rows} rows ({self.n_unique_gene_id} unique gene IDs)."]
        if self.column_mapping:
            renamed = {k: v for k, v in self.column_mapping.items() if k != v}
            if renamed:
                lines.append("Renamed columns: " + ", ".join(f"{v}->{k}" for k, v in renamed.items()))
        lines.append(
            f"Tested genes (non-NA padj): {self.n_tested}; "
            f"NA padj (dropped from ORA universe): {self.n_na_padj}."
        )
        if self.biotype_counts:
            bt = ", ".join(f"{k}: {v}" for k, v in self.biotype_counts.items())
            lines.append(f"Biotypes: {bt}.")
        for w in self.warnings:
            lines.append(f"[warn] {w}")
        return "\n".join(lines)


def _sniff_sep(path: str) -> str:
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        first = fh.readline()
    if first.count("\t") >= first.count(","):
        return "\t"
    return ","


def _build_rename_map(columns: list[str]) -> dict[str, str]:
    """Map observed columns -> canonical names using the alias table."""
    lower = {c.lower().strip(): c for c in columns}
    rename: dict[str, str] = {}
    for canonical, aliases in config.CANONICAL_COLUMNS.items():
        # exact canonical present?
        if canonical.lower() in lower:
            rename[lower[canonical.lower()]] = canonical
            continue
        for alias in aliases:
            if alias in lower:
                rename[lower[alias]] = canonical
                break
    return rename


def load_deseq2(
    path_or_buffer,
    contrast_name: Optional[str] = None,
    manual_map: Optional[dict[str, str]] = None,
) -> tuple[pd.DataFrame, LoadReport]:
    """Load a DESeq2 results table into a canonical dataframe.

    Parameters
    ----------
    path_or_buffer : str | file-like
        Path or an uploaded file buffer (Streamlit ``UploadedFile``).
    contrast_name : str, optional
        Label attached to every row (used by the multi-contrast comparison).
    manual_map : dict, optional
        Explicit ``{observed_column: canonical_name}`` overrides, applied after
        automatic detection (lets the app offer a "fix my columns" widget).

    Returns
    -------
    (DataFrame, LoadReport)
    """
    report = LoadReport()

    if isinstance(path_or_buffer, str):
        sep = _sniff_sep(path_or_buffer)
        df = pd.read_csv(path_or_buffer, sep=sep)
    else:
        # Uploaded buffer: let pandas infer, fall back to tab.
        try:
            df = pd.read_csv(path_or_buffer, sep=None, engine="python")
        except Exception:
            path_or_buffer.seek(0)
            df = pd.read_csv(path_or_buffer, sep="\t")

    # If the first column looks like an unnamed index carrying gene IDs, name it.
    if df.columns[0] in ("Unnamed: 0", "") or str(df.columns[0]).startswith("Unnamed"):
        df = df.rename(columns={df.columns[0]: "gene_id"})

    rename = _build_rename_map(list(df.columns))
    if manual_map:
        rename.update(manual_map)
    df = df.rename(columns=rename)
    report.column_mapping = {v: k for k, v in rename.items()}

    # Required columns check.
    missing = [c for c in config.REQUIRED_CANONICAL if c not in df.columns]
    report.missing_required = missing
    if missing:
        report.warnings.append(
            "Missing required columns: " + ", ".join(missing) +
            ". Use the column mapper to point at the right columns."
        )
        return df, report

    # Coerce numerics; DESeq2 writes 'NA' strings for filtered padj.
    for col in ("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Normalise ID dtypes to string (entrez can arrive as float with trailing .0).
    if "entrez_id" in df.columns:
        df["entrez_id"] = (
            df["entrez_id"].astype("string").str.replace(r"\.0$", "", regex=True)
        )
        df.loc[df["entrez_id"].isin(["nan", "<NA>", ""]), "entrez_id"] = pd.NA
    df["gene_id"] = df["gene_id"].astype("string")

    if contrast_name:
        df["contrast"] = contrast_name
    elif "contrast" not in df.columns:
        df["contrast"] = "contrast_1"

    # --- Reporting ---
    report.n_rows = len(df)
    report.n_unique_gene_id = int(df["gene_id"].nunique())
    for id_col in ("gene_id", "entrez_id", "gene_name"):
        if id_col in df.columns:
            report.coverage[id_col] = int(df[id_col].notna().sum())
    if "gene_biotype" in df.columns:
        report.biotype_counts = df["gene_biotype"].value_counts(dropna=False).to_dict()
    report.n_na_padj = int(df["padj"].isna().sum())
    report.n_tested = int(df["padj"].notna().sum())

    if report.n_unique_gene_id < report.n_rows:
        report.warnings.append(
            f"{report.n_rows - report.n_unique_gene_id} duplicate gene IDs; "
            "duplicates are collapsed by max |ranking metric| during GSEA."
        )
    if "entrez_id" in df.columns:
        n_missing_entrez = int(df["entrez_id"].isna().sum())
        if n_missing_entrez:
            report.warnings.append(
                f"{n_missing_entrez} genes lack an Entrez ID; g:Profiler will "
                "still map Ensembl IDs, but coverage may drop."
            )
    if "gene_biotype" in df.columns:
        non_pc = int((df["gene_biotype"] != "protein_coding").sum())
        if non_pc:
            report.warnings.append(
                f"{non_pc} non-protein-coding genes (e.g. lncRNA) are mostly "
                "unannotated in GO/KEGG and will drop out of enrichment."
            )

    return df, report
