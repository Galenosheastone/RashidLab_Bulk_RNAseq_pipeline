"""Command-line runner for reproducible batch enrichment.

Runs the same pipeline the app uses, but headless — suited to Tempest/HPC where
you process the real 83-sample multi-contrast design and want committed outputs.

Example
-------
    python -m deseq2_enrich.cli \\
        --input results/Sacral_vs_Caudal.tsv \\
        --name Sacral_vs_Caudal \\
        --outdir out/ \\
        --sources GO:BP KEGG REAC WP \\
        --gsea-libs MSigDB_Hallmark_2020 Reactome_2022 \\
        --rank stat --padj 0.05 --lfc 1.0

Writes: <name>_ORA.csv, <name>_GSEA.csv, and interactive HTML figures.
"""
from __future__ import annotations

import argparse
import os

from . import config, degs, plots
from .pipeline import run_contrast


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="DESeq2 -> ORA/GSEA enrichment")
    p.add_argument("--input", required=True, help="DESeq2 results TSV/CSV")
    p.add_argument("--name", default="contrast_1", help="contrast label")
    p.add_argument("--outdir", default="enrich_out")
    p.add_argument("--organism", default=config.ORGANISM)
    p.add_argument("--padj", type=float, default=config.PADJ_THRESHOLD)
    p.add_argument("--lfc", type=float, default=config.LFC_THRESHOLD)
    p.add_argument("--rank", default=config.RANK_METRIC, choices=config.RANK_METRICS)
    p.add_argument("--id-col", default="gene_id")
    p.add_argument("--sources", nargs="+", default=config.ORA_DEFAULT_SOURCES)
    p.add_argument("--gsea-libs", nargs="+", default=config.GSEA_DEFAULT_LIBRARIES)
    p.add_argument("--custom-gmt", default=None, help="optional custom .gmt path")
    p.add_argument("--permutations", type=int, default=config.GSEA_PERMUTATIONS)
    p.add_argument("--no-ora", action="store_true")
    p.add_argument("--no-gsea", action="store_true")
    p.add_argument("--gsea-mode", default=config.GSEA_DEFAULT_MODE,
                   choices=list(config.GSEA_MODES),
                   help="'auto' (default) scores GO/Reactome/WikiPathways on "
                        "native chicken gene sets and uses orthologs only for "
                        "human-only collections; 'native_chicken' forces native; "
                        "'ortholog' is the legacy human-mapped route")
    p.add_argument("--chicken-gmt-keying", default=config.CHICKEN_GMT_KEYING,
                   choices=list(config.CHICKEN_GMT_KEYINGS),
                   help="ID space of the native chicken GMT; 'name' (symbols) is "
                        "correct for most exports, since the ENSG-keyed GMT ships "
                        "GRCg7b IDs that will not match a galgal6 table")
    p.add_argument("--strict-1to1", action="store_true",
                   help="keep only unambiguous chicken<->human ortholog pairs "
                        "(ortholog route only)")
    p.add_argument("--min-size", type=int, default=config.GSEA_MIN_SIZE)
    p.add_argument("--max-size", type=int, default=config.GSEA_MAX_SIZE)
    p.add_argument("--static", action="store_true",
                   help="also write SVG/PDF (needs kaleido+chrome)")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)

    custom = None
    if args.custom_gmt:
        from . import genesets
        custom = genesets.load_gmt(args.custom_gmt)

    print(f"[deseq2_enrich] running contrast '{args.name}' ...")
    res = run_contrast(
        args.input, contrast_name=args.name,
        padj_threshold=args.padj, lfc_threshold=args.lfc,
        id_col=args.id_col, rank_metric=args.rank,
        ora_sources=args.sources, gsea_libraries=args.gsea_libs,
        custom_gmt=custom, organism=args.organism,
        do_ora=not args.no_ora, do_gsea=not args.no_gsea,
        gsea_permutations=args.permutations,
        # The CLI is the supported route for reporting-grade runs: 1000
        # permutations needs ~2 GB, which does not fit the deployed app's
        # ~1 GB box. Defaults here match the app so results are comparable.
        gsea_mode=args.gsea_mode,
        chicken_gmt_keying=args.chicken_gmt_keying,
        strict_one_to_one=args.strict_1to1,
        min_size=args.min_size, max_size=args.max_size,
    )
    print(res.report.as_text())
    print("DEG counts:", res.deg_sets.counts)

    base = os.path.join(args.outdir, args.name)

    if res.ora is not None and not res.ora.empty:
        res.ora.to_csv(f"{base}_ORA.csv", index=False)
        _save(plots.ora_dotplot(res.ora), f"{base}_ORA_dotplot", args.static)
        print(f"  wrote {base}_ORA.csv ({len(res.ora)} terms)")
    if res.gsea is not None:
        res.gsea.table.to_csv(f"{base}_GSEA.csv", index=False)
        _save(plots.gsea_bar(res.gsea.table), f"{base}_GSEA_bar", args.static)
        top = res.gsea.table.iloc[0]["term"]
        _save(plots.gsea_running_plot(res.gsea, top), f"{base}_GSEA_running", args.static)
        print(f"  wrote {base}_GSEA.csv ({len(res.gsea.table)} terms)")

    # DE figures
    _save(plots.volcano(res.df, args.padj, args.lfc), f"{base}_volcano", args.static)
    _save(plots.ma_plot(res.df, args.padj, args.lfc), f"{base}_MA", args.static)

    if res.errors:
        print("[warnings]", res.errors)
    print(f"[deseq2_enrich] done -> {args.outdir}/")
    return 0


def _save(fig, base: str, static: bool) -> None:
    """Always write interactive HTML; optionally SVG/PDF if kaleido+chrome."""
    fig.write_html(f"{base}.html", include_plotlyjs="cdn")
    if static:
        try:
            fig.write_image(f"{base}.svg")
            fig.write_image(f"{base}.pdf")
        except Exception as exc:
            print(f"  [static export skipped for {base}: {exc}]")


if __name__ == "__main__":
    raise SystemExit(main())
