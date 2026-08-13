# DESeq2 → ORA / GSEA enrichment (chicken-first, pure Python)

A pipeline **and** a Streamlit app that take a DESeq2 differential-expression
table and run over-representation analysis (ORA) and pre-ranked GSEA, then plot
the results. Built pure-Python so it deploys to **Streamlit Community Cloud**
with no R/Bioconductor.

**Native chicken end to end.** Both ORA and GSEA run on *Gallus gallus*
annotations by default, so no gene is dropped for lacking a human counterpart
and the ranked list stays your real gene universe. Human orthologs are used
only for the two collections that have no chicken equivalent (MSigDB Hallmark
and Oncogenic Signatures), and only if you select them — in which case they run
as a separate prerank with their own FDR and a visible mapping rate.

---

## What it does

| Stage | Engine | Notes |
|-------|--------|-------|
| Load & validate | pandas | auto-detects DESeq2 columns; reports ID coverage, biotypes, NA-padj |
| DEG selection | — | default `padj < 0.05` with no fold-change cutoff; up/down run separately; **universe = tested genes** |
| **ORA** | **g:Profiler** (`ggallus`) | native chicken GO/KEGG/Reactome/WikiPathways; **custom background** |
| Ortholog map | g:Profiler `orth` | chicken → human; **only** for Hallmark / Oncogenic, keyed on stable Ensembl IDs, with a mapping-rate guardrail |
| Ranking | — | DESeq2 **Wald `stat`** (default), `sign(FC)·-log10p`, or `log2FC`; duplicates collapsed by max magnitude |
| **GSEA** | **gseapy** prerank | **native chicken GO / Reactome / WikiPathways** (g:Profiler GMT) + your custom `.gmt`; Hallmark / Oncogenic via orthologs on request |
| Plots | Plotly + matplotlib | volcano, MA, ORA dotplot, GSEA bar + running-ES, leading-edge heatmap, enrichment map, UpSet |

Two design choices are enforced because they are the usual ORA pitfalls: the
background is the **tested** gene set (non-NA `padj`), not the genome; and up-
and down-regulated genes are analysed **separately**.

---

## Expected input

Any DESeq2 results table (TSV/CSV). Columns are auto-detected from common
aliases. Required (after mapping): `gene_id`, `log2FoldChange`, `pvalue`,
`padj`. Recommended: `entrez_id`, `gene_name`, `stat`, `baseMean`,
`gene_biotype`. If the gene ID sits in an unnamed first column it is picked up
automatically. `padj = NA` rows (independent filtering) are dropped from the
ORA universe as they should be.

DEG threshold defaults. The pipeline defaults to `padj < 0.05` with no
fold-change cutoff. A fold-change cutoff is a filter, not a statistical test;
apply it only when you have a specific biological reason. Adjust `|log2FC|` in
the sidebar for your contrast.

The bundled demo table (`sample_data/DESeq2_chicken_demo.tsv`) is a real
DESeq2 export from GEO accession GSE230804. It compares HH stage 18 chick
embryo hindbrain boundary cells sorted as CSPG-positive (4 biological
replicates) versus CSPG-negative (6 biological replicates), using the
GEO-supplied DESeq2 workbook.

---

## Run locally

```bash
pip install -e .
streamlit run app/streamlit_app.py
```

Toggle **Use bundled sample data** in the sidebar, or upload one file per
contrast, set thresholds, and click **Run enrichment**.

---

## Deploy to Streamlit Community Cloud

1. Push this folder to a **public GitHub repo** (see structure below).
2. On <https://share.streamlit.io> → **New app** → pick the repo.
3. Set **Main file path** to `app/streamlit_app.py`.
4. In **Advanced settings**, choose Python **3.12**.
5. Deploy. `requirements.txt` is picked up automatically.

Notes:
* Nothing licensed is committed — MSigDB/Reactome gene sets are fetched at
  runtime from Enrichr-hosted libraries, so the public repo stays clean.
* No `packages.txt` is needed for the Streamlit app. Plotly figures render
  interactively in the browser, which keeps Community Cloud startup much
  faster and avoids installing Chromium.
* The box is ~1 GB RAM / 1 CPU. Keep **Quick mode** on for previews.
* If the logs show Python 3.14, switch the app's Python version back to 3.12.
  The scientific stack installs much more predictably on 3.12.

---

## Batch use on HPC / Tempest (the CLI)

Same engine, headless, for the real multi-contrast design:

```bash
python -m deseq2_enrich.cli \
  --input results/Sacral_vs_Caudal.tsv \
  --name Sacral_vs_Caudal \
  --outdir out/ \
  --sources GO:BP KEGG REAC WP \
  --gsea-libs MSigDB_Hallmark_2020 Reactome_2022 \
  --rank stat --padj 0.05 --lfc 1.0
```

Writes `<name>_ORA.csv`, `<name>_GSEA.csv`, and figures (HTML always; SVG/PDF
if you additionally install `kaleido` plus Chrome/Chromium and pass
`--static`). Loop over contrasts in a Slurm array; the g:Profiler/Enrichr calls
are memoised within a process.

---

## Custom gene-set modules

Upload a `.gmt` in the sidebar (or `--custom-gmt`) to score your curated
modules (cGAS-STING, necroptosis, apoptosis, inflammasome, osteoclast,
osteoblast) through the identical GSEA path. Format:
`module_name<TAB>description<TAB>GENE1<TAB>GENE2…`. Genes must match the ranking
namespace — under the default native route that is **chicken gene symbols**
(the ranking is built from `gene_name`/`gene_id`), not human symbols. Case is
normalised on both sides, so casing will not silently zero the overlap.

---

## Adapting to another organism

Change `ORGANISM` (and, if needed, `ORTHOLOG_TARGET`) in
`deseq2_enrich/config.py` to any g:Profiler code (`hsapiens`, `mmusculus`, …).
For a species already covered by MSigDB you can skip the ortholog step by
ranking on symbols directly. The native-chicken GMT route is chicken-specific:
`GPROFILER_GMT_URLS` points at g:Profiler's per-organism files, so another
organism needs its own URL and a check that the GMT's ID space matches your
table (see the note on assembly-version mismatch in `config.py`).

---

## Project structure

```
deseq2_enrich/
├── app/streamlit_app.py        # the app (Main file path on Cloud)
├── deseq2_enrich/              # the pipeline package
│   ├── config.py  io.py  degs.py  rank.py
│   ├── ortho.py   ora.py  genesets.py  gsea.py
│   ├── plots.py   pipeline.py  cli.py
├── sample_data/                # bundled demo table
├── requirements.txt  .streamlit/config.toml
└── README.md
```

## Caveats

* ORA and the MSigDB library fetch require internet (both keyless services).
* Non-protein-coding genes (lncRNA) are largely unannotated in GO/KEGG and drop
  out of enrichment; the coverage report tells you how many.
* GSEA significance is permutation-based — use the full 1000 permutations for
  anything you report, not Quick mode.
* g:Profiler `g:SCS` is the multiple-testing correction for ORA; GSEA uses
  gseapy's FDR q-value. Don't compare the two thresholds directly.

### Large collections: chunked GSEA and its FDR

Native GO is ~12.7k chicken terms (~4.4k scoreable). Scoring that in one
`gseapy` call peaks at **~3.7 GB** — it cannot run on Streamlit Community
Cloud's ~1 GB. Collections above `GSEA_CHUNK_THRESHOLD` (800 sets) are
therefore scored in blocks of `GSEA_CHUNK_SIZE` (100), freeing each block
before the next.

Measured on a 16.3k-gene table, full GO:BP, 4,372 sets scored:

| permutations | whole run | chunked |
|---|---|---|
| 100 | ~3.8 GB | **~750–880 MB** (~40 s) |
| 1000 | ~3.7 GB | ~1.8 GB (~5 min) |

So **GO:BP now runs in-app at 100 permutations**. Peak memory scales with
permutation count, so 1000-permutation GO still belongs in a local CLI run —
the app says so when you ask for it.

Two things to know about the numbers this produces:

* **ES, NES and nominal p are unchanged** by chunking — bit-identical to a
  whole run on the same seed and ranking, because gseapy seeds its
  permutations per gene set. There is a regression test for this.
* **FDR is different.** A chunked run reports **Benjamini–Hochberg over the
  pooled permutation p-values** (the `fdr_method` column and the export
  manifest say so), not gseapy's single-run NES-based q-value, which cannot be
  reconstructed without pooling each set's null-NES distribution across blocks.
  BH-adjusted permutation p is the same choice `fgsea`'s `padj` makes. The
  nominal p is floored at `1/(nperm+1)` *before* adjustment, so a term cannot
  report a q below what the multiple-testing burden allows.
