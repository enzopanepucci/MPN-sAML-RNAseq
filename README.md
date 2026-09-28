# MPN–sAML paired transcriptomic analysis

This repository contains a reproducible R Markdown workflow for paired bulk
RNA-seq analysis of disease progression from chronic-phase MPN to secondary AML
(sAML).

## Analysis overview

- Differential expression is performed with DESeq2 using the 12 paired
  progressor patients only (24 samples).
- The paired model is `~ Patient_ID + Condition`, with the post-progression
  sample contrasted against the pre-progression sample.
- Genes are retained when they have at least 10 counts in at least six
  progressor samples.
- X/Y genes are excluded, while genes without a gene-symbol annotation remain
  in the differential-expression results.
- GSEA is performed for KEGG and Hallmark gene sets using the documented
  hybrid ranking strategy.
- CIBERSORTx cell-state proportions are summarized descriptively and tested
  within paired progressor samples.
- Non-progressor samples are used only for descriptive visualizations and do
  not enter the paired differential-expression contrast.
- `Amostra_18` is retained in the supplied count matrix when present, but is
  excluded from the analyzed sample set because it has no corresponding
  metadata.

## Project structure

`analysis.Rmd` is the single source of truth. It reads the files in `inputs/`,
creates the DESeq2, GSEA and CIBERSORTx objects, and writes figures, tables and
checkpoints to `outputs/` and `checkpoints/`.

The only auxiliary script is `scripts/install_packages.R`. It installs missing
CRAN and Bioconductor dependencies when needed; it does not contain analysis
logic.

## Running the workflow

Open `MPN_sAML.Rproj` in RStudio and render `analysis.Rmd`. The setup chunk
checks for required packages and calls `scripts/install_packages.R` only when a
dependency is missing.

To force recalculation instead of loading saved checkpoints, set
`reuse_checkpoints: false` in `config.yml` and render the document again.

The generated files include publication-ready PNG/PDF figures, supplementary
figures, CSV tables, an HTML report, package/session information and input
checksums. The bundled DESeq2/GSEA checkpoints and outputs are provided for
inspection. Rendering `analysis.Rmd` creates a fresh HTML report and refreshes
the plot-object checkpoint with the current figure styling. `outputs/` and
`checkpoints/` are ignored by Git so generated files are not committed
automatically.

## Inputs

The workflow expects the supplied count matrix, sample metadata, CIBERSORTx
results, gene-module definition and frozen annotation/gene-set files under
`inputs/`.
