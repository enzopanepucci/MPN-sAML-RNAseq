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
- `Amostra_18` has no corresponding metadata and is excluded from the
  analyzed sample set.

## Project structure

`analysis.Rmd` is the single source of truth. It creates the DESeq2, GSEA and
CIBERSORTx plotting objects, and writes figures, tables and checkpoints to
`outputs/` and `checkpoints/`. Analytical details, including the exact ranking
formula, are in [`SUPPLEMENTARY_METHODS.md`](SUPPLEMENTARY_METHODS.md).

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
figures, CSV tables, an HTML report, package/session information and checksums.
Rendering `analysis.Rmd` creates a fresh HTML report and refreshes the plot-object
checkpoint with the current figure styling. Generated results and checkpoints
are ignored by Git.

## Software and citations

The workflow uses **R 4.3.3** and these directly loaded packages (versions from
the recorded analysis session): DESeq2 1.42.0, fgsea 1.28.0, dplyr 1.1.4,
tidyr 1.3.1, tibble 3.2.1, ggplot2 3.4.4, ggrepel 0.9.5, pheatmap 1.0.12,
cowplot 1.1.3, patchwork 1.2.0, scales 1.3.0, Rtsne 0.17, yaml 2.3.8,
digest 0.6.34, data.table 1.14.10, matrixStats 1.2.0, knitr 1.45 and
rmarkdown 2.25. The cell-state proportions were calculated with **CIBERSORTx**;
the R workflow reads its precomputed results and does not rerun deconvolution.

The references below are the output of `citation()` supplied by the analyst.
Package-version numbers in these references reflect the installation used to
generate the citations; several differ from the versions in the recorded
analysis session above.

### R package references

- **DESeq2:** Love MI, Huber W, Anders S (2014). “Moderated estimation of fold
  change and dispersion for RNA-seq data with DESeq2.” *Genome Biology* 15,
  550. [doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8).
- **fgsea:** Korotkevich G, Sukhov V, Sergushichev A (2019). “Fast gene set
  enrichment analysis.” *bioRxiv*.
  [doi:10.1101/060012](https://doi.org/10.1101/060012).
- **dplyr:** Wickham H, François R, Henry L, Müller K, Vaughan D (2026).
  *dplyr: A Grammar of Data Manipulation*. R package version 1.2.1.
  <https://CRAN.R-project.org/package=dplyr>.
- **tidyr:** Wickham H, Vaughan D, Girlich M (2025). *tidyr: Tidy Messy Data*.
  R package version 1.3.2. <https://CRAN.R-project.org/package=tidyr>.
- **tibble:** Müller K, Wickham H (2023). *tibble: Simple Data Frames*.
  R package version 3.2.1. <https://CRAN.R-project.org/package=tibble>.
- **ggplot2:** Wickham H (2016). *ggplot2: Elegant Graphics for Data Analysis*.
  Springer-Verlag New York. ISBN 978-3-319-24277-4.
  <https://ggplot2.tidyverse.org>.
- **ggrepel:** Slowikowski K (2026). *ggrepel: Automatically Position
  Non-Overlapping Text Labels with 'ggplot2'*. R package version 0.9.8.
  <https://CRAN.R-project.org/package=ggrepel>.
- **pheatmap:** Kolde R (2019). *pheatmap: Pretty Heatmaps*. R package version
  1.0.12. <https://CRAN.R-project.org/package=pheatmap>.
- **cowplot:** Wilke C (2025). *cowplot: Streamlined Plot Theme and Plot
  Annotations for 'ggplot2'*. R package version 1.2.0.
  <https://CRAN.R-project.org/package=cowplot>.
- **patchwork:** Pedersen T (2025). *patchwork: The Composer of Plots*. R package
  version 1.3.2. <https://CRAN.R-project.org/package=patchwork>.
- **scales:** Wickham H, Pedersen T, Seidel D (2025). *scales: Scale Functions
  for Visualization*. R package version 1.4.0.
  <https://CRAN.R-project.org/package=scales>.
- **Rtsne:** Krijthe JH (2015). *Rtsne: T-Distributed Stochastic Neighbor
  Embedding using Barnes-Hut Implementation*. R package version 0.17.
  <https://github.com/jkrijthe/Rtsne>.
  - van der Maaten L, Hinton G (2008). “Visualizing High-Dimensional Data Using
    t-SNE.” *Journal of Machine Learning Research* 9, 2579–2605.
  - van der Maaten L (2014). “Accelerating t-SNE using Tree-Based Algorithms.”
    *Journal of Machine Learning Research* 15, 3221–3245.
- **yaml:** Garbett SP, Stephens J, Simonov K, Xie Y, Dong Z, Wickham H,
  Horner J, reikoch, Beasley W, O'Connor B, Warnes GR, Quinn M, Kamvar ZN,
  Gao C (2024). *yaml: Methods to Convert R Data to YAML and Back*. R package
  version 2.3.10. <https://CRAN.R-project.org/package=yaml>.
- **digest:** Eddelbuettel D (2024). *digest: Create Compact Hash Digests of R
  Objects*. R package version 0.6.37.
  <https://CRAN.R-project.org/package=digest>.
- **data.table:** Barrett T, Dowle M, Srinivasan A, Gorecki J, Chirico M,
  Hocking T, Schwendinger B, Krylov I (2026). *data.table: Extension of
  `data.frame`*. R package version 1.18.4.
  <https://CRAN.R-project.org/package=data.table>.
- **matrixStats:** Bengtsson H (2025). *matrixStats: Functions that Apply to
  Rows and Columns of Matrices (and to Vectors)*. R package version 1.5.0.
  <https://CRAN.R-project.org/package=matrixStats>.
- **knitr:** Xie Y (2025). *knitr: A General-Purpose Package for Dynamic Report
  Generation in R*. R package version 1.51. <https://yihui.org/knitr/>.
  - Xie Y (2015). *Dynamic Documents with R and knitr*, 2nd edition. Chapman
    and Hall/CRC, Boca Raton, Florida. ISBN 978-1498716963.
    <https://yihui.org/knitr/>.
  - Xie Y (2014). “knitr: A Comprehensive Tool for Reproducible Research in R.”
    In Stodden V, Leisch F, Peng RD (eds.), *Implementing Reproducible
    Computational Research*. Chapman and Hall/CRC. ISBN 978-1466561595.
- **rmarkdown:** Allaire J, Xie Y, Dervieux C, McPherson J, Luraschi J, Ushey K,
  Atkins A, Wickham H, Cheng J, Chang W, Iannone R (2026). *rmarkdown: Dynamic
  Documents for R*. R package version 2.31.
  <https://github.com/rstudio/rmarkdown>.
  - Xie Y, Allaire J, Grolemund G (2018). *R Markdown: The Definitive Guide*.
    Chapman and Hall/CRC, Boca Raton, Florida. ISBN 9781138359338.
    <https://yihui.org/rmarkdown/>.
  - Xie Y, Dervieux C, Riederer E (2020). *R Markdown Cookbook*. Chapman and
    Hall/CRC, Boca Raton, Florida. ISBN 9780367563837.
    <https://yihui.org/rmarkdown-cookbook>.

**CIBERSORTx:** Newman AM, Steen CB, Liu CL, et al. Determining cell type
abundance and expression from bulk tissues with digital cytometry.
*Nature Biotechnology*. 2019;37:773–782.
[doi:10.1038/s41587-019-0114-2](https://doi.org/10.1038/s41587-019-0114-2).

The original RNA-seq data were deposited under restricted access by de Castro
et al. at [Zenodo](https://doi.org/10.5281/zenodo.1212754); access is governed
by the source study's data-sharing conditions.
