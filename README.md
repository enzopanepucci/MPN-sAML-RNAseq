# MPN–sAML paired transcriptomic analysis

This repository contains a reproducible R workflow for paired bulk RNA-seq analysis of disease progression from chronic-phase myeloproliferative neoplasms (MPN) to secondary acute myeloid leukaemia (sAML).

## Analysis overview

- Differential expression is performed with DESeq2 using the 12 paired progressor patients only (24 samples).
- The paired model is `~ Patient_ID + Condition`, with the post-progression sample contrasted against the pre-progression sample.
- Genes are retained when they have at least 10 counts in at least six progressor samples.
- Genes annotated on chromosomes X/Y are excluded. Genes without a gene-symbol annotation remain in the differential-expression analysis.
- Differentially expressed genes are selected using adjusted *P* < 0.05 and |log2 fold change| > 0.5.
- Preranked GSEA uses MSigDB 7.5.1 KEGG and Hallmark gene sets and the documented hybrid ranking strategy.
- Precomputed CIBERSORTx cell-state proportions are summarized across groups and tested within paired progressor samples.
- First and second samples from five non-progressor patients are included as descriptive comparators in the heatmap and deconvolution analyses; they do not enter differential expression or GSEA.
- `Amostra_18` has no corresponding metadata and is excluded from the analyzed sample set. Its count summary is retained in sample quality control.

## Project structure

`analysis.R` is the only executable script. It contains the analytical parameters, runs the analysis, creates plotting objects and exports all figures and tables. Detailed procedures, including the exact ranking formula, are documented in [`SUPPLEMENTARY_METHODS.md`](SUPPLEMENTARY_METHODS.md).

| File or directory | Purpose |
| --- | --- |
| `analysis.R` | Complete analysis and figure/table exports. |
| `environment.yml` | Exported conda environment with fixed dependency versions and builds for Linux x86_64. |
| `MPN_sAML.Rproj` | RStudio project file, if using RStudio. |
| `SUPPLEMENTARY_METHODS.md` | Detailed analytical methods and descriptions of supplementary results. |
| `outputs/figures/` | Figures 1A–D and 2A–B in PDF and PNG formats. |
| `outputs/supplementary_figures/` | Figures S1A–B, S2, S3, S4A–B and S5 in PDF and PNG formats. |
| `outputs/tables/` | Tables S1–S11, filtering summaries, gene rankings and PCA/t-SNE coordinates. |
| `outputs/reports/` | Session information, installed package versions, parameters used, input checksums and figure dimensions. |

Generated outputs are excluded from version control through `.gitignore`.

## Software and citations

The workflow was validated using **R 4.5.3**. The following package versions match `environment.yml`; that file also fixes transitive dependencies and conda builds.

| R package | Version |
| --- | --- |
| DESeq2 | 1.50.2 |
| fgsea | 1.36.2 |
| dplyr | 1.2.1 |
| tidyr | 1.3.2 |
| tibble | 3.3.1 |
| ggplot2 | 4.0.3 |
| ggrepel | 0.9.8 |
| pheatmap | 1.0.13 |
| cowplot | 1.2.0 |
| patchwork | 1.3.2 |
| scales | 1.4.0 |
| Rtsne | 0.17 |
| yaml | 2.3.12 |
| digest | 0.6.39 |
| data.table | 1.18.6.1 |
| matrixStats | 1.5.0 |

The cell-state proportions were calculated with **CIBERSORTx**. The R workflow reads its precomputed results and does not rerun deconvolution.

### Analysis and visualization references

- **DESeq2:** Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15, 550. [doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8).
- **fgsea:** Korotkevich G, Sukhov V, Sergushichev A (2019). Fast gene set enrichment analysis. *bioRxiv*. [doi:10.1101/060012](https://doi.org/10.1101/060012).
- **ggplot2:** Wickham H (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag New York. ISBN 978-3-319-24277-4. [Package documentation](https://ggplot2.tidyverse.org).
- **Rtsne:** Krijthe JH (2015). *Rtsne: T-Distributed Stochastic Neighbor Embedding using Barnes-Hut Implementation*. [Package repository](https://github.com/jkrijthe/Rtsne).
- **t-SNE:** van der Maaten L, Hinton G (2008). Visualizing High-Dimensional Data Using t-SNE. *Journal of Machine Learning Research* 9, 2579–2605.
- **Barnes-Hut t-SNE:** van der Maaten L (2014). Accelerating t-SNE using Tree-Based Algorithms. *Journal of Machine Learning Research* 15, 3221–3245.
- **CIBERSORTx:** Newman AM, Steen CB, Liu CL, et al. (2019). Determining cell type abundance and expression from bulk tissues with digital cytometry. *Nature Biotechnology* 37, 773–782. [doi:10.1038/s41587-019-0114-2](https://doi.org/10.1038/s41587-019-0114-2).

### Additional R software references

Version numbers are listed above. Package documentation and citation information are available through the following links:

- **dplyr:** Wickham H, François R, Henry L, Müller K, Vaughan D. *dplyr: A Grammar of Data Manipulation*. [CRAN](https://CRAN.R-project.org/package=dplyr).
- **tidyr:** Wickham H, Vaughan D, Girlich M. *tidyr: Tidy Messy Data*. [CRAN](https://CRAN.R-project.org/package=tidyr).
- **tibble:** Müller K, Wickham H. *tibble: Simple Data Frames*. [CRAN](https://CRAN.R-project.org/package=tibble).
- **ggrepel:** Slowikowski K. *ggrepel: Automatically Position Non-Overlapping Text Labels with 'ggplot2'*. [CRAN](https://CRAN.R-project.org/package=ggrepel).
- **pheatmap:** Kolde R. *pheatmap: Pretty Heatmaps*. [CRAN](https://CRAN.R-project.org/package=pheatmap).
- **cowplot:** Wilke C. *cowplot: Streamlined Plot Theme and Plot Annotations for 'ggplot2'*. [CRAN](https://CRAN.R-project.org/package=cowplot).
- **patchwork:** Pedersen T. *patchwork: The Composer of Plots*. [CRAN](https://CRAN.R-project.org/package=patchwork).
- **scales:** Wickham H, Pedersen T, Seidel D. *scales: Scale Functions for Visualization*. [CRAN](https://CRAN.R-project.org/package=scales).
- **yaml:** Garbett SP, Stephens J, Simonov K, Xie Y, Dong Z, Wickham H, Horner J, reikoch, Beasley W, O'Connor B, Warnes GR, Quinn M, Kamvar ZN, Gao C. *yaml: Methods to Convert R Data to YAML and Back*. [CRAN](https://CRAN.R-project.org/package=yaml).
- **digest:** Eddelbuettel D. *digest: Create Compact Hash Digests of R Objects*. [CRAN](https://CRAN.R-project.org/package=digest).
- **data.table:** Barrett T, Dowle M, Srinivasan A, Gorecki J, Chirico M, Hocking T, Schwendinger B, Krylov I. *data.table: Extension of `data.frame`*. [CRAN](https://CRAN.R-project.org/package=data.table).
- **matrixStats:** Bengtsson H. *matrixStats: Functions that Apply to Rows and Columns of Matrices (and to Vectors)*. [CRAN](https://CRAN.R-project.org/package=matrixStats).

## Data availability

The RNA-sequencing data re-analysed in this study were generated by de Castro et al. and are deposited under restricted access at [Zenodo](https://doi.org/10.5281/zenodo.1212754). Access is governed by the source study's data-sharing conditions. The analysis code and shareable derived results are available at [MPN-sAML-RNAseq](https://github.com/enzopanepucci/MPN-sAML-RNAseq).
