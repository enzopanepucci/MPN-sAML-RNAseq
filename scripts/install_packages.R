# Rode uma vez. Instala apenas pacotes ausentes, sem atualizar todo o ambiente.
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
cran <- c("BiocManager", "dplyr", "tidyr", "tibble", "ggplot2", "ggrepel",
          "pheatmap", "cowplot", "patchwork", "scales", "Rtsne", "yaml",
          "digest", "data.table", "matrixStats", "rmarkdown", "knitr", "renv")
missing <- setdiff(cran, rownames(installed.packages()))
if (length(missing)) install.packages(missing)
bioc <- c("DESeq2", "fgsea")
missing <- setdiff(bioc, rownames(installed.packages()))
if (length(missing)) BiocManager::install(missing, ask = FALSE, update = FALSE)
message("Pacotes necessários disponíveis.")
