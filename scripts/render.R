# Execute a partir da raiz do projeto (abra MPN_sAML.Rproj no RStudio).
source("scripts/install_packages.R")
dir.create("outputs/reports", recursive = TRUE, showWarnings = FALSE)
# analysis.Rmd é a fonte única; analysis.R é uma cópia gerada para uso no Console.
knitr::purl("analysis.Rmd", output = "analysis.R", documentation = 1, quiet = TRUE)
rmarkdown::render("analysis.Rmd", output_dir = "outputs/reports",
                  envir = globalenv(), clean = TRUE)
