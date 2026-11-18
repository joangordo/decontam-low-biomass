# Install the packages used in this project.
#
# Recommended: run renv::restore() from the project root instead. It installs
# the exact versions recorded in renv.lock, so results match the reports.
# This script is the fallback: it installs the latest versions available.

cran <- c("dplyr", "tidyr", "ggplot2", "knitr", "rmarkdown", "scales",
          "testthat", "renv", "BiocManager")
install.packages(setdiff(cran, rownames(installed.packages())))

# decontam is on Bioconductor; phyloseq is only needed for
# data-raw/extract_karstens.R (the processed CSVs are already in data/).
BiocManager::install(c("decontam", "phyloseq"), update = FALSE)
