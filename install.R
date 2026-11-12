# Install the packages used in this project.
# Afterwards, run renv::init() once to pin exact versions in renv.lock.

cran <- c("dplyr", "tidyr", "ggplot2", "knitr", "rmarkdown", "scales",
          "testthat", "renv", "BiocManager")
install.packages(setdiff(cran, rownames(installed.packages())))

# decontam is on Bioconductor; phyloseq is only needed for
# data-raw/extract_karstens.R (the processed CSVs are already in data/).
BiocManager::install(c("decontam", "phyloseq"), update = FALSE)
