# -----------------------------------------------------------------------------
# Extract the Karstens et al. (2019) mock-community dilution series into plain
# CSV files, so the rest of the project only depends on {decontam} and the
# tidyverse (phyloseq is only needed for this one-off step).
#
# Run once from the project root:
#   Rscript data-raw/extract_karstens.R
#
# Source: https://github.com/lakarstens/ControllingContaminants16S
# Paper:  Karstens L. et al. (2019) mSystems 4(4):e00290-19.
#         doi:10.1128/mSystems.00290-19  (CC BY 4.0)
# -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(phyloseq))

base_url <- paste0(
  "https://raw.githubusercontent.com/lakarstens/",
  "ControllingContaminants16S/master/Analyses/"
)

env <- new.env()
for (f in c("mockDilutions.RData", "mock_taxa.RData")) {
  tmp <- tempfile(fileext = ".RData")
  download.file(paste0(base_url, f), tmp, mode = "wb", quiet = TRUE)
  load(tmp, envir = env)
}

ps <- env$ps                 # 9 dilutions + 1 extraction blank
expected <- env$mock_taxa    # ASVs belonging to the mock community

stopifnot(!taxa_are_rows(ps), all(expected %in% taxa_names(ps)))

# Drop ASVs with zero reads across all samples
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

counts <- as(otu_table(ps), "matrix")            # samples x ASVs
meta <- data.frame(sample_data(ps), check.names = FALSE)

samples <- data.frame(
  sample_id   = rownames(meta),
  description = meta$SampleDescription,
  dna_conc    = meta$DNA_conc,                     # NanoDrop, ng/uL
  is_neg      = meta$is.neg,
  dilution    = suppressWarnings(              # "Blank" -> NA
                  as.integer(sub("^D", "", rownames(meta)))),
  total_reads = rowSums(counts[rownames(meta), , drop = FALSE])
)
samples <- samples[order(is.na(samples$dilution), samples$dilution), ]

taxonomy <- data.frame(asv_id = taxa_names(ps),
                       as(tax_table(ps), "matrix"),
                       row.names = NULL)
taxonomy$expected <- taxonomy$asv_id %in% expected

counts_df <- data.frame(sample_id = rownames(counts), counts,
                        check.names = FALSE, row.names = NULL)
counts_df <- counts_df[match(samples$sample_id, counts_df$sample_id), ]

dir.create("data", showWarnings = FALSE)
write.csv(counts_df, "data/counts.csv", row.names = FALSE)
write.csv(samples,   "data/samples.csv", row.names = FALSE)
write.csv(taxonomy,  "data/taxonomy.csv", row.names = FALSE)

message(sprintf("Wrote %d samples x %d ASVs (%d expected).",
                nrow(counts), ncol(counts), sum(taxonomy$expected)))
