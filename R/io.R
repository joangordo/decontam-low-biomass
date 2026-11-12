# Data loading and validation ------------------------------------------------

#' Load a count table plus sample and taxonomy metadata, and check that they
#' line up.
#'
#' A common shortcut is to assume that the rows of the metadata and the
#' rows of the OTU table are in the same order, without checking it. If they
#' are not, every downstream result is silently wrong. Here the counts are
#' re-ordered by ID and the function fails loudly if anything is missing.
#'
#' @param dir Folder with counts.csv, samples.csv and taxonomy.csv.
#' @return A list with `counts` (samples x ASVs integer matrix), `samples`
#'   and `taxonomy` (data frames in the same order as the matrix).
load_dataset <- function(dir = "data") {
  counts <- utils::read.csv(file.path(dir, "counts.csv"),
                            row.names = 1, check.names = FALSE)
  samples <- utils::read.csv(file.path(dir, "samples.csv"))
  taxonomy <- utils::read.csv(file.path(dir, "taxonomy.csv"))

  counts <- as.matrix(counts)

  missing_s <- setdiff(samples$sample_id, rownames(counts))
  missing_t <- setdiff(taxonomy$asv_id, colnames(counts))
  if (length(missing_s) > 0) {
    stop("Samples in metadata but not in counts: ",
         paste(missing_s, collapse = ", "))
  }
  if (length(missing_t) > 0) {
    stop("ASVs in taxonomy but not in counts: ",
         paste(utils::head(missing_t), collapse = ", "))
  }

  counts <- counts[samples$sample_id, taxonomy$asv_id, drop = FALSE]

  validate_counts(counts)
  stopifnot(
    identical(rownames(counts), samples$sample_id),
    identical(colnames(counts), taxonomy$asv_id),
    is.logical(samples$is_neg), !anyNA(samples$is_neg)
  )

  list(counts = counts, samples = samples, taxonomy = taxonomy)
}

#' Basic sanity checks on a count matrix (samples in rows, taxa in columns).
validate_counts <- function(counts) {
  if (!is.matrix(counts) || !is.numeric(counts)) {
    stop("`counts` must be a numeric matrix.")
  }
  if (anyNA(counts)) stop("`counts` contains NA values.")
  if (any(counts < 0)) stop("`counts` contains negative values.")
  if (any(counts != round(counts))) stop("`counts` must be integer counts.")
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("`counts` needs sample IDs as row names and taxon IDs as column names.")
  }
  invisible(TRUE)
}
