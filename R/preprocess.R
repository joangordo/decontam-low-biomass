# Pre-processing and simple (non-decontam) filters ----------------------------
#
# All `flag_*()` functions return a named logical vector over taxa where TRUE
# means "remove as a likely contaminant", so every method can be evaluated
# with the same code in R/evaluate.R.

#' Row-normalise counts to relative abundances. Samples with zero reads get 0
#' instead of NaN.
relative_abundance <- function(counts) {
  validate_counts(counts)
  totals <- rowSums(counts)
  rel <- counts / ifelse(totals == 0, 1, totals)
  rel
}

#' Abundance pre-filter: flag taxa that never reach `min_rel` relative
#' abundance in at least `min_samples` *true* samples.
#'
#' Computing this over all samples, controls included, would keep a
#' taxon that is abundant only in a negative control. Here controls
#' are excluded from the calculation.
#'
#' Default: >1% in at least one sample.
flag_low_abundance <- function(counts, is_neg, min_rel = 0.01, min_samples = 1) {
  stopifnot(length(is_neg) == nrow(counts))
  rel <- relative_abundance(counts)[!is_neg, , drop = FALSE]
  keep <- colSums(rel > min_rel) >= min_samples
  stats::setNames(!keep, colnames(counts))
}

#' Naive filter: flag every taxon with at least one read in any negative
#' control. Known to remove genuine taxa through cross-contamination
#' (well-to-well leakage puts real sample DNA into the blanks).
flag_in_negative <- function(counts, is_neg) {
  stopifnot(length(is_neg) == nrow(counts), any(is_neg))
  flag <- colSums(counts[is_neg, , drop = FALSE]) > 0
  stats::setNames(flag, colnames(counts))
}

#' Flag taxa that dominate the controls: higher prevalence, higher maximum
#' and higher mean relative abundance in controls than in true samples.
#'
#' A middle ground between ignoring the controls and removing everything seen
#' in them: a taxon is flagged only if the controls dominate it on all three
#' counts.
flag_control_dominant <- function(counts, is_neg) {
  stopifnot(length(is_neg) == nrow(counts), any(is_neg), any(!is_neg))
  rel <- relative_abundance(counts)
  ctrl <- rel[is_neg, , drop = FALSE]
  samp <- rel[!is_neg, , drop = FALSE]

  prev_ctrl <- colMeans(ctrl > 0)
  prev_samp <- colMeans(samp > 0)
  max_ctrl <- apply(ctrl, 2, max)
  max_samp <- apply(samp, 2, max)
  mean_ctrl <- colMeans(ctrl)
  mean_samp <- colMeans(samp)

  flag <- prev_ctrl > prev_samp & max_ctrl > max_samp & mean_ctrl > mean_samp
  stats::setNames(flag, colnames(counts))
}

#' Flag samples with too few reads to be trusted.
#'
#' Two common rules:
#' * "fixed": fewer than `min_reads` reads (suggested 1000 for low-biomass
#'   samples such as biopsies and swabs).
#' * "neg_median": fewer reads than the median of the negative controls
#'   (suggested for very-low-biomass samples such as blood or ascitic fluid,
#'   where a fixed cut-off would discard most samples).
#'
#' Negative controls themselves are never flagged.
#' @return Named logical vector over samples, TRUE = drop.
flag_low_read_samples <- function(counts, is_neg,
                                  rule = c("fixed", "neg_median"),
                                  min_reads = 1000) {
  rule <- match.arg(rule)
  stopifnot(length(is_neg) == nrow(counts))
  totals <- rowSums(counts)
  cutoff <- switch(rule,
    fixed = min_reads,
    neg_median = {
      if (!any(is_neg)) stop("rule = 'neg_median' needs negative controls.")
      stats::median(totals[is_neg])
    }
  )
  flag <- totals < cutoff & !is_neg
  attr(flag, "cutoff") <- cutoff
  stats::setNames(flag, rownames(counts))
}
