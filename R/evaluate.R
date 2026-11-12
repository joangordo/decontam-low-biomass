# Evaluation against the mock-community ground truth --------------------------
#
# Convention: "positive" = contaminant (any ASV not in the mock community),
# so sensitivity = share of contaminants removed and specificity = share of
# genuine mock ASVs kept.

#' Confusion-matrix metrics at the ASV level.
#'
#' @param flagged Named logical vector, TRUE = removed by the method.
#' @param expected Named logical vector, TRUE = genuine mock-community ASV.
classification_metrics <- function(flagged, expected) {
  stopifnot(identical(names(flagged), names(expected)), !anyNA(flagged))
  contaminant <- !expected
  tp <- sum(flagged & contaminant)
  fp <- sum(flagged & !contaminant)
  fn <- sum(!flagged & contaminant)
  tn <- sum(!flagged & !contaminant)
  data.frame(
    tp = tp, fp = fp, fn = fn, tn = tn,
    sensitivity = safe_div(tp, tp + fn),
    specificity = safe_div(tn, tn + fp),
    precision   = safe_div(tp, tp + fp)
  )
}

#' Read-level view: per sample, the share of reads that are contaminants
#' before and after removal, and the share of genuine reads lost.
#'
#' ASV-level and read-level performance can disagree sharply: many
#' contaminant ASVs are rare singletons, while a few abundant ones carry most
#' of the contaminant reads.
read_level_metrics <- function(counts, flagged, expected) {
  stopifnot(identical(colnames(counts), names(flagged)),
            identical(names(flagged), names(expected)))
  total <- rowSums(counts)
  contam_reads <- rowSums(counts[, !expected, drop = FALSE])
  mock_reads <- rowSums(counts[, expected, drop = FALSE])

  kept <- counts[, !flagged, drop = FALSE]
  kept_total <- rowSums(kept)
  kept_contam <- rowSums(kept[, !expected[!flagged], drop = FALSE])
  removed_mock <- rowSums(counts[, flagged & expected, drop = FALSE])

  data.frame(
    sample_id = rownames(counts),
    pct_contam_before = 100 * safe_div(contam_reads, total),
    pct_contam_after = 100 * safe_div(kept_contam, kept_total),
    pct_contam_reads_removed = 100 * safe_div(contam_reads - kept_contam,
                                              contam_reads),
    pct_mock_reads_lost = 100 * safe_div(removed_mock, mock_reads),
    row.names = NULL
  )
}

#' Rank-based AUC (probability that a random contaminant gets a lower decontam
#' score than a random genuine ASV). Unscored taxa (NA) are dropped, so report
#' how many were scored alongside it.
score_auc <- function(score, expected) {
  ok <- !is.na(score)
  s <- score[ok]
  pos <- !expected[ok]                 # contaminants
  n_pos <- sum(pos)
  n_neg <- sum(!pos)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  # Lower score = more contaminant-like, so rank on -score
  r <- rank(-s)
  (sum(r[pos]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

safe_div <- function(num, den) ifelse(den == 0, NA_real_, num / den)
