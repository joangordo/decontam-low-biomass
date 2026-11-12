# Wrappers around decontam::isContaminant() ----------------------------------

#' Run decontam and return a tidy per-taxon result.
#'
#' * frequency: fitted on true samples only, using their DNA concentration.
#'   Negative controls are excluded because the method models how a taxon's
#'   frequency changes with the input DNA of real samples.
#' * prevalence: compares presence/absence in negative controls vs samples.
#' * combined / minimum / either / both: decontam's ways of mixing the two
#'   scores (see ?decontam::isContaminant).
#'
#' Taxa present in fewer than two of the samples used cannot be scored by the
#' frequency method; decontam returns NA for them and never flags them. That
#' matters a lot in low-biomass data, where most contaminants appear in a
#' single sample, so `scored` is returned explicitly.
#'
#' @param counts Samples x taxa count matrix.
#' @param conc DNA concentration per sample (same order as `counts`).
#' @param is_neg Logical, TRUE for negative controls.
#' @param method One of "frequency", "prevalence", "combined", "minimum",
#'   "either", "both".
#' @param threshold Score threshold; a taxon is flagged when its score is
#'   below it. Note that decontam's default is 0.1 (not 0.5).
#' @return Data frame with `asv_id`, `score`, `scored`, `contaminant`.
run_decontam <- function(counts, conc = NULL, is_neg = NULL,
                         method = c("frequency", "prevalence", "combined",
                                    "minimum", "either", "both"),
                         threshold = 0.1) {
  method <- match.arg(method)
  validate_counts(counts)

  needs_conc <- method != "prevalence"
  needs_neg <- method != "frequency"
  if (needs_conc && is.null(conc)) stop("`conc` is required for ", method, ".")
  if (needs_neg && (is.null(is_neg) || !any(is_neg))) {
    stop("`is_neg` with at least one negative control is required for ",
         method, ".")
  }

  if (method == "frequency") {
    use <- if (is.null(is_neg)) rep(TRUE, nrow(counts)) else !is_neg
    use <- use & !is.na(conc)
    res <- decontam::isContaminant(counts[use, , drop = FALSE],
                                   conc = conc[use],
                                   method = "frequency",
                                   threshold = threshold)
  } else {
    use <- if (needs_conc) !is.na(conc) else rep(TRUE, nrow(counts))
    if (needs_conc && any(!use)) {
      message("Dropping ", sum(!use), " sample(s) with missing DNA concentration.")
    }
    res <- decontam::isContaminant(counts[use, , drop = FALSE],
                                   conc = if (needs_conc) conc[use] else NULL,
                                   neg = is_neg[use],
                                   method = method,
                                   threshold = threshold)
  }

  data.frame(
    asv_id = rownames(res),
    score = res$p,
    scored = !is.na(res$p),
    contaminant = res$contaminant & !is.na(res$p),
    row.names = NULL
  )
}
