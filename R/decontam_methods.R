# Wrappers around decontam::isContaminant() ----------------------------------

#' Run decontam and return a tidy per-taxon result.
#'
#' The two scores are computed separately and then combined with decontam's
#' own rules, so every method sees exactly the same evidence:
#'
#' * frequency: fitted on true samples only, using their DNA concentration.
#'   Negative controls are excluded because the method models how a taxon's
#'   frequency changes with the input DNA of real samples (decontam does the
#'   same internally).
#' * prevalence: compares presence/absence in negative controls vs samples,
#'   using every sample, whether or not its DNA concentration is known.
#' * combined / minimum / either / both: decontam's ways of mixing the two
#'   scores (see ?decontam::isContaminant). For "either" and "both" decontam
#'   reports no single score (`p = NA`) and makes the call itself; here the
#'   score is the smaller (either) or larger (both) of the two, which gives
#'   exactly decontam's call when the same threshold is used for both.
#'
#' Taxa present in fewer than two of the samples used cannot be scored by the
#' frequency method; decontam returns NA for them and never flags them. That
#' matters a lot in low-biomass data, where most contaminants appear in a
#' single sample, so `scored` is returned explicitly.
#'
#' decontam needs positive DNA concentrations. Values of zero or below (for
#' example, "below the detection limit") are treated as missing and those
#' samples are left out of the frequency score, with a message. To keep them,
#' replace them before calling, for example by half the detection limit.
#'
#' @param counts Samples x taxa count matrix.
#' @param conc DNA concentration per sample (same order as `counts`).
#' @param is_neg Logical, TRUE for negative controls.
#' @param method One of "frequency", "prevalence", "combined", "minimum",
#'   "either", "both".
#' @param threshold Score threshold (one value); a taxon is flagged when its
#'   score is below it. Note that decontam's default is 0.1 (not 0.5).
#' @return Data frame with `asv_id`, `score`, `score_freq`, `score_prev`,
#'   `scored` and `contaminant`. For "either", `scored` means at least one of
#'   the two scores exists; for every other combination, both must exist.
run_decontam <- function(counts, conc = NULL, is_neg = NULL,
                         method = c("frequency", "prevalence", "combined",
                                    "minimum", "either", "both"),
                         threshold = 0.1) {
  method <- match.arg(method)
  validate_counts(counts)
  if (!is.numeric(threshold) || length(threshold) != 1 || is.na(threshold) ||
      threshold < 0 || threshold > 1) {
    stop("`threshold` must be a single number between 0 and 1.")
  }

  needs_conc <- method != "prevalence"
  needs_neg <- method != "frequency"
  if (needs_conc && is.null(conc)) stop("`conc` is required for ", method, ".")
  if (needs_neg && (is.null(is_neg) || !any(is_neg))) {
    stop("`is_neg` with at least one negative control is required for ",
         method, ".")
  }
  if (is.null(is_neg)) is_neg <- rep(FALSE, nrow(counts))
  if (length(is_neg) != nrow(counts) || anyNA(is_neg)) {
    stop("`is_neg` must be TRUE/FALSE for every sample.")
  }

  n_taxa <- ncol(counts)
  p_freq <- rep(NA_real_, n_taxa)
  p_prev <- rep(NA_real_, n_taxa)

  if (needs_conc) {
    if (length(conc) != nrow(counts)) {
      stop("`conc` must have one value per sample.")
    }
    bad <- !is.na(conc) & conc <= 0
    if (any(bad)) {
      message("Treating ", sum(bad), " DNA concentration(s) of zero or below ",
              "as missing: decontam needs positive values.")
      conc[bad] <- NA
    }
    missing_true <- !is_neg & is.na(conc)
    if (any(missing_true)) {
      message("Frequency score: leaving out ", sum(missing_true),
              " sample(s) without a usable DNA concentration.")
    }
    use <- !is_neg & !is.na(conc)
    if (sum(use) < 2) {
      stop("The frequency score needs at least two true samples with a ",
           "positive DNA concentration.")
    }
    p_freq <- unname(decontam::isContaminant(counts[use, , drop = FALSE],
                                             conc = conc[use],
                                             method = "frequency",
                                             threshold = threshold)$p)
  }

  if (needs_neg) {
    p_prev <- unname(decontam::isContaminant(counts, neg = is_neg,
                                             method = "prevalence",
                                             threshold = threshold)$p)
  }

  # Same formulas as decontam::isContaminant(); NA never counts as a call.
  score <- switch(
    method,
    frequency = p_freq,
    prevalence = p_prev,
    minimum = pmin(p_freq, p_prev),
    combined = stats::pchisq(-2 * log(p_freq * p_prev), df = 4,
                             lower.tail = FALSE),
    either = pmin(p_freq, p_prev, na.rm = TRUE),
    both = pmax(p_freq, p_prev)
  )

  data.frame(
    asv_id = colnames(counts),
    score = score,
    score_freq = p_freq,
    score_prev = p_prev,
    scored = !is.na(score),
    contaminant = !is.na(score) & score < threshold,
    row.names = NULL
  )
}
