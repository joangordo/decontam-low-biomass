# End-to-end cleaning workflow -------------------------------------------------
#
# Follows the usual structure of a low-biomass cleaning workflow: each sample type is
# cleaned separately, against its own negative controls, with a read-depth
# rule suited to its biomass. Every decision is returned so it can be
# reported.

#' Clean a low-biomass count table, one sample type at a time.
#'
#' @param counts Samples x taxa count matrix.
#' @param samples Data frame with `sample_id`, `is_neg`, `dna_conc` and the
#'   column named in `group` (negative controls carry the group they belong
#'   to).
#' @param group Column of `samples` defining sample types.
#' @param method decontam method (default "frequency", see Part 3 of the
#'   report). If a group has no negative controls and the method needs them, the
#'   frequency method is used for it and a note is added to the summary.
#' @param threshold decontam threshold: one value, or a named vector with one
#'   value per group (unnamed groups use the first value).
#' @param low_read_rule Named character vector: rule per group ("fixed" or
#'   "neg_median"); unnamed groups use "fixed".
#' @param min_reads Cut-off for the "fixed" rule.
#' @return A list with `counts` (cleaned true samples), `taxa` (per group and
#'   taxon: score and decision), `dropped_samples` and `summary`.
clean_low_biomass <- function(counts, samples, group = "sample_type",
                              method = "frequency", threshold = 0.1,
                              low_read_rule = character(0),
                              min_reads = 1000) {
  validate_counts(counts)
  stopifnot(identical(rownames(counts), samples$sample_id),
            group %in% names(samples))

  out_counts <- list(); out_taxa <- list(); out_drop <- list(); out_sum <- list()

  for (g in unique(samples[[group]])) {
    idx <- samples[[group]] == g
    cg <- counts[idx, , drop = FALSE]
    sg <- samples[idx, , drop = FALSE]
    note <- ""

    rule <- if (g %in% names(low_read_rule)) low_read_rule[[g]] else "fixed"
    if (rule == "neg_median" && !any(sg$is_neg)) {
      rule <- "fixed"
      note <- paste(note, "No controls: fixed read rule used.")
    }
    low <- flag_low_read_samples(cg, sg$is_neg, rule = rule, min_reads = min_reads)
    out_drop[[g]] <- data.frame(group = rep(g, sum(low)),
                                sample_id = names(low)[low],
                                reads = unname(rowSums(cg)[low]))
    cg <- cg[!low, , drop = FALSE]
    sg <- sg[!low, , drop = FALSE]

    m <- method
    if (!any(sg$is_neg) && m != "frequency") {
      m <- "frequency"
      note <- paste(note, "No controls: frequency method used.")
    }
    thr <- if (g %in% names(threshold)) threshold[[g]] else threshold[[1]]
    res <- suppressWarnings(run_decontam(
      cg, conc = sg$dna_conc, is_neg = sg$is_neg, method = m,
      threshold = thr))
    flagged <- stats::setNames(res$contaminant, res$asv_id)

    true_cg <- cg[!sg$is_neg, , drop = FALSE]
    before <- rowSums(true_cg)
    true_cg[, flagged] <- 0
    after <- rowSums(true_cg)

    out_counts[[g]] <- true_cg
    out_taxa[[g]] <- data.frame(group = g, res, row.names = NULL)
    out_sum[[g]] <- data.frame(
      group = g,
      method = m,
      threshold = thr,
      read_rule = paste0(rule, " (cut-off ", round(attr(low, "cutoff")), ")"),
      samples_in = sum(!samples$is_neg[idx]),
      samples_dropped = sum(low),
      negatives = sum(sg$is_neg),
      taxa_scored = sum(res$scored),
      taxa_removed = sum(flagged),
      median_pct_reads_removed = stats::median(100 * (before - after) / before),
      notes = trimws(note)
    )
  }

  bind <- function(x) {
    out <- do.call(rbind, unname(x))
    rownames(out) <- NULL
    out
  }
  list(counts = do.call(rbind, unname(out_counts)),
       taxa = bind(out_taxa),
       dropped_samples = bind(out_drop),
       summary = bind(out_sum))
}
