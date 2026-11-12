# Simulating low-biomass 16S studies with a known truth -----------------------
#
# Generative model (per study):
#
# * A "true" community of `n_true` taxa with skewed abundances. Each subject
#   carries each taxon with a taxon-specific probability, and relative
#   abundances vary between subjects (gamma/Dirichlet noise), so genuine taxa
#   are not present everywhere, as in real clinical samples.
# * A reagent "kitome" of `n_contam` contaminant taxa contributing roughly the
#   same absolute amount of DNA to every tube (`contam_dna`, with small
#   tube-to-tube variation). Each contaminant is sporadic: it reaches a given
#   tube with a taxon-specific probability.
# * Each true sample has its own biomass B ~ logNormal(log(biomass_ratio *
#   contam_dna), biomass_sdlog). The ratio of biomass to contaminant DNA sets
#   how contaminated the study is.
# * Measured DNA concentration = (B + C) x logNormal(0, conc_noise): small
#   `conc_noise` mimics precise fluorometric/qPCR quantification, large values
#   mimic imprecise spectrophotometry.
# * Negative controls receive only contaminant DNA, plus a fraction `leak` of
#   their reads from two random true samples (well-to-well cross-contamination).
# * Reads are drawn as Multinomial(depth, composition).
#
# Simplifications (stated in the report): no taxon is both genuine and a
# contaminant, a single batch, and a single kitome composition.

#' Simulate one study.
#'
#' @return A list with `counts` (samples x taxa), `samples` (sample_id,
#'   is_neg, dna_conc, biomass, contam_dna, total_reads) and `taxonomy`
#'   (asv_id, expected).
simulate_study <- function(n_samples = 40, n_neg = 4,
                           n_true = 150, n_contam = 100,
                           biomass_ratio = 10, biomass_sdlog = 1,
                           contam_dna = 1, contam_sdlog = 0.2,
                           conc_noise = 0.1,
                           depth_meanlog = log(30000), depth_sdlog = 0.3,
                           neg_depth_factor = 0.5,
                           leak = 0.005,
                           subject_heterogeneity = 5,
                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(n_samples >= 2, n_neg >= 0, leak >= 0, leak < 1)

  # Taxon-level parameters -----------------------------------------------------
  p_true <- normalise(stats::rlnorm(n_true, 0, 2))
  presence_prob <- stats::runif(n_true, 0.3, 1)
  p_contam <- normalise(stats::rlnorm(n_contam, 0, 1.5))

  taxa <- c(sprintf("T%03d", seq_len(n_true)), sprintf("C%03d", seq_len(n_contam)))
  n_tax <- length(taxa)

  # True samples ---------------------------------------------------------------
  biomass <- stats::rlnorm(n_samples, log(biomass_ratio * contam_dna), biomass_sdlog)
  c_dna <- contam_dna * stats::rlnorm(n_samples + n_neg, 0, contam_sdlog)

  true_comp <- t(vapply(seq_len(n_samples), function(i) {
    present <- stats::runif(n_true) < presence_prob
    if (!any(present)) present[which.max(p_true)] <- TRUE
    g <- stats::rgamma(n_true, shape = subject_heterogeneity * p_true * n_true / 10 + 1e-3)
    normalise(g * present)
  }, numeric(n_true)))

  contam_presence <- stats::runif(n_contam, 0.3, 1)
  contam_comp <- t(vapply(seq_len(n_samples + n_neg), function(i) {
    present <- stats::runif(n_contam) < contam_presence
    if (!any(present)) present[which.max(p_contam)] <- TRUE
    g <- stats::rgamma(n_contam, shape = 50 * p_contam * n_contam / 10 + 1e-3)
    normalise(g * present)
  }, numeric(n_contam)))

  comp <- matrix(0, n_samples + n_neg, n_tax)
  for (i in seq_len(n_samples)) {
    tot <- biomass[i] + c_dna[i]
    comp[i, ] <- c(biomass[i] * true_comp[i, ], c_dna[i] * contam_comp[i, ]) / tot
  }
  for (j in seq_len(n_neg)) {
    i <- n_samples + j
    # Cross-contamination comes from neighbouring wells: two random samples
    neighbours <- sample(seq_len(n_samples), min(2, n_samples))
    leaked <- normalise(colMeans(comp[neighbours, , drop = FALSE]))
    comp[i, ] <- (1 - leak) * c(rep(0, n_true), contam_comp[i, ]) + leak * leaked
  }

  # Sequencing -----------------------------------------------------------------
  depth <- round(stats::rlnorm(n_samples + n_neg, depth_meanlog, depth_sdlog))
  depth[n_samples + seq_len(n_neg)] <-
    round(depth[n_samples + seq_len(n_neg)] * neg_depth_factor)
  counts <- t(vapply(seq_len(n_samples + n_neg), function(i) {
    as.numeric(stats::rmultinom(1, depth[i], comp[i, ]))
  }, numeric(n_tax)))

  ids <- c(sprintf("S%02d", seq_len(n_samples)), sprintf("NEG%02d", seq_len(n_neg)))
  dimnames(counts) <- list(ids, taxa)

  dna_true <- c(biomass, rep(0, n_neg)) + c_dna
  samples <- data.frame(
    sample_id = ids,
    is_neg = rep(c(FALSE, TRUE), c(n_samples, n_neg)),
    dna_conc = dna_true * stats::rlnorm(n_samples + n_neg, 0, conc_noise),
    biomass = c(biomass, rep(0, n_neg)),
    contam_dna = c_dna,
    total_reads = rowSums(counts)
  )
  taxonomy <- data.frame(asv_id = taxa,
                         expected = rep(c(TRUE, FALSE), c(n_true, n_contam)))

  list(counts = counts, samples = samples, taxonomy = taxonomy)
}

normalise <- function(x) {
  s <- sum(x)
  if (s == 0) return(x)
  x / s
}

#' Score every method on one dataset.
#'
#' decontam is called once with method = "combined", which returns the
#' frequency, prevalence and combined scores together; each method is then
#' thresholded at every value in `thresholds`. This gives exactly the same
#' calls as running isContaminant() separately per method and threshold
#' (the frequency score is computed on true samples only in both cases).
#'
#' @return One row per method x threshold with ASV- and read-level metrics,
#'   computed over taxa observed in at least one true sample.
evaluate_methods <- function(sim, thresholds = c(0.1, 0.5)) {
  counts <- sim$counts
  s <- sim$samples
  expected <- stats::setNames(sim$taxonomy$expected, sim$taxonomy$asv_id)
  is_neg <- s$is_neg

  in_samples <- colSums(counts[!is_neg, , drop = FALSE]) > 0
  samp_counts <- counts[!is_neg, in_samples, drop = FALSE]
  exp_u <- expected[in_samples]

  flags <- list()
  if (any(is_neg)) {
    res <- suppressWarnings(decontam::isContaminant(
      counts, conc = s$dna_conc, neg = is_neg, method = "combined"))
    p_freq <- res$p.freq[in_samples]
    p_prev <- res$p.prev[in_samples]
    p_comb <- res$p[in_samples]
    for (t in thresholds) {
      flags[[paste0("frequency|", t)]] <- lt(p_freq, t)
      flags[[paste0("prevalence|", t)]] <- lt(p_prev, t)
      flags[[paste0("combined|", t)]] <- lt(p_comb, t)
      flags[[paste0("either|", t)]] <- lt(p_freq, t) | lt(p_prev, t)
    }
    flags[["in negative control|NA"]] <- flag_in_negative(counts, is_neg)[in_samples]
  } else {
    res <- suppressWarnings(decontam::isContaminant(
      counts[!is_neg, ], conc = s$dna_conc[!is_neg], method = "frequency"))
    for (t in thresholds) {
      flags[[paste0("frequency|", t)]] <- lt(res$p[in_samples], t)
    }
  }
  flags[["abundance >1%|NA"]] <- flag_low_abundance(counts, is_neg, 0.01, 1)[in_samples]

  do.call(rbind, lapply(names(flags), function(key) {
    f <- stats::setNames(flags[[key]], names(exp_u))
    cm <- classification_metrics(f, exp_u)
    contam_reads <- sum(samp_counts[, !exp_u])
    true_reads <- sum(samp_counts[, exp_u])
    parts <- strsplit(key, "|", fixed = TRUE)[[1]]
    data.frame(
      method = parts[1],
      threshold = suppressWarnings(as.numeric(parts[2])),
      sensitivity = cm$sensitivity,
      specificity = cm$specificity,
      contam_reads_removed = 100 * sum(samp_counts[, f & !exp_u]) / contam_reads,
      true_reads_lost = 100 * sum(samp_counts[, f & exp_u]) / true_reads,
      true_taxa_lost = cm$fp
    )
  }))
}

# TRUE where score < threshold; unscored (NA) taxa are never flagged
lt <- function(p, t) !is.na(p) & p < t
