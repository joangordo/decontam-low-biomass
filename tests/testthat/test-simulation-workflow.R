root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)

test_that("simulate_study is reproducible and well formed", {
  a <- simulate_study(n_samples = 10, n_neg = 3, n_true = 30, n_contam = 20, seed = 42)
  b <- simulate_study(n_samples = 10, n_neg = 3, n_true = 30, n_contam = 20, seed = 42)
  expect_identical(a, b)
  expect_equal(dim(a$counts), c(13, 50))
  expect_identical(rownames(a$counts), a$samples$sample_id)
  expect_identical(colnames(a$counts), a$taxonomy$asv_id)
  expect_equal(sum(a$samples$is_neg), 3)
  expect_equal(sum(a$taxonomy$expected), 30)
  expect_true(all(a$samples$dna_conc > 0))
})

test_that("without leakage, negative controls contain no genuine taxa", {
  s <- simulate_study(n_samples = 10, n_neg = 4, n_true = 30, n_contam = 20,
                      leak = 0, seed = 1)
  expect_equal(sum(s$counts[s$samples$is_neg, s$taxonomy$expected]), 0)
  s2 <- simulate_study(n_samples = 10, n_neg = 4, n_true = 30, n_contam = 20,
                       leak = 0.05, seed = 1)
  expect_gt(sum(s2$counts[s2$samples$is_neg, s2$taxonomy$expected]), 0)
})

test_that("lower biomass means a larger contaminant share", {
  share <- function(r) {
    s <- simulate_study(n_samples = 30, n_neg = 2, biomass_ratio = r, seed = 7)
    tc <- s$counts[!s$samples$is_neg, ]
    stats::median(rowSums(tc[, !s$taxonomy$expected]) / rowSums(tc))
  }
  expect_gt(share(1), share(100))
})

test_that("evaluate_methods returns one row per method x threshold", {
  s <- simulate_study(n_samples = 15, n_neg = 3, n_true = 30, n_contam = 20, seed = 3)
  ev <- evaluate_methods(s, thresholds = c(0.1, 0.5))
  expect_equal(nrow(ev), 4 * 2 + 2)
  expect_true(all(ev$sensitivity >= 0 & ev$sensitivity <= 1, na.rm = TRUE))
  expect_true(all(ev$true_reads_lost >= 0 & ev$true_reads_lost <= 100))
})

test_that("the combined call gives the same flags as separate calls", {
  s <- simulate_study(n_samples = 15, n_neg = 3, n_true = 30, n_contam = 20, seed = 5)
  sep <- run_decontam(s$counts, s$samples$dna_conc, s$samples$is_neg,
                      "frequency", 0.3)
  comb <- suppressWarnings(decontam::isContaminant(
    s$counts, conc = s$samples$dna_conc, neg = s$samples$is_neg,
    method = "combined"))
  expect_equal(sep$contaminant, lt(comb$p.freq, 0.3))
})

test_that("clean_low_biomass handles groups, controls and fallbacks", {
  a <- simulate_study(n_samples = 12, n_neg = 3, n_true = 30, n_contam = 20, seed = 11)
  b <- simulate_study(n_samples = 12, n_neg = 0, n_true = 30, n_contam = 20, seed = 12)
  rownames(b$counts) <- b$samples$sample_id <- paste0("B_", b$samples$sample_id)
  counts <- rbind(a$counts, b$counts)
  samples <- rbind(cbind(a$samples, sample_type = "swab"),
                   cbind(b$samples, sample_type = "blood"))
  counts["S01", ] <- 0; counts["S01", 1] <- 10          # shallow sample
  samples$total_reads <- rowSums(counts)

  res <- clean_low_biomass(counts, samples, method = "combined",
                           threshold = 0.1, low_read_rule = c(swab = "fixed"))
  expect_equal(res$summary$method, c("combined", "frequency"))
  expect_true("S01" %in% res$dropped_samples$sample_id)
  expect_false(any(grepl("^NEG", rownames(res$counts))))
  expect_equal(nrow(res$counts), 12 - 1 + 12)
  expect_true(all(res$counts >= 0))
})

test_that("every method gives decontam's own calls", {
  # Regression test: for "either" and "both" decontam returns p = NA and
  # decides internally, which an earlier version of run_decontam() missed.
  s <- simulate_study(n_samples = 15, n_neg = 4, n_true = 30, n_contam = 20,
                      leak = 0.005, seed = 21)
  for (m in c("frequency", "prevalence", "combined", "minimum",
              "either", "both")) {
    for (t in c(0.1, 0.5)) {
      ours <- run_decontam(s$counts, s$samples$dna_conc, s$samples$is_neg,
                           m, t)
      ref <- suppressMessages(suppressWarnings(decontam::isContaminant(
        s$counts, conc = s$samples$dna_conc, neg = s$samples$is_neg,
        method = m, threshold = t)))
      expect_identical(ours$contaminant, unname(ref$contaminant),
                       label = paste(m, t))
    }
  }
  either <- run_decontam(s$counts, s$samples$dna_conc, s$samples$is_neg,
                         "either", 0.5)
  both <- run_decontam(s$counts, s$samples$dna_conc, s$samples$is_neg,
                       "both", 0.5)
  expect_gt(sum(either$contaminant), 0)
  expect_true(all(either$contaminant[both$contaminant]))   # both implies either
  expect_true(all(either$scored >= both$scored))
})

test_that("clean_low_biomass removes taxa with the either method", {
  s <- simulate_study(n_samples = 15, n_neg = 4, n_true = 30, n_contam = 20,
                      leak = 0.005, seed = 21)
  samples <- cbind(s$samples, sample_type = "swab")
  res <- clean_low_biomass(s$counts, samples, method = "either",
                           threshold = 0.5, min_reads = 0)
  expect_gt(res$summary$taxa_removed, 0)
})

test_that("non-positive DNA concentrations are treated as missing", {
  s <- simulate_study(n_samples = 15, n_neg = 3, n_true = 30, n_contam = 20,
                      seed = 8)
  conc <- s$samples$dna_conc
  conc[s$samples$is_neg] <- 0                       # blanks below detection
  conc[which(!s$samples$is_neg)[1]] <- 0            # one sample too
  expect_message(
    res <- run_decontam(s$counts, conc, s$samples$is_neg, "combined", 0.1),
    "zero or below")
  keep <- !s$samples$is_neg & conc > 0
  ref <- run_decontam(s$counts[keep, ], conc[keep], NULL, "frequency", 0.1)
  expect_equal(res$score_freq, ref$score)           # same frequency evidence
  prev <- run_decontam(s$counts, NULL, s$samples$is_neg, "prevalence", 0.1)
  expect_equal(res$score_prev, prev$score)          # blanks still used
  expect_error(run_decontam(s$counts, rep(0, nrow(s$counts)), s$samples$is_neg,
                            "frequency"), "at least two")
})

test_that("methods_text describes exactly what clean_low_biomass did", {
  s <- simulate_study(n_samples = 15, n_neg = 3, n_true = 30, n_contam = 20, seed = 5)
  samples <- cbind(s$samples, sample_type = "swab")
  res <- clean_low_biomass(s$counts, samples, method = "frequency", threshold = 0.1,
                           min_reads = 0)
  txt <- methods_text(res)
  sm <- res$summary
  expect_type(txt, "character")
  expect_length(txt, 1)
  expect_match(txt, "decontam frequency method (threshold 0.10)", fixed = TRUE)
  expect_match(txt, sprintf("%d of %d scored taxa", sm$taxa_removed, sm$taxa_scored),
               fixed = TRUE)
  expect_match(txt, "3 negative controls", fixed = TRUE)
  expect_match(txt, "DNA concentration", fixed = TRUE)
  expect_false(grepl("presence in negative controls", txt))
})
