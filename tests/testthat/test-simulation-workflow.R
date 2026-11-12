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
