# Small hand-checkable examples for every helper.
# Run from the project root with: testthat::test_dir("tests/testthat")

root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)

toy <- matrix(
  c(90, 10,  0,  0,
    80, 15,  5,  0,
     0,  0, 50, 50),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("s1", "s2", "neg"), c("A", "B", "C", "D"))
)
is_neg <- c(FALSE, FALSE, TRUE)
expected <- c(A = TRUE, B = TRUE, C = FALSE, D = FALSE)

test_that("validate_counts rejects bad input", {
  bad <- toy; bad[1, 1] <- -1
  expect_error(validate_counts(bad), "negative")
  bad <- toy; bad[1, 1] <- NA
  expect_error(validate_counts(bad), "NA")
  expect_error(validate_counts(unname(toy)), "row names")
})

test_that("relative_abundance handles empty samples", {
  z <- rbind(toy, empty = 0)
  rel <- relative_abundance(z)
  expect_equal(unname(rowSums(rel[1:3, ])), c(1, 1, 1))
  expect_equal(unname(rel["empty", ]), rep(0, 4))
})

test_that("abundance pre-filter ignores negative controls", {
  # C and D are abundant only in the negative control -> flagged
  # (a filter computed over all samples would have kept them)
  f <- flag_low_abundance(toy, is_neg, min_rel = 0.10)
  expect_equal(f, c(A = FALSE, B = FALSE, C = TRUE, D = TRUE))
})

test_that("negative-control filters behave as documented", {
  expect_equal(flag_in_negative(toy, is_neg),
               c(A = FALSE, B = FALSE, C = TRUE, D = TRUE))
  expect_equal(flag_control_dominant(toy, is_neg),
               c(A = FALSE, B = FALSE, C = TRUE, D = TRUE))
})

test_that("low-read rules never drop controls", {
  # every sample has 100 reads
  f <- flag_low_read_samples(toy, is_neg, rule = "fixed", min_reads = 101)
  expect_equal(names(f), c("s1", "s2", "neg"))
  expect_equal(as.vector(f), c(TRUE, TRUE, FALSE))
  g <- flag_low_read_samples(toy, is_neg, rule = "neg_median")
  expect_equal(attr(g, "cutoff"), 100)
  expect_equal(as.vector(g), c(FALSE, FALSE, FALSE))
})

test_that("classification metrics count correctly", {
  m <- classification_metrics(c(A = FALSE, B = TRUE, C = TRUE, D = FALSE),
                              expected)
  expect_equal(c(m$tp, m$fp, m$fn, m$tn), c(1, 1, 1, 1))
  expect_equal(m$sensitivity, 0.5)
  expect_equal(m$specificity, 0.5)
})

test_that("read-level metrics are consistent", {
  r <- read_level_metrics(toy[1:2, ], c(A = FALSE, B = FALSE, C = TRUE, D = TRUE),
                          expected)
  expect_equal(r$pct_contam_before, c(0, 5))
  expect_equal(r$pct_contam_after, c(0, 0))
  expect_equal(r$pct_mock_reads_lost, c(0, 0))
})

test_that("score_auc is 1 for perfect separation and 0.5 for ties", {
  expect_equal(score_auc(c(0.9, 0.8, 0.1, 0.05), expected), 1)
  expect_equal(score_auc(rep(0.5, 4), expected), 0.5)
  expect_equal(score_auc(c(0.9, NA, 0.1, NA), expected), 1)
})

test_that("load_dataset re-orders and validates the real data", {
  d <- load_dataset(file.path(root, "data"))
  expect_identical(rownames(d$counts), d$samples$sample_id)
  expect_identical(colnames(d$counts), d$taxonomy$asv_id)
  expect_equal(sum(d$taxonomy$expected), 9)
})
