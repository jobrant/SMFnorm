# tests/testthat/test-find-shared-sites.R

test_that("find_shared_sites returns correct number of shared sites", {
  samples <- make_test_dataset(n_samples = 3, n_sites = 100,
                                shared_frac = 0.8, 
                               group_ids = c("A", "A", "B"), seed = 42)

  shared_ids <- find_shared_sites(samples, filter = FALSE, quiet = TRUE)

  # Should find exactly the 80 shared sites
  expect_equal(length(shared_ids), 80)
})


test_that("find_shared_sites filters all samples to same size", {
  samples <- make_test_dataset(n_samples = 3, n_sites = 100,
                                shared_frac = 0.8, 
                               group_ids = c("A", "A", "B"),seed = 42)

  filtered <- find_shared_sites(samples, filter = TRUE, quiet = TRUE)

  # All samples should have the same number of rows
  row_counts <- sapply(filtered, nrow)
  expect_true(all(row_counts == row_counts[1]))
  expect_equal(unname(row_counts[1]), 80)
})


test_that("find_shared_sites preserves attributes", {
  samples <- make_test_dataset(n_samples = 3, n_sites = 50, 
                               group_ids = c("A", "A", "B"),seed = 42)
  original_meta <- attr(samples, "sample_metadata")

  filtered <- find_shared_sites(samples, filter = TRUE, quiet = TRUE)

  expect_equal(attr(filtered, "sample_metadata"), original_meta)
})


test_that("find_shared_sites with 100% overlap returns all sites", {
  samples <- make_test_dataset(n_samples = 2, n_sites = 50,
                                shared_frac = 1.0, seed = 42)

  shared_ids <- find_shared_sites(samples, filter = FALSE, quiet = TRUE)

  expect_equal(length(shared_ids), 50)
})


test_that("find_shared_sites errors on bad input", {
  expect_error(find_shared_sites(list()), "non-empty list")
  expect_error(find_shared_sites("not a list"), "non-empty list")

  # Missing uniqueID column
  bad_sample <- list(data.table(chr = "1", pos = 1, mc = 1, cov = 10))
  expect_error(find_shared_sites(bad_sample), "uid")
})
