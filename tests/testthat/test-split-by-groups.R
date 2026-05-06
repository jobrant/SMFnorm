# tests/testthat/test-split-by-groups.R

test_that("split_by_groups creates correct group structure", {
  samples <- make_test_dataset(
    n_samples = 4, n_sites = 50,
    group_ids = c("Treated", "Treated", "Control", "Control"),
    seed = 42
  )

  split_data <- split_by_groups(samples)

  expect_true("Treated" %in% names(split_data))
  expect_true("Control" %in% names(split_data))
  expect_equal(length(split_data$Treated), 2)
  expect_equal(length(split_data$Control), 2)
})


test_that("split_by_groups filters to requested groups", {
  samples <- make_test_dataset(
    n_samples = 6, n_sites = 50,
    group_ids = c("A", "A", "B", "B", "C", "C"),
    seed = 42
  )

  split_data <- split_by_groups(samples, group_names = c("A", "C"))

  expect_equal(length(names(split_data)), 2)
  expect_true("A" %in% names(split_data))
  expect_true("C" %in% names(split_data))
  expect_false("B" %in% names(split_data))
})


test_that("split_by_groups preserves sample metadata", {
  samples <- make_test_dataset(n_samples = 4, n_sites = 50, seed = 42)

  split_data <- split_by_groups(samples)

  expect_false(is.null(attr(split_data, "sample_metadata")))
})


test_that("split_by_groups errors without metadata", {
  samples <- list(
    S1 = make_test_sample(50, seed = 1),
    S2 = make_test_sample(50, seed = 2)
  )
  # No sample_metadata attribute

  expect_error(split_by_groups(samples), "sample metadata")
})


test_that("split_by_groups warns on missing groups", {
  samples <- make_test_dataset(n_samples = 4, n_sites = 50, seed = 42)

  expect_warning(
    split_by_groups(samples, group_names = c("A", "NonExistent")),
    "not found"
  )
})
