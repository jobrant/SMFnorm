# tests/testthat/test-data-integrity.R

# ── Rate and count consistency ───────────────────────────────────────────────

test_that("rate always equals mc/cov after normalization", {
  split_data <- make_test_split_data(
    n_reps_per_group = 3, n_sites = 300, seed = 42
  )

  normalized <- normalize_methylation_data(
    split_data,
    do_coverage_norm = TRUE,
    normalize_rates = TRUE,
    rate_within_groups = TRUE,
    rate_between_groups = TRUE,
    within_alpha = 0.3,
    between_alpha = 0.5,
    min_coverage = 1,
    diagnostics = FALSE
  )

  for (g in names(normalized)) {
    for (r in names(normalized[[g]])) {
      dt <- normalized[[g]][[r]]
      # mc should be within 1 of rate * cov (integer rounding)
      expect_true(all(abs(dt$mc - dt$rate * dt$cov) <= 1),
                  label = paste("mc ≈ rate * cov for", g, r))
      # mc must be non-negative integer not exceeding cov
      expect_true(all(dt$mc >= 0))
      expect_true(all(dt$mc <= dt$cov))
    }
  }
})


test_that("mc values are non-negative after normalization", {
  split_data <- make_test_split_data(
    n_reps_per_group = 3, n_sites = 300, seed = 42
  )

  normalized <- normalize_methylation_data(
    split_data,
    do_coverage_norm = TRUE,
    normalize_rates = TRUE,
    within_alpha = 0.3,
    min_coverage = 1,
    diagnostics = FALSE
  )

  for (g in names(normalized)) {
    for (r in names(normalized[[g]])) {
      dt <- normalized[[g]][[r]]
      expect_true(all(dt$mc >= 0),
                  label = paste("mc >= 0 for", g, r))
    }
  }
})


test_that("cov values are positive after normalization", {
  split_data <- make_test_split_data(
    n_reps_per_group = 2, n_sites = 100, seed = 42
  )

  normalized <- normalize_methylation_data(
    split_data,
    do_coverage_norm = TRUE,
    normalize_rates = FALSE,
    min_coverage = 1,
    diagnostics = FALSE
  )

  for (g in names(normalized)) {
    for (r in names(normalized[[g]])) {
      dt <- normalized[[g]][[r]]
      expect_true(all(dt$cov > 0),
                  label = paste("cov > 0 for", g, r))
    }
  }
})


# ── No NaN/NA/Inf values ────────────────────────────────────────────────────

test_that("no NaN, NA, or Inf values after full normalization", {
  split_data <- make_test_split_data(
    n_reps_per_group = 3, n_sites = 300, seed = 42
  )

  normalized <- normalize_methylation_data(
    split_data,
    do_coverage_norm = TRUE,
    normalize_rates = TRUE,
    rate_within_groups = TRUE,
    rate_between_groups = TRUE,
    within_alpha = 0.3,
    between_alpha = 0.5,
    min_coverage = 1,
    diagnostics = FALSE
  )

  for (g in names(normalized)) {
    for (r in names(normalized[[g]])) {
      dt <- normalized[[g]][[r]]
      expect_false(any(is.na(dt$rate)),
                   label = paste("No NA rates for", g, r))
      expect_false(any(is.nan(dt$rate)),
                   label = paste("No NaN rates for", g, r))
      expect_false(any(is.infinite(dt$rate)),
                   label = paste("No Inf rates for", g, r))
      expect_false(any(is.na(dt$mc)),
                   label = paste("No NA mc for", g, r))
      expect_false(any(is.na(dt$cov)),
                   label = paste("No NA cov for", g, r))
    }
  }
})


# ── Edge cases ───────────────────────────────────────────────────────────────

test_that("normalization handles sites with rate = 0", {
  split_data <- make_test_split_data(
    n_reps_per_group = 2, n_sites = 100, seed = 42
  )

  # Set first 10 sites to rate = 0 in all samples
  for (g in names(split_data)) {
    for (r in names(split_data[[g]])) {
      split_data[[g]][[r]][1:10, `:=`(mc = 0L, rate = 0.0)]
    }
  }

  # Should not error
  normalized <- normalize_methylation_data(
    split_data,
    do_coverage_norm = TRUE,
    normalize_rates = TRUE,
    within_alpha = 0.3,
    min_coverage = 1,
    diagnostics = FALSE
  )

  # Rates should still be non-negative
  for (g in names(normalized)) {
    for (r in names(normalized[[g]])) {
      expect_true(all(normalized[[g]][[r]]$rate >= 0))
    }
  }
})


test_that("normalization handles sites with rate = 1", {
  split_data <- make_test_split_data(
    n_reps_per_group = 2, n_sites = 100, seed = 42
  )

  # Set first 10 sites to rate = 1 in all samples
  for (g in names(split_data)) {
    for (r in names(split_data[[g]])) {
      covs <- split_data[[g]][[r]]$cov[1:10]
      split_data[[g]][[r]][1:10, `:=`(mc = covs, rate = 1.0)]
    }
  }

  # Should not error
  normalized <- normalize_methylation_data(
    split_data,
    do_coverage_norm = TRUE,
    normalize_rates = TRUE,
    within_alpha = 0.3,
    min_coverage = 1,
    diagnostics = FALSE
  )

  # Rates should not exceed 1
  for (g in names(normalized)) {
    for (r in names(normalized[[g]])) {
      expect_true(all(normalized[[g]][[r]]$rate <= 1))
    }
  }
})


test_that("normalization handles two samples per group", {
  split_data <- make_test_split_data(
    n_reps_per_group = 2, n_sites = 100, seed = 42
  )

  # Minimum viable design — should not error
  expect_no_error({
    normalized <- normalize_methylation_data(
      split_data,
      do_coverage_norm = TRUE,
      normalize_rates = TRUE,
      within_alpha = 0.3,
      min_coverage = 1,
      diagnostics = FALSE
    )
  })
})


# ── Structure preservation ───────────────────────────────────────────────────

test_that("normalization preserves number of samples and groups", {
  split_data <- make_test_split_data(
    n_reps_per_group = 3, n_sites = 100,
    group_names = c("X", "Y", "Z"), seed = 42
  )

  normalized <- normalize_methylation_data(
    split_data,
    do_coverage_norm = TRUE,
    normalize_rates = TRUE,
    within_alpha = 0.3,
    min_coverage = 1,
    diagnostics = FALSE
  )

  expect_equal(length(normalized), 3)
  expect_equal(names(normalized), c("X", "Y", "Z"))
  for (g in names(normalized)) {
    expect_equal(length(normalized[[g]]), 3)
  }
})


test_that("normalization preserves genomic coordinate columns", {
  split_data <- make_test_split_data(
    n_reps_per_group = 2, n_sites = 50, seed = 42
  )

  # Record original coordinates
  orig_chr <- split_data[[1]][[1]]$chr
  orig_pos <- split_data[[1]][[1]]$pos
  orig_strand <- split_data[[1]][[1]]$strand

  normalized <- normalize_methylation_data(
    split_data,
    do_coverage_norm = TRUE,
    normalize_rates = TRUE,
    within_alpha = 0.3,
    min_coverage = 1,
    diagnostics = FALSE
  )

  # chr, pos, strand should be unchanged
  expect_equal(normalized[[1]][[1]]$chr, orig_chr)
  expect_equal(normalized[[1]][[1]]$pos, orig_pos)
  expect_equal(normalized[[1]][[1]]$strand, orig_strand)
})
