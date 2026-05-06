# tests/testthat/test-normalization.R

# ── Coverage normalization ───────────────────────────────────────────────────

test_that("normalize_coverage equalizes mean coverage within groups", {
  split_data <- make_test_split_data(
    n_reps_per_group = 3, n_sites = 200, seed = 42
  )

  normalized <- normalize_coverage(split_data, diagnostics = FALSE)

  # Within each group, mean coverage should be approximately equal
  for (g in names(normalized)) {
    mean_covs <- sapply(normalized[[g]], function(dt) mean(dt$cov))
    # All means within 1% of each other
    relative_range <- diff(range(mean_covs)) / mean(mean_covs)
    expect_lt(relative_range, 0.01,
              label = paste("Coverage range within group", g))
  }
})


test_that("normalize_coverage preserves rates", {
  split_data <- make_test_split_data(
    n_reps_per_group = 2, n_sites = 100, seed = 42
  )

  # Record original rates
  orig_rates <- lapply(split_data, function(g) {
    lapply(g, function(dt) dt$rate)
  })

  normalized <- normalize_coverage(split_data, diagnostics = FALSE)

  # Rates should be unchanged (coverage norm scales mc and cov together)
  for (g in names(normalized)) {
    for (r in names(normalized[[g]])) {
      new_rates <- normalized[[g]][[r]]$rate
      old_rates <- orig_rates[[g]][[r]]
      expect_equal(new_rates, old_rates, tolerance = 1e-10,
                   label = paste("Rates preserved for", g, r))
    }
  }
})


test_that("normalize_coverage mc never exceeds cov", {
  split_data <- make_test_split_data(
    n_reps_per_group = 2, n_sites = 200, seed = 42
  )

  normalized <- normalize_coverage(split_data, diagnostics = FALSE)

  for (g in names(normalized)) {
    for (r in names(normalized[[g]])) {
      dt <- normalized[[g]][[r]]
      expect_true(all(dt$mc <= dt$cov),
                  label = paste("mc <= cov for", g, r))
    }
  }
})


# ── Rate normalization ───────────────────────────────────────────────────────

test_that("normalize_methylation_data clamps rates to [0, 1]", {
  split_data <- make_test_split_data(
    n_reps_per_group = 3, n_sites = 500, seed = 42
  )

  normalized <- normalize_methylation_data(
    split_data,
    do_coverage_norm = TRUE,
    normalize_rates = TRUE,
    rate_within_groups = TRUE,
    rate_between_groups = FALSE,
    within_alpha = 0.3,
    min_coverage = 1,
    diagnostics = FALSE
  )

  for (g in names(normalized)) {
    for (r in names(normalized[[g]])) {
      dt <- normalized[[g]][[r]]
      expect_true(all(dt$rate >= 0),
                  label = paste("rate >= 0 for", g, r))
      expect_true(all(dt$rate <= 1),
                  label = paste("rate <= 1 for", g, r))
    }
  }
})


test_that("normalize_methylation_data mc never exceeds cov after normalization", {
  split_data <- make_test_split_data(
    n_reps_per_group = 3, n_sites = 500, seed = 42
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
      expect_true(all(dt$mc <= dt$cov),
                  label = paste("mc <= cov for", g, r))
    }
  }
})


test_that("rate normalization reduces within-group variance", {
  # Create data with intentionally different rate distributions
  set.seed(42)
  split_data <- make_test_split_data(
    n_reps_per_group = 3, n_sites = 500, seed = 42
  )

  # Measure variance before
  var_before <- sapply(split_data, function(g) {
    rates <- do.call(cbind, lapply(g, function(dt) dt$rate))
    mean(apply(rates, 1, var))
  })

  normalized <- normalize_methylation_data(
    split_data,
    do_coverage_norm = FALSE,
    normalize_rates = TRUE,
    rate_within_groups = TRUE,
    rate_between_groups = FALSE,
    within_alpha = 0.3,
    min_coverage = 1,
    diagnostics = FALSE
  )

  # Measure variance after
  var_after <- sapply(normalized, function(g) {
    rates <- do.call(cbind, lapply(g, function(dt) dt$rate))
    mean(apply(rates, 1, var))
  })

  # Variance should decrease
  for (g in names(split_data)) {
    expect_lt(var_after[g], var_before[g],
              label = paste("Variance reduced for group", g))
  }
})


test_that("alpha = 0 collapses all rates to quantile mean", {
  split_data <- make_test_split_data(
    n_reps_per_group = 3, n_sites = 200, seed = 42
  )

  normalized <- normalize_methylation_data(
    split_data,
    do_coverage_norm = FALSE,
    normalize_rates = TRUE,
    rate_within_groups = TRUE,
    rate_between_groups = FALSE,
    within_alpha = 0.0,
    min_coverage = 1,
    diagnostics = FALSE
  )

  # With alpha = 0, all samples within a group should have identical rates
  for (g in names(normalized)) {
    rates <- do.call(cbind, lapply(normalized[[g]], function(dt) dt$rate))
    row_vars <- apply(rates, 1, var)
    expect_true(all(row_vars < 1e-10),
                label = paste("Rates identical at alpha=0 for group", g))
  }
})


test_that("alpha = 1 preserves original rate structure", {
  split_data <- make_test_split_data(
    n_reps_per_group = 3, n_sites = 200, seed = 42
  )

  # Record original rates
  orig_rates <- lapply(split_data, function(g) {
    lapply(g, function(dt) dt$rate)
  })

  normalized <- normalize_methylation_data(
    split_data,
    do_coverage_norm = FALSE,
    normalize_rates = TRUE,
    rate_within_groups = TRUE,
    rate_between_groups = FALSE,
    within_alpha = 1.0,
    min_coverage = 1,
    diagnostics = FALSE
  )

  # With alpha = 1, sample-specific deviation is fully preserved
  # Rates should be shifted (mean adjusted) but the relative differences
  # between samples at each site should be preserved
  for (g in names(normalized)) {
    for (r in names(normalized[[g]])) {
      old <- orig_rates[[g]][[r]]
      new <- normalized[[g]][[r]]$rate
      # Correlation should be very high (structure preserved)
      cor_val <- cor(old, new)
      expect_gt(cor_val, 0.95,
                label = paste("Rate structure preserved at alpha=1 for", g, r))
    }
  }
})


test_that("min_coverage filter removes low-coverage sites", {
  split_data <- make_test_split_data(
    n_reps_per_group = 2, n_sites = 200, seed = 42
  )

  # Artificially set some sites to low coverage
  split_data[[1]][[1]][1:10, cov := 2L]
  split_data[[1]][[1]][1:10, mc := 1L]
  split_data[[1]][[1]][1:10, rate := mc / cov]

  normalized <- normalize_methylation_data(
    split_data,
    do_coverage_norm = FALSE,
    normalize_rates = FALSE,
    min_coverage = 5,
    diagnostics = FALSE
  )

  # Should have fewer sites after filtering
  for (g in names(normalized)) {
    for (r in names(normalized[[g]])) {
      expect_true(all(normalized[[g]][[r]]$cov >= 5),
                  label = paste("Min coverage enforced for", g, r))
    }
  }
})


# ── Between-group normalization ──────────────────────────────────────────────

test_that("between-group normalization reduces group-level differences", {
  split_data <- make_test_split_data(
    n_reps_per_group = 3, n_sites = 500,
    group_names = c("High", "Low"), seed = 42
  )

  # Make one group systematically higher
  for (r in names(split_data$High)) {
    split_data$High[[r]][, rate := pmin(1, rate + 0.15)]
    split_data$High[[r]][, mc := as.integer(round(cov * rate))]
  }

  # Measure group means before
  mean_before <- sapply(split_data, function(g) {
    mean(sapply(g, function(dt) mean(dt$rate)))
  })
  diff_before <- abs(diff(mean_before))

  normalized <- normalize_methylation_data(
    split_data,
    do_coverage_norm = FALSE,
    normalize_rates = TRUE,
    rate_within_groups = TRUE,
    rate_between_groups = TRUE,
    within_alpha = 0.3,
    between_alpha = 0.5,
    min_coverage = 1,
    diagnostics = FALSE
  )

  # Measure group means after
  mean_after <- sapply(normalized, function(g) {
    mean(sapply(g, function(dt) mean(dt$rate)))
  })
  diff_after <- abs(diff(mean_after))

  # Between-group difference should be reduced
  expect_lt(diff_after, diff_before,
            label = "Between-group difference reduced")
})
