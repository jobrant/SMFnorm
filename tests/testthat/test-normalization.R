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

  normalized <- quiet_within_only(normalize_methylation_data(
    split_data,
    do_coverage_norm = TRUE,
    normalize_rates = TRUE,
    rate_within_groups = TRUE,
    rate_between_groups = FALSE,
    within_alpha = 0.3,
    min_coverage = 1,
    diagnostics = FALSE
  ))

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

  normalized <- quiet_within_only(normalize_methylation_data(
    split_data,
    do_coverage_norm = FALSE,
    normalize_rates = TRUE,
    rate_within_groups = TRUE,
    rate_between_groups = FALSE,
    within_alpha = 0.3,
    min_coverage = 1,
    diagnostics = FALSE
  ))

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

  normalized <- quiet_within_only(normalize_methylation_data(
    split_data,
    do_coverage_norm = FALSE,
    normalize_rates = TRUE,
    rate_within_groups = TRUE,
    rate_between_groups = FALSE,
    within_alpha = 0.0,
    min_coverage = 1,
    diagnostics = FALSE
  ))

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

  normalized <- quiet_within_only(normalize_methylation_data(
    split_data,
    do_coverage_norm = FALSE,
    normalize_rates = TRUE,
    rate_within_groups = TRUE,
    rate_between_groups = FALSE,
    within_alpha = 1.0,
    min_coverage = 1,
    diagnostics = FALSE
  ))

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


# ── What alpha does ──────────────────────────────────────────────────────────
#
# Within each rate bin, .normalize_methylation_within_set() sets
#   new = mbar + alpha * (x - m_s)
# where m_s is the sample's own bin mean and mbar the across-sample bin mean.
# Replacing m_s by mbar IS the normalization and is applied in full at any
# alpha; alpha only scales site-level deviations. alpha = 1 is therefore not
# the identity. These tests pin that down: an earlier reading ("alpha = 1
# leaves rates unchanged", "higher between_alpha preserves more of the group
# difference") led to defaults that shrank noise and inflated false positives.

.alpha_test_set <- function(effs = c(1.0, 0.85, 0.7), n = 4000, seed = 7) {
  set.seed(seed)
  p <- stats::runif(n, 0.25, 0.75)          # away from 0/1 so clipping is moot
  setNames(lapply(seq_along(effs), function(i) {
    cov <- rep(200L, n)
    rate <- p * effs[i] + stats::rnorm(n, 0, 0.01)
    data.table(cov = cov, rate = rate, mc = as.integer(round(cov * rate)))
  }), paste0("S", seq_along(effs)))
}

.norm_rates <- function(set, alpha) {
  out <- .normalize_methylation_within_set(set, alpha = alpha,
                                           sites_per_quantile = 200,
                                           diagnostics = FALSE)
  do.call(cbind, lapply(out, `[[`, "rate"))
}


test_that("defaults apply no shrinkage (alpha = 1)", {
  expect_identical(formals(normalize_methylation_data)$within_alpha, 1)
  expect_identical(formals(normalize_methylation_data)$between_alpha, 1)
  expect_identical(formals(normalize_methylation_rates)$within_alpha, 1)
  expect_identical(formals(normalize_methylation_rates)$between_alpha, 1)
})


test_that("alpha = 1 is not the identity: it removes per-sample efficiency", {
  set <- .alpha_test_set()
  before <- vapply(set, function(d) mean(d$rate), numeric(1))
  after  <- colMeans(.norm_rates(set, alpha = 1))

  # Samples start ~30% apart in mean rate and end essentially equal.
  expect_gt(diff(range(before)) / mean(before), 0.25)
  expect_lt(diff(range(after)) / mean(after), 0.01)
})


test_that("the efficiency correction does not depend on alpha", {
  set <- .alpha_test_set()
  for (a in c(1, 0.5, 0.3)) {
    after <- colMeans(.norm_rates(set, alpha = a))
    expect_lt(diff(range(after)) / mean(after), 0.01,
              label = paste("sample means aligned at alpha =", a))
  }
})


test_that("alpha only scales site-level deviations from the bin mean", {
  set <- .alpha_test_set()
  r1  <- .norm_rates(set, alpha = 1)
  r05 <- .norm_rates(set, alpha = 0.5)

  # new_a = mbar + a * (x - m_s), so 2 * new_0.5 - new_1 = mbar: identical
  # across samples at every site.
  mbar <- 2 * r05 - r1
  expect_lt(max(apply(mbar, 1, function(v) diff(range(v)))), 1e-9)

  # And deviations around mbar at alpha = 0.5 are exactly half those at 1.
  expect_equal(r05 - mbar, 0.5 * (r1 - mbar), tolerance = 1e-9)
})


test_that("alpha < 1 shrinks replicate variation; alpha = 1 does not", {
  set <- .alpha_test_set(effs = c(1, 1, 1))     # no efficiency difference
  raw <- do.call(cbind, lapply(set, `[[`, "rate"))
  sd_of <- function(m) mean(apply(m, 1, stats::sd))

  expect_equal(sd_of(.norm_rates(set, alpha = 1)), sd_of(raw), tolerance = 0.02)
  expect_equal(sd_of(.norm_rates(set, alpha = 0.5)), 0.5 * sd_of(raw),
               tolerance = 0.02)
})


test_that("between-group step at alpha = 1 removes a global group shift", {
  split_data <- make_test_split_data(
    n_reps_per_group = 3, n_sites = 500,
    group_names = c("High", "Low"), seed = 42
  )
  for (r in names(split_data$Low)) {           # group-wide efficiency loss
    split_data$Low[[r]][, rate := rate * 0.7]
    split_data$Low[[r]][, mc := as.integer(round(cov * rate))]
  }
  gm <- function(d) sapply(d, function(g) mean(sapply(g, function(dt) mean(dt$rate))))
  before <- abs(diff(gm(split_data)))

  normalized <- normalize_methylation_data(
    split_data, do_coverage_norm = FALSE, normalize_rates = TRUE,
    rate_within_groups = TRUE, rate_between_groups = TRUE,
    within_alpha = 1, between_alpha = 1, min_coverage = 1,
    diagnostics = FALSE
  )
  after <- abs(diff(gm(normalized)))

  expect_lt(after, 0.1 * before)
})


# ── Within-group-only normalization with unequal group levels ────────────────
#
# Within-group normalization keeps each group's level but removes replicate
# spread, so an efficiency difference between groups looks MORE significant
# than in the raw data (M1 real null: 22 -> 907 false positives). A warning
# flags this when the between-group step is off.

.groups_with_gap <- function(gap) {
  split_data <- make_test_split_data(
    n_reps_per_group = 3, n_sites = 500,
    group_names = c("A", "B"), seed = 42
  )
  # make_test_split_data() gives the groups different mean rates, so first
  # rescale B to A's level, then impose the requested relative gap.
  gm <- function(g) mean(sapply(g, function(dt) mean(dt$rate)))
  scale <- gm(split_data$A) / gm(split_data$B) * (1 - gap)
  for (r in names(split_data$B)) {
    split_data$B[[r]][, rate := pmin(1, rate * scale)]
    split_data$B[[r]][, mc := as.integer(round(cov * rate))]
  }
  split_data
}

.norm_quiet <- function(d, rbg) {
  normalize_methylation_data(
    d, do_coverage_norm = FALSE, normalize_rates = TRUE,
    rate_within_groups = TRUE, rate_between_groups = rbg,
    min_coverage = 1, diagnostics = FALSE
  )
}

test_that("within-only normalization warns when group levels differ", {
  expect_warning(.norm_quiet(.groups_with_gap(0.3), rbg = FALSE),
                 "only within-group rate normalization")
})

test_that("no warning when the between-group step is on", {
  expect_no_warning(.norm_quiet(.groups_with_gap(0.3), rbg = TRUE))
})

test_that("no warning when group levels match", {
  expect_no_warning(.norm_quiet(.groups_with_gap(0), rbg = FALSE))
})

test_that("within-only normalization keeps the group gap and shrinks spread", {
  d <- .groups_with_gap(0.3)
  gm <- function(x) sapply(x, function(g) mean(sapply(g, function(dt) mean(dt$rate))))
  sp <- function(x) mean(sapply(x, function(g) {
    mean(apply(do.call(cbind, lapply(g, `[[`, "rate")), 1, stats::sd))
  }))
  after <- suppressWarnings(.norm_quiet(d, rbg = FALSE))

  expect_equal(abs(diff(gm(after))), abs(diff(gm(d))), tolerance = 0.05)
  expect_lte(sp(after), sp(d))
})
