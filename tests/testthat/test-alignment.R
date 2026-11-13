# tests/testthat/test-alignment.R
#
# Normalization works row by row across samples, so row i must be the same
# site in every sample. These tests cover the checks and the behaviour that
# depend on that:
#   - misaligned samples are an error, not a silent wrong answer;
#   - the min_coverage filter keeps the same sites in every sample;
#   - each sample comes back under its own name, holding its own data;
#   - the sample_metadata attribute survives normalization.

.norm <- function(d, ...) {
  suppressMessages(normalize_methylation_data(
    d, diagnostics = FALSE, ...
  ))
}

.all_samples <- function(d) unlist(unname(d), recursive = FALSE)


# ── Alignment check ──────────────────────────────────────────────────────────

test_that("rows in a different order in one sample are an error", {
  d <- make_test_split_data(n_reps_per_group = 2, n_sites = 300, seed = 1)
  d$GroupA$GroupA_Rep2 <- d$GroupA$GroupA_Rep2[rev(seq_len(.N))]

  expect_error(.norm(d, min_coverage = 0), "same sites in the same order")
  expect_error(.norm(d, min_coverage = 0), "GroupA_Rep2")
})


test_that("the same number of rows but different sites is an error", {
  d <- make_test_split_data(n_reps_per_group = 2, n_sites = 300, seed = 1)
  d$GroupB$GroupB_Rep1[1, pos := pos + 1L]

  expect_error(.norm(d, min_coverage = 0), "same sites in the same order")
})


test_that("different numbers of sites is an error naming the counts", {
  d <- make_test_split_data(n_reps_per_group = 2, n_sites = 300, seed = 1)
  d$GroupA$GroupA_Rep1 <- d$GroupA$GroupA_Rep1[-1]

  expect_error(.norm(d, min_coverage = 0), "different numbers of sites")
  expect_error(.norm(d, min_coverage = 0), "GroupA_Rep1 = 299")
})


test_that("the rate step on its own also rejects misaligned samples", {
  d <- make_test_split_data(n_reps_per_group = 2, n_sites = 300, seed = 1)
  d$GroupA$GroupA_Rep2 <- d$GroupA$GroupA_Rep2[rev(seq_len(.N))]

  expect_error(
    suppressMessages(normalize_methylation_rates(d, diagnostics = FALSE)),
    "same sites in the same order"
  )
})


test_that("aligned samples pass, with chr as factor in one and character in another", {
  d <- make_test_split_data(n_reps_per_group = 2, n_sites = 300, seed = 1)
  d$GroupA$GroupA_Rep1[, chr := factor(chr)]

  expect_no_error(.norm(d, min_coverage = 0))
})


# ── Coverage filter ──────────────────────────────────────────────────────────

test_that("min_coverage keeps exactly the sites that pass in every sample", {
  d <- make_test_split_data(n_reps_per_group = 3, n_sites = 2000, seed = 7)
  thr <- 12

  # Samples differ in which sites pass, so a per-sample filter would leave
  # them with different rows.
  pass <- lapply(.all_samples(d), function(dt) dt$cov >= thr)
  expect_gt(length(unique(vapply(pass, sum, integer(1)))), 1)
  keep <- Reduce(`&`, pass)
  expect_gt(sum(keep), 100)
  expect_lt(sum(keep), 2000)

  out <- .norm(d, min_coverage = thr, normalize_rates = FALSE)

  expected_pos <- d[[1]][[1]]$pos[keep]
  for (dt in .all_samples(out)) {
    expect_identical(dt$pos, expected_pos)
    expect_true(all(dt$cov >= thr))
  }
})


test_that("normalization runs after a filter that removes different sites per sample", {
  d <- make_test_split_data(n_reps_per_group = 3, n_sites = 2000, seed = 7)

  out <- .norm(d, min_coverage = 12)

  n <- vapply(.all_samples(out), nrow, integer(1))
  expect_length(unique(n), 1)
  ref <- .all_samples(out)[[1]]
  for (dt in .all_samples(out)) {
    expect_identical(dt$pos, ref$pos)
    expect_identical(dt$chr, ref$chr)
  }
})


test_that("a threshold no site reaches in every sample is an error", {
  d <- make_test_split_data(n_reps_per_group = 2, n_sites = 200, seed = 1)

  expect_error(.norm(d, min_coverage = 1e6), "No sites have coverage")
})


test_that("min_coverage = 0 keeps every site", {
  d <- make_test_split_data(n_reps_per_group = 2, n_sites = 200, seed = 1)

  out <- .norm(d, min_coverage = 0)

  expect_true(all(vapply(.all_samples(out), nrow, integer(1)) == 200L))
})


# ── Names, sample identity and metadata ──────────────────────────────────────

test_that("group and sample names are unchanged by normalization", {
  d <- make_test_split_data(n_reps_per_group = 3, n_sites = 500, seed = 3)

  for (between in c(TRUE, FALSE)) {
    out <- quiet_within_only(.norm(d, min_coverage = 0,
                                   rate_between_groups = between))
    expect_identical(names(out), names(d))
    for (g in names(d)) {
      expect_identical(names(out[[g]]), names(d[[g]]))
    }
  }
})


test_that("each sample comes back holding its own data", {
  d <- make_test_split_data(n_reps_per_group = 3, n_sites = 500, seed = 3)

  out <- .norm(d, min_coverage = 0)

  # Rate normalization does not change coverage, and every sample has its own
  # coverage values, so coverage identifies which sample a table came from.
  covs <- lapply(.all_samples(d), `[[`, "cov")
  expect_identical(anyDuplicated(covs), 0L)
  for (g in names(d)) {
    for (s in names(d[[g]])) {
      expect_identical(out[[g]][[s]]$cov, d[[g]][[s]]$cov,
                       label = paste("coverage of", s))
    }
  }
})


test_that("the sample_metadata attribute survives both normalization entry points", {
  d <- make_test_split_data(n_reps_per_group = 2, n_sites = 500, seed = 3)
  meta <- attr(d, "sample_metadata")
  expect_false(is.null(meta))

  full <- suppressMessages(normalize_methylation_data(d, min_coverage = 12,
                                                      diagnostics = FALSE))
  expect_identical(attr(full, "sample_metadata"), meta)

  rates_only <- suppressMessages(normalize_methylation_rates(d, diagnostics = FALSE))
  expect_identical(attr(rates_only, "sample_metadata"), meta)
})
