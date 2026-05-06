# tests/testthat/test-create-bigwig.R

# Skip all tests if rtracklayer is not available
skip_if_not_installed("rtracklayer")

# ── Helper: create test data with chr prefix (BigWig requires it) ────────────

make_bigwig_test_data <- function(seed = 42) {
  set.seed(seed)

  samples <- make_test_dataset(
    n_samples = 4,
    n_sites = 50,
    shared_frac = 1.0,
    group_ids = c("GroupA", "GroupA", "GroupB", "GroupB"),
    seed = seed
  )

  # BigWig requires chr prefix
  samples <- lapply(samples, function(dt) {
    dt[, chr := paste0("chr", chr)]
    dt[, uniqueID := paste(chr, pos, site, sep = "_")]
    dt
  })
  attr(samples, "sample_metadata") <- data.table(
    group_id = c("GroupA", "GroupA", "GroupB", "GroupB"),
    replicate = 1:4,
    sample_id = names(samples),
    file_name = paste0(names(samples), ".tsv.gz")
  )

  return(samples)
}


# ── Single sample BigWig tests ───────────────────────────────────────────────

test_that("create_bigwig produces a file for single sample meth type", {
  test_data <- make_bigwig_test_data()
  out_dir <- tempdir()

  result_path <- create_bigwig(
    meth_data = test_data,
    out = out_dir,
    genome = "hg38",
    sample_name = names(test_data)[1],
    type = "meth"
  )

  expect_true(file.exists(result_path))
  file.remove(result_path)
})


test_that("create_bigwig produces a file for single sample cov type", {
  test_data <- make_bigwig_test_data()
  out_dir <- tempdir()

  result_path <- create_bigwig(
    meth_data = test_data,
    out = out_dir,
    genome = "hg38",
    sample_name = names(test_data)[1],
    type = "cov"
  )

  expect_true(file.exists(result_path))
  file.remove(result_path)
})


test_that("single sample meth BigWig values match original rates", {
  test_data <- make_bigwig_test_data()
  out_dir <- tempdir()
  sample_name <- names(test_data)[1]

  result_path <- create_bigwig(
    meth_data = test_data,
    out = out_dir,
    genome = "hg38",
    sample_name = sample_name,
    type = "meth"
  )

  bw <- rtracklayer::import.bw(result_path)

  # Pick a position from the test data and verify
  orig <- test_data[[sample_name]]
  test_row <- orig[1]

  bw_match <- bw[GenomicRanges::seqnames(bw) == test_row$chr &
                    GenomicRanges::start(bw) == test_row$pos]

  if (length(bw_match) > 0) {
    expect_equal(bw_match$score, test_row$rate, tolerance = 1e-6)
  }

  file.remove(result_path)
})


test_that("single sample cov BigWig values match original coverage", {
  test_data <- make_bigwig_test_data()
  out_dir <- tempdir()
  sample_name <- names(test_data)[1]

  result_path <- create_bigwig(
    meth_data = test_data,
    out = out_dir,
    genome = "hg38",
    sample_name = sample_name,
    type = "cov"
  )

  bw <- rtracklayer::import.bw(result_path)

  orig <- test_data[[sample_name]]
  test_row <- orig[1]

  bw_match <- bw[GenomicRanges::seqnames(bw) == test_row$chr &
                    GenomicRanges::start(bw) == test_row$pos]

  if (length(bw_match) > 0) {
    expect_equal(bw_match$score, test_row$cov, tolerance = 1e-6)
  }

  file.remove(result_path)
})


# ── Aggregated BigWig tests ──────────────────────────────────────────────────

test_that("create_bigwig produces a file for aggregated replicates", {
  test_data <- make_bigwig_test_data()
  out_dir <- tempdir()

  result_path <- create_bigwig(
    meth_data = test_data,
    out = out_dir,
    genome = "hg38",
    aggregate_replicates = TRUE,
    group_name = "GroupA",
    type = "meth"
  )

  expect_true(file.exists(result_path))
  file.remove(result_path)
})


test_that("aggregated meth BigWig values are correct weighted means", {
  test_data <- make_bigwig_test_data()
  out_dir <- tempdir()
  metadata <- attr(test_data, "sample_metadata")
  group_samples <- metadata[group_id == "GroupA", sample_id]

  result_path <- create_bigwig(
    meth_data = test_data,
    out = out_dir,
    genome = "hg38",
    aggregate_replicates = TRUE,
    group_name = "GroupA",
    type = "meth"
  )

  bw <- rtracklayer::import.bw(result_path)

  # Compute expected aggregated rate at first position
  test_pos <- test_data[[group_samples[1]]]$pos[1]
  test_chr <- test_data[[group_samples[1]]]$chr[1]

  combined <- rbindlist(lapply(group_samples, function(s) {
    test_data[[s]][chr == test_chr & pos == test_pos]
  }))
  expected_rate <- sum(combined$mc) / sum(combined$cov)

  bw_match <- bw[GenomicRanges::seqnames(bw) == test_chr &
                    GenomicRanges::start(bw) == test_pos]

  if (length(bw_match) > 0) {
    expect_equal(bw_match$score, expected_rate, tolerance = 1e-6)
  }

  file.remove(result_path)
})


test_that("aggregated cov BigWig values are summed coverage", {
  test_data <- make_bigwig_test_data()
  out_dir <- tempdir()
  metadata <- attr(test_data, "sample_metadata")
  group_samples <- metadata[group_id == "GroupA", sample_id]

  result_path <- create_bigwig(
    meth_data = test_data,
    out = out_dir,
    genome = "hg38",
    aggregate_replicates = TRUE,
    group_name = "GroupA",
    type = "cov"
  )

  bw <- rtracklayer::import.bw(result_path)

  test_pos <- test_data[[group_samples[1]]]$pos[1]
  test_chr <- test_data[[group_samples[1]]]$chr[1]

  combined <- rbindlist(lapply(group_samples, function(s) {
    test_data[[s]][chr == test_chr & pos == test_pos]
  }))
  expected_cov <- sum(combined$cov)

  bw_match <- bw[GenomicRanges::seqnames(bw) == test_chr &
                    GenomicRanges::start(bw) == test_pos]

  if (length(bw_match) > 0) {
    expect_equal(bw_match$score, expected_cov, tolerance = 1e-6)
  }

  file.remove(result_path)
})


# ── All groups export test ───────────────────────────────────────────────────

test_that("create_bigwig works for all groups in a dataset", {
  test_data <- make_bigwig_test_data()
  out_dir <- tempdir()
  metadata <- attr(test_data, "sample_metadata")

  created_files <- character()

  for (group in unique(metadata$group_id)) {
    result_path <- create_bigwig(
      meth_data = test_data,
      out = out_dir,
      genome = "hg38",
      aggregate_replicates = TRUE,
      group_name = group,
      type = "meth"
    )
    expect_true(file.exists(result_path))
    created_files <- c(created_files, result_path)
  }

  # Clean up
  file.remove(created_files)
})


# ── Validation with real test data (if available) ────────────────────────────

test_that("create_bigwig works with real test data", {
  test_data_dir <- "test-data"
  sample_sheet <- file.path(test_data_dir, "M-Series_batches.csv")

  skip_if_not(
    file.exists(sample_sheet),
    "Real test data not available (test-data/M-Series_batches.csv)"
  )

  all_samples <- load_data(
    dir_path = test_data_dir,
    sample_sheet = sample_sheet
  )

  # Subset to small region
  test_samples <- lapply(all_samples, function(dt) {
    dt_copy <- copy(dt[chr == "1" & pos < 1000000])
    dt_copy[, chr := paste0("chr", chr)]
    dt_copy[, uniqueID := paste(chr, pos, site, sep = "_")]
    return(dt_copy)
  })
  attr(test_samples, "sample_metadata") <- attr(all_samples, "sample_metadata")

  out_dir <- tempdir()

  result_path <- create_bigwig(
    meth_data = test_samples,
    out = out_dir,
    genome = "hg38",
    sample_name = names(test_samples)[1],
    type = "meth"
  )

  expect_true(file.exists(result_path))

  bw <- rtracklayer::import.bw(result_path)
  expect_gt(length(bw), 0)

  file.remove(result_path)
})
