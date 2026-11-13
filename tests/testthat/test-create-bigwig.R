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


# ── Depth-scaled coverage tracks ─────────────────────────────────────────────

# Read a coverage BigWig back as a vector aligned to the rows of `dt`.
bw_scores_for <- function(path, dt) {
  bw <- rtracklayer::import.bw(path)
  key_bw <- paste(as.character(GenomicRanges::seqnames(bw)),
                  GenomicRanges::start(bw))
  key_dt <- paste(dt$chr, dt$pos)
  bw$score[match(key_dt, key_bw)]
}

test_that("normalize_coverage divides a sample's coverage by total / mean(totals)", {
  test_data <- make_bigwig_test_data()
  out_dir <- tempfile("bw_"); dir.create(out_dir)
  totals <- vapply(test_data, function(dt) sum(dt$cov), numeric(1))
  expect_gt(diff(range(totals)) / mean(totals), 0.05)   # depths really differ

  for (s in names(test_data)) {
    path <- suppressMessages(create_bigwig(
      meth_data = test_data, out = out_dir, sample_name = s,
      type = "cov", normalize_coverage = TRUE
    ))
    dt <- test_data[[s]]
    got <- bw_scores_for(path, dt)
    expect_false(anyNA(got))
    expect_equal(got, dt$cov / (totals[[s]] / mean(totals)), tolerance = 1e-5)
    # After scaling, every sample's track sums to the same depth.
    expect_equal(sum(got), mean(totals), tolerance = 1e-5)
  }
})

test_that("normalized coverage goes to its own file and leaves the raw track alone", {
  test_data <- make_bigwig_test_data()
  out_dir <- tempfile("bw_"); dir.create(out_dir)
  s <- names(test_data)[1]

  raw <- suppressMessages(create_bigwig(test_data, type = "cov", out = out_dir,
                                        sample_name = s))
  norm <- suppressMessages(create_bigwig(test_data, type = "cov", out = out_dir,
                                         sample_name = s,
                                         normalize_coverage = TRUE))

  expect_match(basename(raw), "_cov\\.bw$")
  expect_match(basename(norm), "_covnorm\\.bw$")
  expect_equal(bw_scores_for(raw, test_data[[s]]), test_data[[s]]$cov,
               tolerance = 1e-6)
})

test_that("normalize_coverage with aggregated replicates scales by group totals", {
  test_data <- make_bigwig_test_data()
  out_dir <- tempfile("bw_"); dir.create(out_dir)
  meta <- attr(test_data, "sample_metadata")
  totals <- vapply(test_data, function(dt) sum(dt$cov), numeric(1))
  group_totals <- tapply(totals[meta$sample_id], meta$group_id, sum)

  path <- suppressMessages(create_bigwig(
    meth_data = test_data, out = out_dir, type = "cov",
    aggregate_replicates = TRUE, group_name = "GroupA",
    normalize_coverage = TRUE
  ))

  ids <- meta[group_id == "GroupA", sample_id]
  summed <- Reduce(`+`, lapply(test_data[ids], `[[`, "cov"))
  got <- bw_scores_for(path, test_data[[ids[1]]])
  expect_equal(got, summed / (group_totals[["GroupA"]] / mean(group_totals)),
               tolerance = 1e-5)
})

test_that("normalize_coverage does not modify the input data", {
  test_data <- make_bigwig_test_data()
  before <- lapply(test_data, function(dt) dt$cov + 0L)
  out_dir <- tempfile("bw_"); dir.create(out_dir)

  suppressMessages(create_bigwig(test_data, type = "cov", out = out_dir,
                                 sample_name = names(test_data)[1],
                                 normalize_coverage = TRUE))

  expect_identical(lapply(test_data, `[[`, "cov"), before)
})

test_that("normalize_coverage is rejected for rate tracks and bad values", {
  test_data <- make_bigwig_test_data()
  s <- names(test_data)[1]

  expect_error(
    create_bigwig(test_data, type = "meth", out = tempdir(), sample_name = s,
                  normalize_coverage = TRUE),
    "only applies to type = \"cov\""
  )
  expect_error(
    create_bigwig(test_data, type = "cov", out = tempdir(), sample_name = s,
                  normalize_coverage = NA),
    "must be TRUE or FALSE"
  )
})

test_that("coverage normalization is no longer part of the normalization API", {
  expect_false("normalize_coverage" %in% getNamespaceExports("SMFnorm"))
  expect_false(any(c("do_coverage_norm", "coverage_between_groups") %in%
                     names(formals(normalize_methylation_data))))
})

test_that("normalize_methylation_data leaves coverage as observed integer depth", {
  d <- make_test_split_data(n_reps_per_group = 2, n_sites = 500, seed = 5)

  out <- suppressMessages(normalize_methylation_data(d, min_coverage = 0,
                                                     diagnostics = FALSE))

  for (g in names(d)) {
    for (s in names(d[[g]])) {
      expect_identical(out[[g]][[s]]$cov, d[[g]][[s]]$cov)
      expect_type(out[[g]][[s]]$mc, "integer")
    }
  }
})
