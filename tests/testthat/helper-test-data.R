# tests/testthat/helper-test-data.R
# Helper functions for generating synthetic test data
# This file is automatically sourced by testthat before any test files.

library(data.table)

#' Create a small synthetic methylation sample
#'
#' @param n_sites Number of sites.
#' @param chrs Character vector of chromosome names to use.
#' @param mean_cov Mean coverage (Poisson lambda).
#' @param mean_rate Mean methylation rate (beta distribution centered here).
#' @param seed Random seed.
#' @return data.table mimicking allc-format data.
make_test_sample <- function(n_sites = 100,
                              chrs = c("1", "2", "3"),
                              mean_cov = 20,
                              mean_rate = 0.3,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  dt <- data.table(
    chr  = sample(chrs, n_sites, replace = TRUE),
    pos  = sort(sample(1e6, n_sites)),
    strand = sample(c("+", "-"), n_sites, replace = TRUE),
    site = rep("GCH", n_sites)
  )

  dt[, cov := rpois(n_sites, lambda = mean_cov)]
  dt[cov == 0, cov := 1L]  # no zero coverage
  dt[, mc := rbinom(n_sites, size = cov, prob = mean_rate)]
  dt[, rate := mc / cov]
  dt[, uniqueID := paste(chr, pos, site, sep = "_")]
  .create_uid(dt)

  return(dt)
}


#' Create a list of test samples (mimics load_data output)
#'
#' @param n_samples Number of samples.
#' @param n_sites Number of sites per sample.
#' @param shared_frac Fraction of sites shared across all samples.
#' @param group_ids Character vector of group assignments.
#' @param seed Random seed.
#' @return Named list of data.tables with sample_metadata attribute.
make_test_dataset <- function(n_samples = 4,
                               n_sites = 100,
                               shared_frac = 0.8,
                               group_ids = c("A", "A", "B", "B"),
                               seed = 42) {
  set.seed(seed)

  # Generate a reference set of sites
  n_shared <- floor(n_sites * shared_frac)
  n_unique <- n_sites - n_shared

  ref_dt <- data.table(
    chr  = sample(c("1", "2", "3"), n_shared, replace = TRUE),
    pos  = sort(sample(1e6, n_shared)),
    strand = sample(c("+", "-"), n_shared, replace = TRUE),
    site = rep("GCH", n_shared)
  )
  ref_dt[, uniqueID := paste(chr, pos, site, sep = "_")]
  .create_uid(ref_dt)

  samples <- list()
  sample_ids <- paste0("S", seq_len(n_samples))

  for (i in seq_len(n_samples)) {
    # Start with shared sites
    dt <- copy(ref_dt)

    # Add some unique sites
    if (n_unique > 0) {
      unique_dt <- data.table(
        chr  = sample(c("1", "2", "3"), n_unique, replace = TRUE),
        pos  = sort(sample((1e6 + 1):(2e6), n_unique)),
        strand = sample(c("+", "-"), n_unique, replace = TRUE),
        site = rep("GCH", n_unique)
      )
      unique_dt[, uniqueID := paste(chr, pos, site, sep = "_")]
      .create_uid(unique_dt)
      dt <- rbind(dt, unique_dt)
    }
    
    # Recreate uid after rbind since dt is a new object
    .create_uid(dt)

    # Generate counts with sample-specific variation
    dt[, cov := rpois(.N, lambda = sample(15:30, 1))]
    dt[cov == 0, cov := 1L]
    dt[, mc := rbinom(.N, size = cov, prob = runif(1, 0.2, 0.5))]
    dt[, rate := mc / cov]

    samples[[sample_ids[i]]] <- dt
  }

  # Attach metadata
  metadata <- data.table(
    group_id = group_ids,
    replicate = seq_len(n_samples),
    sample_id = sample_ids,
    file_name = paste0(sample_ids, ".tsv.gz")
  )
  attr(samples, "sample_metadata") <- metadata

  return(samples)
}


#' Create a split (nested) test dataset (mimics split_by_groups output)
#'
#' @param n_reps_per_group Number of replicates per group.
#' @param n_sites Number of sites (all shared).
#' @param group_names Character vector of group names.
#' @param seed Random seed.
#' @return Nested list: group → replicate → data.table.
make_test_split_data <- function(n_reps_per_group = 2,
                                  n_sites = 100,
                                  group_names = c("GroupA", "GroupB"),
                                  seed = 42) {
  set.seed(seed)

  # Reference positions (shared by all)
  ref_dt <- data.table(
    chr  = sample(c("1", "2", "3"), n_sites, replace = TRUE),
    pos  = sort(sample(1e6, n_sites)),
    strand = sample(c("+", "-"), n_sites, replace = TRUE),
    site = rep("GCH", n_sites)
  )
  ref_dt[, uniqueID := paste(chr, pos, site, sep = "_")]
  .create_uid(ref_dt)

  result <- list()
  sample_ids <- character()
  group_ids <- character()

  for (g in group_names) {
    group_list <- list()
    for (r in seq_len(n_reps_per_group)) {
      dt <- copy(ref_dt)
      dt[, cov := rpois(.N, lambda = sample(15:30, 1))]
      dt[cov == 0, cov := 1L]
      dt[, mc := rbinom(.N, size = cov, prob = runif(1, 0.2, 0.5))]
      dt[, rate := mc / cov]

      rep_name <- paste0(g, "_Rep", r)
      group_list[[rep_name]] <- dt
      sample_ids <- c(sample_ids, rep_name)
      group_ids <- c(group_ids, g)
    }
    result[[g]] <- group_list
  }

  metadata <- data.table(
    group_id = group_ids,
    replicate = seq_along(sample_ids),
    sample_id = sample_ids,
    file_name = paste0(sample_ids, ".tsv.gz")
  )
  attr(result, "sample_metadata") <- metadata

  return(result)
}
