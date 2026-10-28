#' Normalize Methylation Rates
#'
#' Normalizes methylation rates using quantile-based normalization to account for
#' non-linear enzyme efficiency differences.
#'
#' @param data_list List of normalized coverage data
#' @param group_names Optional character vector of group names. If NULL (default),
#'   will be extracted from the metadata attached to data_list.
#' @param within_groups Logical indicating whether to normalize within groups (default = TRUE)
#' @param between_groups Logical indicating whether to normalize between groups (default = TRUE)
#' @param within_alpha Numeric between 0 and 1. Shrinkage of site-level
#'   deviations in the within-group step (default = 1, no shrinkage). It does
#'   not control how much efficiency correction is applied; see Details.
#' @param between_alpha Numeric between 0 and 1. Shrinkage of site-level
#'   deviations in the between-group step (default = 1, no shrinkage). It does
#'   not control how much of the group difference is removed; see Details.
#' @param sites_per_quantile Target number of sites per quantile (default = 1000).
#'        For large datasets (>100k sites), this provides a good balance between
#'        normalization granularity and statistical power.
#' @param max_quantiles Maximum number of quantiles to use (default = 50)
#' @param diagnostics Logical indicating whether to print diagnostic information
#'
#' @details This function performs a two-step normalization:
#'   \enumerate{
#'     \item Within-group: the replicates of each group are normalized to one
#'       another, removing per-replicate labeling-efficiency differences.
#'     \item Between-group (if \code{between_groups = TRUE}): all samples of all
#'       groups are normalized together, removing global efficiency differences
#'       between groups.
#'   }
#'
#'   Each step bins sites into equal-count quantiles of their mean rate across
#'   the samples in the set. Within bin \eqn{b}, the rate \eqn{x_{s,i}} of
#'   sample \eqn{s} at site \eqn{i} becomes
#'   \deqn{x'_{s,i} = \bar{m}_b + \alpha (x_{s,i} - m_{s,b})}
#'   where \eqn{m_{s,b}} is sample \eqn{s}'s own mean over the bin's sites and
#'   \eqn{\bar{m}_b} is the average of those means across samples. Results are
#'   clipped to [0, 1].
#'
#'   The two means do different jobs:
#'   \itemize{
#'     \item Replacing \eqn{m_{s,b}} by \eqn{\bar{m}_b} is the normalization. It
#'       shifts each sample, bin by bin, to the common level, and it is applied
#'       in full at every value of alpha. A local difference such as a DMR
#'       survives because its sites are a negligible fraction of any bin.
#'     \item \eqn{\alpha} only scales each site's deviation from its sample's
#'       bin mean. At \eqn{\alpha = 1} those deviations are kept exactly:
#'       \eqn{x' = x + (\bar{m}_b - m_{s,b})}, a per-sample, per-bin shift. This
#'       is not the identity: every sample is still normalized. At
#'       \eqn{\alpha = 0} every sample's sites collapse to the bin's common level.
#'   }
#'
#'   Within a bin of similar rates, the deviations are mostly sampling noise,
#'   so \eqn{\alpha < 1} shrinks noise and makes replicates look more alike
#'   than they are. Downstream tests then overstate precision and call more
#'   false positives; on a real same-condition null (replicates split 2 vs 2),
#'   \code{within_alpha = 0.5} turned 0-4 false-positive DMRs into hundreds.
#'   Values below 1 are therefore not recommended for differential analysis.
#'
#'   The between-group step removes any difference shared by most sites in a
#'   rate stratum, which is what a global efficiency difference looks like. It
#'   cannot tell that from a genuine genome-wide biological difference stratified
#'   by rate, and will remove both. Use \code{between_groups = FALSE} when such a
#'   difference is expected between conditions.
#'
#'   Within-group normalization on its own is not a neutral fallback. It
#'   removes each group's replicate spread but leaves the groups' levels as
#'   they are, so any efficiency difference between groups is tested against
#'   less noise and looks more significant than in the unnormalized data. With
#'   \code{within_groups = TRUE} and \code{between_groups = FALSE}, a warning
#'   is issued when the groups' mean rates differ by more than 5\% of their
#'   average.
#'
#' @return List of lists containing normalized methylation data
#'
#' @export
normalize_methylation_rates <- function(data_list,
                                        group_names = NULL,
                                        within_groups = TRUE,
                                        between_groups = TRUE,
                                        within_alpha = 1,
                                        between_alpha = 1,
                                        sites_per_quantile = 1000,
                                        max_quantiles = 50,
                                        diagnostics = TRUE) {

  # Extract group names from metadata if not provided
  if (is.null(group_names)) {
    # Check for metadata attribute
    if (is.list(data_list) && !is.null(attr(data_list, "sample_metadata"))) {
      sample_metadata <- attr(data_list, "sample_metadata")
      group_names <- unique(sample_metadata$group_id)
      if (diagnostics) {
        message("Using group names extracted from metadata: ",
                paste(group_names, collapse = ", "))
      }
    } else if (is.list(data_list) && all(names(data_list) != "")) {
      # If no metadata but data_list has names, use those as group names
      group_names <- names(data_list)
      if (diagnostics) {
        message("Using list names as group names: ",
                paste(group_names, collapse = ", "))
      }
    } else {
      stop("No group_names provided and unable to extract from metadata")
    }
  }

  # Input validation
  if (!is.list(data_list)) stop("data_list must be a list")
  if (!all(sapply(data_list, is.list))) stop("data_list must be a list of lists")

  result <- data_list

  # Within-group normalization alone removes each group's replicate spread but
  # keeps the groups' levels, so any efficiency difference between groups is
  # tested against less noise and looks MORE significant than in the raw data.
  # On a real same-condition null (M1 split 2 vs 2) this turned 22 false
  # positives into 907. Warn when the groups' levels clearly differ.
  if (within_groups && !between_groups && length(result) > 1) {
    .warn_unequal_group_rates(result)
  }

  # Step 1: Within-group normalization
  if(within_groups) {
    if(diagnostics) cat("\n=== Performing within-group normalization ===\n")
    result <- lapply(names(result), function(group) {
      if(diagnostics) cat(sprintf("\nProcessing group: %s\n", group))
      .normalize_methylation_within_set(
        result[[group]],
        alpha = within_alpha,
        sites_per_quantile = sites_per_quantile,
        max_quantiles = max_quantiles,
        diagnostics = diagnostics
      )
    })
    names(result) <- names(data_list)
  }

  # Step 2: Between-group normalization
  if(between_groups) {
    if(diagnostics) cat("\n=== Performing between-group normalization ===\n")
    # Flatten data for between-group normalization
    all_samples <- unlist(result, recursive = FALSE)

    # The per-bin shift to the common level is applied in full regardless of
    # between_alpha; between_alpha only scales site-level deviations.
    normalized_data <- .normalize_methylation_within_set(
      all_samples,
      alpha = between_alpha,
      sites_per_quantile = sites_per_quantile * 2,  # Fewer quantiles for between-group
      max_quantiles = max_quantiles,
      diagnostics = diagnostics
    )

    # Reconstruct the group structure
    result <- list()
    current_idx <- 1
    for(group in names(data_list)) {
      n_samples <- length(data_list[[group]])
      result[[group]] <- normalized_data[current_idx:(current_idx + n_samples - 1)]
      current_idx <- current_idx + n_samples
    }
  }

  names(result) <- names(data_list)
  return(result)
}
