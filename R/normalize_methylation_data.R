#' Complete Methylation Data Normalization
#'
#' Filters sites by coverage and normalizes methylation rates.
#'
#' @param data_list List of methylation data by group
#' @param group_names Optional character vector of group names. If NULL (default),
#'   will be extracted from the metadata attached to data_list.
#' @param normalize_rates Logical, whether to perform rate normalization
#' @param rate_within_groups Logical, normalize rates within groups
#' @param rate_between_groups Logical, normalize rates between groups
#' @param within_alpha Shrinkage of site-level deviations in the within-group
#'   rate step (default = 1, no shrinkage). The efficiency correction itself is
#'   applied in full at any value; see \code{\link{normalize_methylation_rates}}.
#' @param between_alpha Shrinkage of site-level deviations in the between-group
#'   rate step (default = 1, no shrinkage). See
#'   \code{\link{normalize_methylation_rates}}.
#' @param min_coverage Minimum coverage threshold. A site is kept only if every
#'   sample in every group has at least this coverage there, so all samples keep
#'   the same sites. Use 0 to skip the filter.
#' @param max_quantiles Maximum number of quantiles to use for rate normalization
#' @param sites_per_quantile Target number of sites per quantile for rate normalization
#' @param diagnostics Logical, whether to output diagnostic information
#'
#' @details
#' This function provides a unified interface to the complete methylation normalization
#' workflow: a coverage filter shared by all samples, then methylation rate
#' normalization. Read counts and coverage are not rescaled, so \code{cov}
#' stays the observed sequencing depth. For depth-scaled coverage tracks, see
#' the \code{normalize_coverage} argument of \code{\link{create_bigwig}}.
#'
#' The defaults (\code{within_alpha = 1}, \code{between_alpha = 1}) apply the
#' full efficiency correction without shrinking site-level variation. Alpha
#' values below 1 shrink sampling noise along with everything else and make
#' downstream differential tests overconfident; they are not recommended for
#' differential analysis. See \code{\link{normalize_methylation_rates}} for the
#' formula and what each step removes.
#'
#' \code{\link{normalize_methylation_rates}} can also be called directly.
#'
#' All samples must contain the same sites in the same row order, as returned by
#' \code{\link{find_shared_sites}}; the function stops if they do not.
#'
#' @return Normalized methylation data list with the same group and sample names
#'   as the input, and its \code{sample_metadata} attribute if it had one.
#'
#' @examples
#' \dontrun{
#' # Basic usage with defaults
#' normalized_data <- normalize_methylation_data(
#'   data_list = split_data,
#'   group_names = group_names
#' )
#'
#' # Correct efficiency within groups only, e.g. when the conditions are
#' # expected to differ genome-wide and the between-group step would remove
#' # that difference along with any efficiency difference
#' normalized_data <- normalize_methylation_data(
#'   data_list = split_data,
#'   group_names = group_names,
#'   rate_between_groups = FALSE
#' )
#' }
#'
#' @seealso
#' \code{\link{normalize_methylation_rates}} for rate normalization only
#'
#' @export
#'
normalize_methylation_data <- function(data_list,
                                       group_names = NULL,
                                       normalize_rates = TRUE,
                                       rate_within_groups = TRUE,
                                       rate_between_groups = TRUE,
                                       within_alpha = 1,
                                       between_alpha = 1,
                                       min_coverage = 5,
                                       max_quantiles = 50,
                                       sites_per_quantile = 1000,
                                       diagnostics = FALSE) {

  group_names <- extract_group_names(data_list, diagnostics)

  # Input validation
  if (!is.list(data_list)) stop("data_list must be a list")
  if (!all(vapply(data_list, is.list, logical(1)))) {
    stop("data_list must be a list of lists")
  }
  sample_metadata <- attr(data_list, "sample_metadata")
  group_ids <- names(data_list)

  # Every later step works row by row across samples, so all samples must hold
  # the same sites in the same order.
  .check_sites_aligned(unlist(unname(data_list), recursive = FALSE))

  # Step 0: Filter by minimum coverage. A site is kept only if EVERY sample
  # reaches min_coverage there, so all samples keep the same sites. Filtering
  # each sample on its own coverage would leave them with different rows.
  if (min_coverage > 0) {
    keep <- Reduce(`&`, lapply(unlist(unname(data_list), recursive = FALSE),
                               function(df) df$cov >= min_coverage))
    if (!any(keep)) {
      stop(sprintf("No sites have coverage >= %s in every sample.",
                   format(min_coverage)), call. = FALSE)
    }
    if (diagnostics) {
      message(sprintf(
        "Coverage filter: kept %d of %d sites with coverage >= %s in every sample",
        sum(keep), length(keep), format(min_coverage)))
    }
    if (!all(keep)) {
      data_list <- lapply(data_list, function(group) {
        lapply(group, function(df) df[keep, ])
      })
    }
  }

  # Step 1: Rate normalization
  if (normalize_rates) {
    data_list <- normalize_methylation_rates(data_list,
                                             group_names = group_names,
                                             within_groups = rate_within_groups,
                                             between_groups = rate_between_groups,
                                             within_alpha = within_alpha,
                                             between_alpha = between_alpha,
                                             max_quantiles = max_quantiles,
                                             sites_per_quantile = sites_per_quantile,
                                             diagnostics = diagnostics)
  }

  names(data_list) <- group_ids
  attr(data_list, "sample_metadata") <- sample_metadata
  return(data_list)
}
