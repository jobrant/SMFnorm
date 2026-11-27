#' Complete Methylation Data Normalization
#'
#' Filters sites by coverage and normalizes methylation rates.
#'
#' @param data_list Methylation data: a list of groups, each a list of samples
#'   (as returned by \code{\link{split_by_groups}}); a
#'   \code{RangedSummarizedExperiment} with assays \code{coverage} and
#'   \code{rate} (as made by \code{\link{create_RSE}}); or a \code{BSseq}
#'   object. See "Objects" below.
#' @param group_col For an object: the name of the \code{colData} column that
#'   gives each sample's group. Default \code{"group"}.
#' @param group_names Optional character vector naming the groups to
#'   normalize. Only these groups are used and returned, in the order given;
#'   the others play no part in the coverage filter or the between-group step.
#'   If NULL (default), all groups in \code{data_list} are used. Naming a group
#'   that is not in \code{data_list} is an error.
#' @param normalize_rates Logical, whether to perform rate normalization
#' @param rate_within_groups Logical, normalize rates within groups
#' @param rate_between_groups Logical, normalize rates between groups. There
#'   is no default: with two or more groups it must be given, because neither
#'   setting is safe to apply without thought. See "Choosing
#'   \code{rate_between_groups}" below. Ignored for a single group.
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
#' @section Choosing \code{rate_between_groups}:
#' The within-group step makes the replicates of each group agree. The
#' between-group step then brings all groups to a common level in every rate
#' bin. Whether to run it depends on why the groups differ, which the data
#' alone do not show:
#' \itemize{
#'   \item \code{TRUE} removes a difference in labeling efficiency between the
#'     groups. It equally removes a genuine biological difference that is
#'     broad (shared by a large part of the genome), and pushes the opposite
#'     shift onto the remaining sites. Local differences such as individual
#'     footprints or differentially accessible regions are kept.
#'   \item \code{FALSE} keeps every difference between the groups, technical
#'     or biological. Because the replicates have been made to agree, a
#'     technical difference between the groups is then tested against less
#'     noise and gives more false positives than the unnormalized data would.
#' }
#' Use \code{TRUE} when the groups are expected to have similar overall
#' accessibility and a difference in their levels is more likely technical
#' (different batches, enzyme lots or labeling times). Use \code{FALSE} when
#' a broad biological difference is expected, and then treat a difference in
#' group levels with care. \code{plot_efficiency(data, reference = "all")}
#' shows the difference between the groups before the choice is made.
#'
#' @details
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
#' @return For a list: the normalized list with the same group and sample names
#'   as the input (or the groups in \code{group_names}), and its
#'   \code{sample_metadata} attribute if it had one. For an object: an object
#'   of the same class; see "Objects".
#'
#' @section Objects:
#' A \code{RangedSummarizedExperiment} or \code{BSseq} object is normalized by
#' the same code as a list and returned as the same class, restricted to the
#' sites that pass the coverage filter (and to the groups in
#' \code{group_names}, if given). Rows keep their order.
#' \itemize{
#'   \item \code{RangedSummarizedExperiment}: the \code{rate} assay holds the
#'     normalized rates; \code{coverage} is unchanged.
#'   \item \code{BSseq}: the \code{M} assay holds methylated counts consistent
#'     with the normalized rates, \code{round(Cov * rate)}; \code{Cov} is
#'     unchanged. Smoothed values, if present, are left as they were and
#'     should be recomputed.
#' }
#' Sites with a missing value in any sample are dropped. A site with no reads
#' in a sample has rate 0 there; set \code{min_coverage} to at least 1 to
#' leave such sites out. The settings used are recorded in
#' \code{metadata(x)$SMFnorm}. On-disk (HDF5-backed) assays are read into
#' memory, and the result is held in memory.
#'
#' @examples
#' \dontrun{
#' # Replicates made to agree within each group; group levels left as they are
#' normalized_data <- normalize_methylation_data(
#'   data_list = split_data,
#'   rate_between_groups = FALSE
#' )
#'
#' # Two of the groups only
#' normalized_data <- normalize_methylation_data(
#'   data_list = split_data,
#'   group_names = c("Control", "Treated"),
#'   rate_between_groups = FALSE
#' )
#'
#' # Also bring the groups to a common level, when a difference between them
#' # is expected to be technical
#' normalized_data <- normalize_methylation_data(
#'   data_list = split_data,
#'   rate_between_groups = TRUE
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
                                       rate_between_groups,
                                       within_alpha = 1,
                                       between_alpha = 1,
                                       min_coverage = 5,
                                       max_quantiles = 50,
                                       sites_per_quantile = 1000,
                                       diagnostics = FALSE,
                                       group_col = "group") {

  # A SummarizedExperiment-type object: convert, normalize as a list, write back
  if (.is_se(data_list)) {
    args <- list(normalize_rates = normalize_rates,
                 rate_within_groups = rate_within_groups,
                 within_alpha = within_alpha, between_alpha = between_alpha,
                 max_quantiles = max_quantiles,
                 sites_per_quantile = sites_per_quantile)
    if (!missing(rate_between_groups)) {
      args$rate_between_groups <- rate_between_groups
    }
    return(.normalize_se(
      data_list, group_col, group_names, min_coverage,
      normalize = function(lst) {
        do.call(normalize_methylation_data,
                c(list(data_list = lst, min_coverage = 0,
                       diagnostics = diagnostics), args))
      },
      settings = c(args, list(min_coverage = min_coverage))))
  }

  # Input validation, and the groups to work on. `group_names` used to be
  # overwritten here and had no effect.
  data_list <- .select_groups(data_list, group_names,
                             diagnostics = diagnostics && !normalize_rates)
  if (normalize_rates) {
    rate_between_groups <- .require_between_choice(
      if (missing(rate_between_groups)) NULL else rate_between_groups,
      "rate_between_groups", length(data_list))
  } else if (missing(rate_between_groups)) {
    rate_between_groups <- FALSE
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
