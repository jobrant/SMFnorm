#' Groups to normalize
#'
#' Checks that \code{data_list} is a named list of groups and returns it
#' restricted to \code{group_names}, in the order asked for.
#'
#' @param data_list Named list of groups, each a list of samples.
#' @param group_names Groups to keep, or NULL for all of them.
#' @param diagnostics Whether to report the groups used.
#' @return \code{data_list} with only the requested groups. The
#'   \code{sample_metadata} attribute, if any, is kept and restricted to those
#'   groups.
#' @keywords internal
#' @noRd
.select_groups <- function(data_list, group_names = NULL, diagnostics = FALSE) {
  if (!is.list(data_list) || is.data.frame(data_list)) {
    stop("data_list must be a list")
  }
  is_group <- vapply(data_list, function(x) is.list(x) && !is.data.frame(x),
                     logical(1))
  if (!all(is_group)) {
    stop("data_list must be a list of lists: groups, each holding its samples. ",
         "Use split_by_groups() on the output of load_data().")
  }
  available <- names(data_list)
  if (is.null(available) || any(is.na(available)) || any(available == "") ||
      anyDuplicated(available)) {
    stop("data_list must have a unique name for every group")
  }

  sample_metadata <- attr(data_list, "sample_metadata")
  if (!is.null(group_names)) {
    if (!is.character(group_names) || anyNA(group_names) ||
        anyDuplicated(group_names) || length(group_names) == 0L) {
      stop("group_names must be a character vector of distinct group names")
    }
    missing_groups <- setdiff(group_names, available)
    if (length(missing_groups) > 0) {
      stop("group_names not found in data_list: ",
           paste(missing_groups, collapse = ", "),
           ". Available groups: ", paste(available, collapse = ", "))
    }
    data_list <- data_list[group_names]
    if (!is.null(sample_metadata) && "group_id" %in% names(sample_metadata)) {
      sample_metadata <- sample_metadata[sample_metadata$group_id %in% group_names, ]
    }
    attr(data_list, "sample_metadata") <- sample_metadata
  }
  if (diagnostics) {
    message("Groups: ", paste(names(data_list), collapse = ", "))
  }
  data_list
}


#' The between-group step must be chosen, not defaulted
#'
#' @param value The argument as given, or NULL when it was missing.
#' @param arg Its name, for the message.
#' @param n_groups Number of groups being normalized.
#' @return TRUE or FALSE. A single group has nothing to normalize between, so
#'   the step is off whatever was given.
#' @keywords internal
#' @noRd
.require_between_choice <- function(value, arg, n_groups) {
  if (n_groups < 2L) return(FALSE)
  if (is.null(value)) {
    stop("'", arg, "' must be set to TRUE or FALSE; it has no default.\n",
         "  TRUE  brings all groups to a common level in every rate bin. ",
         "This removes a\n        labeling-efficiency difference between ",
         "groups, and equally a genuine broad\n        biological difference. ",
         "Local differences (footprints, DMRs) are kept.\n",
         "  FALSE keeps every difference between groups. A technical ",
         "difference between\n        groups then gives more false positives ",
         "than unnormalized data.\n",
         "See ?normalize_methylation_data, and plot_efficiency(data, ",
         "reference = \"all\") to look at\nthe difference between your groups.",
         call. = FALSE)
  }
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("'", arg, "' must be TRUE or FALSE.", call. = FALSE)
  }
  value
}
