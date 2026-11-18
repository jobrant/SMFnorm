#' Find and Filter Shared Sites
#'
#' Identify shared sites across all samples and optionally filters data to only include these sites.
#'
#' @param allc A list of samples returned by load_data(): data.tables with
#'   columns `chr`, `pos` and `site` (or, failing that, a `uid` or legacy
#'   `uniqueID` column).
#' @param filter Logical, whether to filter the data to only include shared sites (default = TRUE).
#'  If FALSE, only returns the vector of shared site IDs.
#' @param quiet Logical, whether to suppress progress messages (default = FALSE).
#'
#' @return If filter=TRUE, a new list holding only the shared sites, with every
#'  sample in the same row order (sorted by `uid`). The input tables are not
#'  modified. If filter=FALSE, the vector of shared site identifiers.
#'
#' @details Sites are matched on chromosome, position and context. Chromosome
#'  names are matched as written, except that the standard names (1-22, X, Y,
#'  M) match with or without a "chr" prefix. A sample containing the same site
#'  twice is an error.
#'
#' @export
find_shared_sites <- function(allc, filter = TRUE, quiet = FALSE) {
  # Input validation
  if (!is.list(allc) || is.data.frame(allc) || length(allc) == 0) {
    stop("Input must be a non-empty list of samples, as returned by load_data().")
  }
  ids <- names(allc)
  if (is.null(ids)) ids <- as.character(seq_along(allc))

  # Site keys, comparable across samples. Built from chr/pos/site with one
  # chromosome numbering for all samples; an existing uid column is not
  # trusted when chromosome names are non-standard, because it may have been
  # numbered from another set of names.
  has <- function(cols) all(vapply(allc, function(x) all(cols %in% names(x)),
                                   logical(1)))
  chr_levels <- NULL
  if (has(c("chr", "pos", "site"))) {
    chr_levels <- .chr_levels(allc)
    keys <- lapply(allc, function(x) {
      if (is.null(chr_levels) && "uid" %in% names(x) && !anyNA(x$uid)) {
        x$uid
      } else {
        .create_uid(x, in_place = FALSE, chr_levels = chr_levels)
      }
    })
    key_is_uid <- TRUE
  } else if (has("uid")) {
    keys <- lapply(allc, function(x) x$uid)
    key_is_uid <- TRUE
  } else if (has("uniqueID")) {
    if (!quiet) message("Note: using the legacy string uniqueID.")
    keys <- lapply(allc, function(x) x$uniqueID)
    key_is_uid <- FALSE
  } else {
    stop("Each sample must have 'chr', 'pos' and 'site' columns, or a 'uid' ",
         "column. Reload the data with load_data().")
  }

  for (i in seq_along(keys)) {
    if (anyNA(keys[[i]])) {
      stop("Sample '", ids[i], "' has sites with a missing identifier ",
           "(missing chr, pos or site).", call. = FALSE)
    }
    dup <- anyDuplicated(keys[[i]])
    if (dup > 0) {
      stop("Sample '", ids[i], "' contains the same site more than once ",
           "(first repeat at row ", dup, ").", call. = FALSE)
    }
  }

  if (!quiet) message("Finding shared sites across ", length(allc), " samples...")
  shared <- sort(Reduce(intersect, keys))
  if (!quiet) message("Found ", length(shared), " shared sites.")
  if (length(shared) == 0) {
    warning("No sites are shared by all samples. Check that the samples use ",
            "the same chromosome naming and genome build.", call. = FALSE)
  }

  if (!filter) {
    return(shared)
  }

  if (!quiet) message("Filtering samples to include only shared sites...")
  filtered_data <- vector("list", length(allc))
  for (i in seq_along(allc)) {
    # Row subset: a new table, in the same (sorted) site order for every sample.
    res <- data.table::as.data.table(allc[[i]])[match(shared, keys[[i]])]
    if (key_is_uid) {
      res[, uid := shared]
      data.table::setkey(res, uid)
    }
    filtered_data[[i]] <- res
  }
  if (!quiet) message("Filtering complete.")

  # Preserve attributes
  attributes(filtered_data) <- attributes(allc)
  attr(filtered_data, "chr_levels") <- chr_levels

  return(filtered_data)
}
