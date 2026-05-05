#' Find and Filter Shared Sites
#'
#' Identify shared sites across all samples and optionally filters data to only include these sites.
#'
#' @param allc A list of samples returned by load_data(), with each element
#' being a data.table containing a `uniqueID` column.
#' @param filter Logical, whether to filter the data to only include shared sites (default = TRUE).
#'  If FALSE, only returns the vector of shared site IDs.
#' @param quiet Logical, whether to suppress progress messages (default = FALSE).
#'
#' @return If filter=TRUE, returns the filtered list with only shared sites.
#'  If filter=FALSE, returns a vector of shared site identifiers (uniqueID values).
#'
#' @export
find_shared_sites <- function(allc, filter = TRUE, quiet = FALSE) {
  # Input validation
  if (!is.list(allc) || length(allc) == 0) {
    stop("Input must be a non-empty list of samples, as returned by load_data().")
  }
  
  # Check for uid column, create if missing (backward compatibility)
  if (!("uid" %in% colnames(allc[[1]]))) {
    if ("uniqueID" %in% colnames(allc[[1]])) {
      # Legacy data with string uniqueID — fall back to old behavior
      if (!quiet) message("Note: Using legacy string uniqueID. Consider reloading data for faster processing.")
      return(.find_shared_sites_legacy(allc, filter, quiet))
    }
    stop("Each element must be a data.table with a 'uid' column. Reload data with the latest version.")
  }
  
  # Set numeric key for fast operations
  if (filter) {
    for (i in seq_along(allc)) {
      if (!identical(data.table::key(allc[[i]]), "uid")) {
        data.table::setkey(allc[[i]], uid)
        if (!quiet) message("Setting uid key on sample ", i)
      }
    }
  }
  
  if (!quiet) message("Finding shared sites across ", length(allc), " samples...")
  
  if (length(allc) > 20) {
    # Frequency table approach for many samples
    all_uids <- unlist(lapply(allc, function(x) x$uid), use.names = FALSE)
    uid_counts <- tabulate(match(all_uids, unique(all_uids)))
    unique_uids <- unique(all_uids)
    shared_uids <- unique_uids[uid_counts == length(allc)]
  } else {
    # Reduce/intersect for fewer samples — fast on numeric vectors
    shared_uids <- Reduce(intersect, lapply(allc, function(x) x$uid))
  }
  
  if (!quiet) message("Found ", length(shared_uids), " shared sites.")
  
  if (!filter) {
    return(shared_uids)
  }
  
  # Filter each sample using fast data.table key lookup
  if (!quiet) {
    message("Filtering samples to include only shared sites...")
    pb <- utils::txtProgressBar(min = 0, max = length(allc), style = 3)
  }
  
  # Create a keyed lookup table for fast subsetting
  shared_dt <- data.table::data.table(uid = shared_uids, key = "uid")
  
  filtered_data <- vector("list", length(allc))
  names(filtered_data) <- names(allc)
  
  for (i in seq_along(allc)) {
    # Fast keyed semi-join
    filtered_data[[i]] <- allc[[i]][shared_dt, nomatch = NULL]
    if (!quiet) utils::setTxtProgressBar(pb, i)
  }
  
  if (!quiet) {
    close(pb)
    message("Filtering complete.")
  }
  
  # Preserve attributes
  attributes(filtered_data) <- attributes(allc)
  
  return(filtered_data)
}

