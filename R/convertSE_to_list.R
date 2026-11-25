#' Convert SummarizedExperiment back to a nested list of data tables
#'
#' @description
#' This function performs the inverse operation of \code{create_RSE}. It extracts
#' genomic coordinates and assay data, then reconstructs a hierarchical list
#' from the groups and samples recorded in \code{colData}.
#'
#' @param rse A \code{SummarizedExperiment} or \code{RangedSummarizedExperiment} object.
#' @param verbose Logical. If \code{TRUE} (default), prints progress messages.
#' @param groups Character vector. The names of specific groups to extract.
#' If \code{NULL} (default), all samples are extracted.
#'
#' @details
#' Groups and samples are read from the \code{group} and \code{sample} columns
#' of \code{colData}, which \code{create_RSE} writes. If either column is
#' absent, the column names are split at their first underscore instead.
#'
#' @return A doubly-nested named list of \code{data.table} objects, with groups
#'   and samples named as they were given to \code{create_RSE}.
#'
#' @export

convertSE_to_list <- function(rse, groups = NULL, verbose = TRUE) {

  # Validate input object
  if (!is(rse, "SummarizedExperiment") && !is(rse, "RangedSummarizedExperiment")) {
    stop("Input 'rse' must be a SummarizedExperiment or RangedSummarizedExperiment object.")
  }

  # Step 1: Group and sample of every column. create_RSE() records both in
  # colData. Objects without those columns fall back to splitting the column
  # name at its first underscore, which is only right when group names
  # contain no underscore.
  col_data <- SummarizedExperiment::colData(rse)
  all_col_ids <- rownames(col_data)

  if (all(c("group", "sample") %in% names(col_data))) {
    parsed_samples <- as.character(col_data$group)
    parsed_reps    <- as.character(col_data$sample)
  } else {
    name_split <- stringr::str_split(all_col_ids, pattern = "_", n = 2, simplify = TRUE)
    parsed_samples <- name_split[, 1]
    parsed_reps    <- name_split[, 2]
  }

  # Handle Group Filtering
  if (!is.null(groups)) {
    if (verbose) message("Filtering for specified groups...")
    keep_idx <- which(parsed_samples %in% groups)

    if (length(keep_idx) == 0) {
      stop("None of the specified groups were found. Available groups: ",
           paste(unique(parsed_samples), collapse = ", "))
    }
    missing_groups <- setdiff(groups, parsed_samples)
    if (length(missing_groups) > 0) {
      warning("The following groups are not found in the data: ",
              paste(missing_groups, collapse = ", "))
    }

    # Subset SE or RSE and our parsed vectors
    rse <- rse[, keep_idx]

    # CRITICAL: Drop unused factor levels in colData (e.g., removing Group B labels)
    SummarizedExperiment::colData(rse) <- droplevels(SummarizedExperiment::colData(rse))

    parsed_samples <- parsed_samples[keep_idx]
    parsed_reps    <- parsed_reps[keep_idx]
    all_col_ids    <- all_col_ids[keep_idx]
  }



  # Step 2: Extracting genomic metadata from rowRanges..."
  if (verbose) message("Extracting genomic metadata from rowRanges...")
  # Extract the 'spine' (chr, pos, strand, etc.)
  # mcols() extracts the metadata columns like 'site' and 'uniqueID'
  coords <- as.data.table(SummarizedExperiment::rowRanges(rse))

  # Ensure column names match your original requirements
  # GRanges 'seqnames' becomes 'chr', 'start' becomes 'pos'
  setnames(coords, old = c("seqnames", "start"), new = c("chr", "pos"))
  coords[, chr := as.character(chr)]
  coords[, strand := as.character(strand)]

  # Drop IRanges specific columns like 'end' and 'width' if the width is 1 throughout
  if(length(unique(coords$width))==1){
    coords <- coords[, .(chr, pos, strand, site, uniqueID)]
  }

  # Extract assay data
  if (verbose) message("Extracting assay data...")
  cov_arr <- SummarizedExperiment::assay(rse, "coverage")
  rate_arr <- SummarizedExperiment::assay(rse, "rate")

  # One chromosome numbering for every sample (see find_shared_sites())
  chr_levels <- .chr_levels(list(coords))

  # Create a flat list of data.tables
  flat_list <- lapply(seq_len(ncol(rse)), function(i) {
    # Combine coordinates with measurements for this specific column
    dt <- data.table::copy(coords)
    dt[, cov := cov_arr[, i]]
    dt[, rate := rate_arr[, i]]

    # Methylated counts are whole reads, never more than the coverage
    dt[, mc := pmin(cov, round(cov * rate))]

    .create_uid(dt, chr_levels = chr_levels)
    dt
  })

  names(flat_list) <- colnames(rse)

  if (verbose) message ("Nesting into Group/Sample hierarchy...")

  # Split by Group (Tier 1)
  unique_samples <- unique(parsed_samples)

  lol <- lapply(unique_samples, function(s_name) {
    sample_idx <- which(parsed_samples == s_name)
    sample_reps <- flat_list[sample_idx]
    names(sample_reps) <- parsed_reps[sample_idx]
    return(sample_reps)
  })

  names(lol) <- unique_samples

  if (verbose) message("Conversion complete")

  return(lol)
}
