# Normalizing SummarizedExperiment-type objects
#
# normalize_methylation_data() and normalize_methylation_rates() accept, as
# well as the nested list, a RangedSummarizedExperiment made by create_RSE()
# (assays "coverage" and "rate") or a BSseq object (assays "M" and "Cov"), and
# return an object of the same class. The object is turned into the nested
# list, normalized by the same code, and the result written back.


#' Is this a SummarizedExperiment-type object (including BSseq)?
#' @keywords internal
#' @noRd
.is_se <- function(x) methods::is(x, "SummarizedExperiment")


#' Coverage, rate and (if stored) methylated-count matrices of an object
#'
#' @param x A SummarizedExperiment-type object.
#' @return List with `kind` ("bsseq" or "rse"), and in-memory numeric matrices
#'   `cov` and `rate`. On-disk (DelayedArray) assays are read into memory.
#' @keywords internal
#' @noRd
.se_matrices <- function(x) {
  nms <- SummarizedExperiment::assayNames(x)
  get_assay <- function(name) {
    m <- as.matrix(SummarizedExperiment::assay(x, name, withDimnames = FALSE))
    storage.mode(m) <- "double"
    m
  }
  if (methods::is(x, "BSseq") || all(c("M", "Cov") %in% nms)) {
    if (!all(c("M", "Cov") %in% nms)) {
      stop("A BSseq object needs assays 'M' and 'Cov'; found: ",
           paste(nms, collapse = ", "), call. = FALSE)
    }
    cov <- get_assay("Cov")
    m <- get_assay("M")
    rate <- m / cov
    rate[!is.na(cov) & cov == 0] <- 0      # as load_data(): no reads, rate 0
    return(list(kind = "bsseq", cov = cov, rate = rate))
  }
  if (all(c("coverage", "rate") %in% nms)) {
    return(list(kind = "rse", cov = get_assay("coverage"), rate = get_assay("rate")))
  }
  stop("The object needs assays 'coverage' and 'rate' (as made by create_RSE()) ",
       "or 'M' and 'Cov' (a BSseq object); found: ",
       if (length(nms)) paste(nms, collapse = ", ") else "none", call. = FALSE)
}


#' Group and sample of every column of an object
#'
#' @param x A SummarizedExperiment-type object.
#' @param group_col Name of the colData column holding each sample's group.
#' @return data.frame with `group` and `sample`, one row per column of `x`.
#' @keywords internal
#' @noRd
.se_samples <- function(x, group_col = "group") {
  cd <- SummarizedExperiment::colData(x)
  if (!is.character(group_col) || length(group_col) != 1L || is.na(group_col)) {
    stop("'group_col' must be the name of one colData column.", call. = FALSE)
  }
  if (!group_col %in% names(cd)) {
    stop("colData has no column '", group_col, "' giving each sample's group. ",
         "Columns available: ",
         if (length(names(cd))) paste(names(cd), collapse = ", ") else "none",
         ". Name the column with 'group_col', or add one, e.g. ",
         "colData(x)$group <- c(\"A\", \"A\", \"B\", \"B\").", call. = FALSE)
  }
  group <- as.character(cd[[group_col]])
  if (anyNA(group) || any(group == "")) {
    stop("colData column '", group_col, "' has missing group names.", call. = FALSE)
  }
  sample <- if ("sample" %in% names(cd) && group_col == "group" &&
                !anyNA(cd[["sample"]])) {
    as.character(cd[["sample"]])           # as stored by create_RSE()
  } else {
    colnames(x)
  }
  if (is.null(sample) || anyNA(sample) || any(sample == "")) {
    stop("The object's columns need names (sample identifiers).", call. = FALSE)
  }
  if (anyDuplicated(paste(group, sample, sep = "\r"))) {
    stop("Sample names are not unique within their groups.", call. = FALSE)
  }
  data.frame(group = group, sample = sample, stringsAsFactors = FALSE)
}


#' Nested group/sample list from an object's matrices
#'
#' @param x The object (for its row coordinates).
#' @param mats List from `.se_matrices()`.
#' @param samples data.frame from `.se_samples()`.
#' @return Nested list (group -> sample -> data.table) with a
#'   `sample_metadata` attribute, as split_by_groups() returns.
#' @keywords internal
#' @noRd
.se_group_list <- function(x, mats, samples) {
  n <- nrow(mats$cov)
  if (methods::is(x, "RangedSummarizedExperiment")) {
    rr <- SummarizedExperiment::rowRanges(x)
    coords <- data.table::data.table(
      chr = as.character(GenomicRanges::seqnames(rr)),
      pos = as.integer(GenomicRanges::start(rr)),
      strand = as.character(GenomicRanges::strand(rr)))
    site <- S4Vectors::mcols(rr)[["site"]]
    if (!is.null(site)) coords[, site := as.character(site)]
  } else {
    coords <- data.table::data.table(chr = "row", pos = seq_len(n))
  }

  groups <- unique(samples$group)
  out <- lapply(groups, function(g) {
    cols <- which(samples$group == g)
    lst <- lapply(cols, function(j) {
      dt <- data.table::copy(coords)
      dt[, cov := mats$cov[, j]]
      dt[, rate := mats$rate[, j]]
      dt[, mc := pmin(cov, round(cov * rate))]
      dt
    })
    names(lst) <- samples$sample[cols]
    lst
  })
  names(out) <- groups
  attr(out, "sample_metadata") <- data.table::data.table(
    group_id = samples$group, sample_id = samples$sample)
  out
}


#' Normalize a SummarizedExperiment-type object
#'
#' @param x The object.
#' @param group_col colData column holding the groups.
#' @param group_names Groups to keep, or NULL.
#' @param min_coverage Coverage every sample must reach at a site, or NULL for
#'   no coverage filter.
#' @param normalize Function taking the nested list and returning it normalized.
#' @param settings List of the settings used, stored in the result's metadata.
#' @return `x` restricted to the groups and sites used, with normalized values.
#' @keywords internal
#' @noRd
.normalize_se <- function(x, group_col, group_names, min_coverage, normalize,
                          settings = list()) {
  samples <- .se_samples(x, group_col)
  if (!is.null(group_names)) {
    if (!is.character(group_names) || anyNA(group_names) ||
        anyDuplicated(group_names) || length(group_names) == 0L) {
      stop("group_names must be a character vector of distinct group names",
           call. = FALSE)
    }
    missing_groups <- setdiff(group_names, samples$group)
    if (length(missing_groups)) {
      stop("group_names not found in colData column '", group_col, "': ",
           paste(missing_groups, collapse = ", "), ". Available groups: ",
           paste(unique(samples$group), collapse = ", "), call. = FALSE)
    }
    # Columns in the order of group_names, as the list interface returns them
    keep_cols <- unlist(lapply(group_names, function(g) which(samples$group == g)))
    x <- x[, keep_cols]
    samples <- samples[keep_cols, , drop = FALSE]
  }

  mats <- .se_matrices(x)

  # All samples must hold a value at every site used
  keep <- !(rowSums(is.na(mats$cov)) > 0 | rowSums(is.na(mats$rate)) > 0)
  n_na <- sum(!keep)
  if (n_na > 0) {
    message(sprintf(paste0("Dropping %d of %d sites with a missing coverage or ",
                           "rate in at least one sample."), n_na, length(keep)))
  }
  if (!is.null(min_coverage) && min_coverage > 0) {
    enough <- rowSums(mats$cov >= min_coverage, na.rm = TRUE) == ncol(mats$cov)
    message(sprintf(
      "Coverage filter: kept %d of %d sites with coverage >= %s in every sample",
      sum(keep & enough), sum(keep), format(min_coverage)))
    keep <- keep & enough
  }
  if (!any(keep)) {
    stop("No sites are left: every site has a missing value or too little ",
         "coverage in at least one sample.", call. = FALSE)
  }
  if (!all(keep)) {
    x <- x[keep, ]
    mats$cov <- mats$cov[keep, , drop = FALSE]
    mats$rate <- mats$rate[keep, , drop = FALSE]
  }

  normalized <- normalize(.se_group_list(x, mats, samples))

  flat <- unlist(unname(normalized), recursive = FALSE)
  rate <- vapply(flat, function(dt) as.numeric(dt$rate), numeric(nrow(mats$cov)))
  dimnames(rate) <- NULL
  if (mats$kind == "bsseq") {
    # Counts consistent with the normalized rates; coverage is not changed
    m <- vapply(flat, function(dt) as.numeric(dt$mc), numeric(nrow(mats$cov)))
    dimnames(m) <- NULL
    storage.mode(m) <- "integer"
    SummarizedExperiment::assay(x, "M", withDimnames = FALSE) <- m
    extra <- setdiff(SummarizedExperiment::assayNames(x), c("M", "Cov"))
    if (length(extra)) {
      warning("Assay(s) ", paste(extra, collapse = ", "), " were computed from ",
              "the unnormalized counts and are returned unchanged; recompute ",
              "them (e.g. smoothing) from the normalized object.", call. = FALSE)
    }
  } else {
    SummarizedExperiment::assay(x, "rate", withDimnames = FALSE) <- rate
  }

  md <- S4Vectors::metadata(x)
  md$SMFnorm <- c(settings, list(
    n_sites = nrow(x),
    version = as.character(utils::packageVersion("SMFnorm"))))
  S4Vectors::metadata(x) <- md
  x
}


#' Nested list from a SummarizedExperiment-type object, for plotting
#' @keywords internal
#' @noRd
.se_as_group_list <- function(x, group_col = "group") {
  mats <- .se_matrices(x)
  keep <- !(rowSums(is.na(mats$cov)) > 0 | rowSums(is.na(mats$rate)) > 0)
  if (!all(keep)) {
    x <- x[keep, ]
    mats$cov <- mats$cov[keep, , drop = FALSE]
    mats$rate <- mats$rate[keep, , drop = FALSE]
  }
  .se_group_list(x, mats, .se_samples(x, group_col))
}
