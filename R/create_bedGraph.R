#' Converts normalized data into bedGraph files
#'
#' @param data Input data, one of:
#'        1. A list of groups, each a named list of samples (as returned by
#'           \code{\link{split_by_groups}} or
#'           \code{\link{normalize_methylation_data}})
#'        2. A named list of samples (one group, or the flat list returned by
#'           \code{\link{load_data}})
#'        3. A single data frame containing one sample
#'        4. A character string specifying path to input file
#' @param type Character string for type of value to include in bedGraph.
#'        Should be one of "cov" or "meth" (for coverage bedGraph or methylation rate bedGraph, respectively.)
#' @param all Logical. Should function be performed on all files in object. Default is TRUE.
#'        Ignored if input is a single sample or file.
#' @param out Directory to write output bedGraph files to. Default is current working directory.
#' @param group_name Name of the group (optional), prefixed to the file name.
#'        Used for a single sample or a list of samples; a list of groups uses
#'        its own group names.
#' @param sample_name Name to use for the output file of a single sample or
#'        file (optional). Defaults to "sample" for a data frame and to the
#'        file name without its extension for a file. Lists use their names.
#' @param use_cpp Logical. Use the C++ writer. Default is TRUE. Falls back to
#'        the R implementation if C++ is not available.
#' @details
#' Converts methylation data to bedGraph format. Can handle multiple input types:
#' - If given a normalized data object (list of lists), processes all samples by default
#' - If given a single sample data frame, processes just that sample
#' - If given a file path, loads and processes that file
#'
#' Files are named \code{<group>_<sample>_<type>.bedGraph} (without the group
#' when there is none). Chromosome names are written in UCSC style, as in
#' \code{\link{create_bigwig}}: when no name starts with "chr" the prefix is
#' added and "MT" becomes "chrM"; names that already have the prefix are
#' written unchanged.
#'
#' @return Invisibly, a character vector of the bedGraph file paths created.
#' @export
create_bedgraph <- function(data,
                            type = c("meth", "cov"),
                            all = TRUE,
                            out = ".",
                            group_name = NULL,
                            sample_name = NULL,
                            use_cpp = TRUE) {

  # Match and check type argument
  type <- match.arg(type)

  # Check and create output directory if needed
  out <- normalizePath(out, mustWork = FALSE)
  if (!dir.exists(out)) {
    dir.create(out, recursive = TRUE)
  }

  # Determine C++ implementation availability
  has_rcpp <- requireNamespace("Rcpp", quietly = TRUE) &&
    is.function(get("createBedgraphCpp", envir = asNamespace("SMFnorm")))

  use_cpp_impl <- use_cpp && has_rcpp

  if (use_cpp && !has_rcpp) {
    message("Rcpp or createBedgraphCpp not available. Using R implementation.")
    message("For faster processing, ensure Rcpp is installed and createBedgraphCpp is exported.")
    use_cpp_impl <- FALSE
  }

  # Store implementation choice
  .options <- list(use_cpp_impl = use_cpp_impl)

  # Handle different input types. A data frame is itself a list, so test for
  # it before testing for lists of samples.
  if (is.character(data)) {
    # Input is a file path
    if (length(data) != 1L) {
      stop("A file path input must be a single path.")
    }
    if (!file.exists(data)) {
      stop("Input file does not exist: ", data)
    }
    if (is.null(sample_name)) {
      sample_name <- sub("(\\.(tsv|txt|bed|allc))?(\\.gz)?$", "", basename(data))
    }
    files <- process_single_file(data, type, out, sample_name, .options)

  } else if (is.data.frame(data)) {
    # Input is a single sample
    files <- process_single_sample(sample_df = data,
                                   type = type,
                                   out = out,
                                   group_name = group_name,
                                   sample_name = sample_name,
                                   .options = .options)

  } else if (is.list(data)) {
    if (!all) {
      stop("When 'all=FALSE', input should be a single sample or file path")
    }
    if (length(data) == 0L) stop("'data' is an empty list.")
    is_sample <- vapply(data, is.data.frame, logical(1))
    is_group  <- vapply(data, function(x) {
      !is.data.frame(x) && is.list(x) && length(x) > 0L &&
        all(vapply(x, is.data.frame, logical(1)))
    }, logical(1))

    if (all(is_sample)) {
      # One group: a list of samples, as in one element of split_by_groups()
      # output or the flat list returned by load_data()
      .check_sample_names(data, "'data'")
      files <- unlist(lapply(names(data), function(sample) {
        process_single_sample(data[[sample]], type, out,
                              group_name = group_name,
                              sample_name = sample,
                              .options = .options)
      }))
    } else if (all(is_group)) {
      # Groups of samples
      .check_sample_names(data, "'data'")
      for (group in names(data)) {
        .check_sample_names(data[[group]], paste0("group '", group, "'"))
      }
      files <- process_all_samples(data, type, out, .options)
    } else {
      stop("A list input must be either a list of samples (data frames) or a ",
           "list of groups, each a list of samples.")
    }
  } else {
    stop("Invalid input type. Must be either a file path, data frame, or list of samples")
  }

  invisible(files)
}


#' Stop unless every element of a list of samples or groups has its own name
#' @keywords internal
#' @noRd
.check_sample_names <- function(x, what) {
  nm <- names(x)
  if (is.null(nm) || any(is.na(nm)) || any(nm == "") || anyDuplicated(nm)) {
    stop("Every element of ", what, " needs a unique name; the names become ",
         "the bedGraph file names.")
  }
  invisible(TRUE)
}
