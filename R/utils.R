#' @title Internal Utility Functions
#' @description
#' This file contains internal utility functions that are not exported

# This is a special note about parameters that only exist in the main function
# and not in these internal helpers, but are documented for clarity:
# - use_cpp is a parameter in load_data() that controls whether to use C++ implementation


# Data Loading Utilities --------------------------------------------------

#' Internal function to load sample sheet
#' @param sample_sheet Path to sample sheet file or data.frame/data.table
#' @return A data.table containing the sample information
#' @keywords internal
#' @noRd
.load_sample_sheet <- function(sample_sheet) {
  if (is.character(sample_sheet) && length(sample_sheet) == 1) {
    message("Loading sample sheet from: ", sample_sheet)
    sample_data <- data.table::fread(sample_sheet)
  } else if (is.data.frame(sample_sheet)) {
    # Convert to data.table if it's a data.frame
    if (!data.table::is.data.table(sample_sheet)) {
      sample_data <- data.table::as.data.table(sample_sheet)
    } else {
      sample_data <- sample_sheet
    }
  } else {
    stop("sample_sheet must be either a file path or a data.frame/data.table")
  }

  # Validate sample sheet format
  required_cols <- c("group_id", "replicate", "sample_id", "file_name")
  missing_cols <- setdiff(required_cols, colnames(sample_data))
  if (length(missing_cols) > 0) {
    stop("Sample sheet is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  return(sample_data)
}

#' Internal function to filter sample data by groups
#' @param sample_data A data.table containing sample information
#' @param groups Vector of group IDs to include
#' @return Filtered sample_data
#' @keywords internal
#' @noRd
.filter_by_groups <- function(sample_data, groups) {
  if (!is.null(groups)) {
    # Check if requested groups exist in the sample sheet
    missing_groups <- setdiff(groups, unique(sample_data$group_id))
    if (length(missing_groups) > 0) {
      warning("The following requested groups are not in the sample sheet: ",
              paste(missing_groups, collapse = ", "))
    }

    # Filter the sample data
    original_samples <- nrow(sample_data)
    sample_data <- sample_data[group_id %in% groups]

    if (nrow(sample_data) == 0) {
      stop("No samples found for the requested groups: ", paste(groups, collapse = ", "))
    }

    message(sprintf("Filtered to %d samples from %d groups",
                    nrow(sample_data), length(intersect(groups, unique(sample_data$group_id)))))
  }

  return(sample_data)
}

#' Internal function to filter sample data by type
#' @param sample_data A data.table containing sample information
#' @param type character string of methylation type, e.g. GCH or HCG
#' @return Filtered sample_data
#' @keywords internal
#' @noRd
.filter_by_type <- function(sample_data, type) {
  if (!is.null(type)) {
    original_samples <- nrow(sample_data)
    sample_data <- sample_data[grep(type, file_name)]
    if (nrow(sample_data) == 0) {
      stop(sprintf("No %s files found in sample sheet.", type))
    }
    message(sprintf("Filtered sample sheet to %d %s files", nrow(sample_data), type))
  }

  return(sample_data)
}

#' Internal function to validate input files
#'
#' A file named in the sample sheet but missing from the directory is an
#' error: carrying on without it would silently drop a replicate.
#' @param files.list Named character vector of file paths
#' @return files.list, unchanged
#' @keywords internal
#' @noRd
.check_files_exist <- function(files.list) {
  if (length(files.list) == 0) {
    stop("No files to load: the sample sheet has no rows after filtering.",
         call. = FALSE)
  }
  missing_files <- files.list[!file.exists(files.list)]
  if (length(missing_files) > 0) {
    stop(sprintf("%d file(s) named in the sample sheet were not found:\n  %s",
                 length(missing_files), paste(missing_files, collapse = "\n  ")),
         call. = FALSE)
  }
  if (length(files.list) < 2) {
    stop("At least two files are required for analysis.", call. = FALSE)
  }
  return(files.list)
}

# Internal function to load using C++ implementation
.load_data_cpp <- function(files.list, cores, min_coverage = 0L) {
  raw_data <- readMethylationFiles(files.list, min_coverage = min_coverage)
  all_samples <- lapply(raw_data, data.table::as.data.table)
  names(all_samples) <- names(files.list)
  return(all_samples)
}

# Internal function as fallback using pure R
.load_data_r <- function(files.list, cores, min_coverage = 0L) {

  if (cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    # Use parallel processing
    all_samples <- .load_data_parallel(files.list, cores,
                                       min_coverage = min_coverage)
  } else {
    # Use sequential processing with progress bar
    all_samples <- .load_data_sequential(files.list,
                                         min_coverage = min_coverage)
  }

  return(all_samples)
}

#' Internal function to load a single file using R
#'
#' Reads the same records as the C++ reader (\code{readMethylationFile}) and
#' fails in the same cases: a first line that is not a record is taken as a
#' header and skipped, blank lines are skipped, and any other line that is not
#' a valid record is an error, as is a file with no records.
#'
#' @return data.table with columns chr, pos, strand, site, mc, cov
#' @keywords internal
#' @noRd
.load_single_file_r <- function(file, min_coverage = 0L) {
  con <- gzfile(file, "r")
  # A truncated gzip stream surfaces as a warning from readLines(); treat it
  # as an error. A missing final newline is harmless.
  lines <- withCallingHandlers(
    tryCatch(readLines(con), finally = close(con)),
    warning = function(w) {
      if (grepl("incomplete final line", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
      stop("Error reading ", file, ": ", conditionMessage(w),
           " (the file may be truncated or corrupt).", call. = FALSE)
    })

  line_no <- seq_along(lines)
  blank <- !nzchar(trimws(lines))
  lines <- lines[!blank]
  line_no <- line_no[!blank]
  if (length(lines) == 0) {
    stop("No records found in ", file, ".", call. = FALSE)
  }

  fields <- data.table::tstrsplit(lines, "\t", fixed = TRUE, fill = NA_character_)
  if (length(fields) < 6) {
    stop(file, ": expected at least 6 tab-separated columns ",
         "(chr, pos, strand, site, mc, cov), found ", length(fields), ".",
         call. = FALSE)
  }
  as_int <- function(x) {
    out <- rep(NA_integer_, length(x))
    ok <- !is.na(x) & grepl("^[+-]?[0-9]+$", x)
    out[ok] <- suppressWarnings(as.integer(x[ok]))
    out
  }
  pos <- as_int(fields[[2]])
  mc  <- as_int(fields[[5]])
  cov <- as_int(fields[[6]])
  bad <- is.na(pos) | is.na(mc) | is.na(cov) | is.na(fields[[1]]) |
    is.na(fields[[3]]) | is.na(fields[[4]])

  keep <- rep(TRUE, length(lines))
  if (bad[1]) {           # header line
    keep[1] <- FALSE
    bad[1] <- FALSE
  }
  if (any(bad)) {
    stop(sprintf(paste0("%s: line %d is not a valid record (expected chr, pos, ",
                        "strand, site, mc, cov separated by tabs)."),
                 file, line_no[which(bad)[1]]), call. = FALSE)
  }
  wrong <- keep & (mc < 0L | cov < 0L | mc > cov)
  if (any(wrong)) {
    stop(sprintf("%s: line %d has impossible counts (mc = %d, cov = %d).",
                 file, line_no[which(wrong)[1]], mc[which(wrong)[1]],
                 cov[which(wrong)[1]]), call. = FALSE)
  }
  if (!any(keep)) {
    stop("No records found in ", file, ".", call. = FALSE)
  }

  df <- data.table::data.table(
    chr = fields[[1]][keep],
    pos = pos[keep],
    strand = fields[[3]][keep],
    site = fields[[4]][keep],
    mc = mc[keep],
    cov = cov[keep]
  )

  # Apply minimum coverage filter
  if (min_coverage > 0L) {
    df <- df[cov >= min_coverage]
  }

  return(df)
}

#' Internal function to sequentially load files using R implementation
#' @return list of data.tables of methylation data
#' @keywords internal
#' @noRd
.load_data_sequential <- function(files.list, min_coverage = 0L) {
  # Create progress bar
  pb <- progress::progress_bar$new(
    format = "  Loading [:bar] :percent in :elapsed",
    total = length(files.list),
    clear = FALSE
  )

  all_samples <- lapply(files.list, function(file) {
    pb$tick()
    .load_single_file_r(file, min_coverage = min_coverage)
  })

  names(all_samples) <- names(files.list)
  return(all_samples)
}

#' Internal function to load files in parallel using R implementation
#'
#' Falls back to sequential loading only if the cluster itself cannot be set
#' up (for example SMFnorm is not installed where the workers can load it).
#' An error from reading a file is not a reason to fall back: it is raised.
#'
#' @return list of data.tables of methylation data
#' @keywords internal
#' @noRd
.load_data_parallel <- function(files.list, cores, min_coverage = 0L) {
  n_workers <- max(1L, min(as.integer(cores), length(files.list),
                           parallel::detectCores() - 1L, na.rm = TRUE))
  if (n_workers < 2L) {
    return(.load_data_sequential(files.list, min_coverage = min_coverage))
  }

  cl <- tryCatch({
    cl <- parallel::makeCluster(n_workers)
    ok <- unlist(parallel::clusterEvalQ(
      cl, requireNamespace("SMFnorm", quietly = TRUE)))
    if (!all(ok)) {
      parallel::stopCluster(cl)
      stop("SMFnorm could not be loaded on the workers")
    }
    cl
  }, error = function(e) {
    message("Could not start parallel workers (", conditionMessage(e),
            "). Loading files sequentially...")
    NULL
  })
  if (is.null(cl)) {
    return(.load_data_sequential(files.list, min_coverage = min_coverage))
  }
  on.exit(parallel::stopCluster(cl), add = TRUE)

  message(sprintf("Loading %d files on %d workers...", length(files.list),
                  n_workers))
  all_samples <- parallel::parLapply(
    cl, unname(files.list),
    function(f, min_coverage) {
      reader <- utils::getFromNamespace(".load_single_file_r", "SMFnorm")
      reader(f, min_coverage = min_coverage)
    },
    min_coverage = min_coverage)
  names(all_samples) <- names(files.list)
  return(all_samples)
}

#' Internal function to post-process loaded samples
#'
#' Adds the rate and the numeric site identifier. Identifiers are assigned
#' across all samples together so they are comparable between samples whatever
#' the chromosome names are (see \code{.chr_levels}).
#'
#' @return list of data.tables of methylation data
#' @keywords internal
#' @noRd
.post_process_samples <- function(all_samples, sample_data, files.list) {
  chr_levels <- .chr_levels(all_samples)
  if (!is.null(chr_levels)) {
    message(sprintf(paste0("Chromosome names other than 1-22, X, Y and M found ",
                           "(%d names in total); site identifiers are numbered ",
                           "from the names in these samples."),
                    length(chr_levels)))
  }
  for (i in seq_along(all_samples)) {
    df <- all_samples[[i]]
    df[, rate := ifelse(cov > 0, mc / cov, 0)]
    .create_uid(df, chr_levels = chr_levels)
    if (nrow(df) == 0) {
      warning("Sample '", names(all_samples)[i], "' has no sites",
              " after loading (check min_coverage).", call. = FALSE)
    }
  }

  attr(all_samples, "sample_metadata") <-
    sample_data[match(names(files.list), sample_data$sample_id)]
  attr(all_samples, "chr_levels") <- chr_levels

  return(all_samples)
}


# Methylation Normalization Utilities -------------------------------------

#' Warn when within-group-only normalization meets unequal group levels
#'
#' Called by \code{\link{normalize_methylation_rates}} when the within-group
#' step runs without the between-group step. Within-group normalization keeps
#' each group's level but removes its replicate spread, so a difference between
#' group levels (for example a labeling-efficiency difference) is tested
#' against less noise and looks more significant than in the raw data.
#'
#' @param data_list Named list of groups, each a list of sample data frames
#'   with a \code{rate} column.
#' @param tolerance Relative difference between the highest and lowest group
#'   mean rate, as a fraction of their average, above which to warn
#'   (default 0.05).
#' @return Invisibly, the named vector of group mean rates.
#' @keywords internal
.warn_unequal_group_rates <- function(data_list, tolerance = 0.05) {
  group_means <- vapply(data_list, function(group) {
    mean(vapply(group, function(df) mean(df$rate, na.rm = TRUE), numeric(1)))
  }, numeric(1))

  spread <- (max(group_means) - min(group_means)) / mean(group_means)
  if (is.finite(spread) && spread > tolerance) {
    warning(sprintf(paste0(
      "Group mean methylation rates differ by %.0f%% (%s), but only ",
      "within-group rate normalization is being applied. Within-group ",
      "normalization removes replicate spread and keeps this difference, so ",
      "it will look more significant in downstream tests than in the raw ",
      "data. If the difference is technical (e.g. labeling efficiency), set ",
      "rate_between_groups = TRUE (between_groups in ",
      "normalize_methylation_rates()). Keep it FALSE only if the groups are ",
      "expected to differ genome-wide."),
      100 * spread,
      paste(sprintf("%s = %.3f", names(group_means), group_means),
            collapse = ", ")),
      call. = FALSE)
  }
  invisible(group_means)
}

#' Check that all samples hold the same sites in the same row order
#'
#' Normalization works row by row across samples, so row i must be the same
#' site in every sample. Compares row counts, then the site coordinates
#' (\code{chr}, \code{pos} and \code{site} where present, otherwise \code{uid}).
#' Tables with none of these columns can only be checked for row count.
#'
#' @param sample_list List of sample data frames.
#' @return Invisibly \code{TRUE}; stops with an error if samples are not aligned.
#' @keywords internal
#' @noRd
.check_sites_aligned <- function(sample_list) {
  if (length(sample_list) < 2) return(invisible(TRUE))
  ids <- names(sample_list)
  if (is.null(ids)) ids <- as.character(seq_along(sample_list))
  how <- paste0("All samples must contain the same sites in the same order. ",
                "Run find_shared_sites() on the loaded data before splitting ",
                "it into groups, and do not filter or reorder samples ",
                "individually afterwards.")

  n_rows <- vapply(sample_list, nrow, integer(1))
  if (length(unique(n_rows)) > 1) {
    stop("Samples have different numbers of sites (",
         paste(ids, n_rows, sep = " = ", collapse = ", "), "). ", how,
         call. = FALSE)
  }

  shared_cols <- Reduce(intersect, lapply(sample_list, names))
  key_cols <- intersect(c("chr", "pos", "site"), shared_cols)
  if (!all(c("chr", "pos") %in% key_cols)) {
    key_cols <- intersect("uid", shared_cols)
  }

  as_key <- function(x) if (is.numeric(x)) as.numeric(x) else as.character(x)
  ref <- sample_list[[1]]
  for (i in seq_along(sample_list)[-1]) {
    for (col in key_cols) {
      if (!identical(as_key(ref[[col]]), as_key(sample_list[[i]][[col]]))) {
        stop("Sample '", ids[i], "' does not have the same sites in the same ",
             "order as sample '", ids[1], "' (column '", col, "' differs). ",
             how, call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

# Decimal places kept in a site's mean rate, and in bin boundaries, when sites
# are assigned to rate bins. Far finer than any real difference in rate, far
# coarser than floating-point rounding error.
.bin_digits <- 10L

#' Internal function to normalize methylation rates within a set of samples
#'
#' Within each quantile bin of mean rate, sample s's rate at site i becomes
#' \code{avg_quant_rate + alpha * (x - mean(x over the sites in the bin))}: the
#' sample's own bin mean is replaced by the across-sample bin mean (the
#' normalization, applied in full at any alpha), and alpha scales the site-level
#' deviations. See \code{\link{normalize_methylation_rates}}.
#'
#' @param sample_list List of replicate data frames
#' @param alpha Numeric between 0 and 1. Shrinkage of site-level deviations
#'   within each bin (default = 1, no shrinkage)
#' @param sites_per_quantile Target number of sites per quantile (default = 1000)
#' @param max_quantiles Maximum number of quantiles to use (default = 50)
#' @param diagnostics Logical indicating whether to print diagnostic information
#' @keywords internal
.normalize_methylation_within_set <- function(sample_list,
                                             alpha = 1,
                                             sites_per_quantile = 1000,
                                             max_quantiles = 50,
                                             diagnostics = TRUE) {

  if(diagnostics) message("Calculating average rates...")
  
  # All samples must hold the same sites in the same order
  .check_sites_aligned(sample_list)
  avg_rates <- rowMeans(do.call(cbind, lapply(sample_list, function(df) df$rate)))

  # Bin membership must depend on a site's mean rate only. Rates from counts
  # take few distinct values, so many sites share a mean rate and some of those
  # values fall exactly on a bin boundary. Without rounding, differences in the
  # last floating-point digit (from how the rate was computed, or from the
  # platform) decide which side of the boundary each tied site lands on, and
  # that moved normalized rates by up to ~0.015. Rounding makes sites with the
  # same mean rate exactly equal, so they always share a bin.
  avg_rates <- round(avg_rates, .bin_digits)

  # Dynamic quantile calculation
  n_sites <- length(avg_rates)
  n_quantiles <- as.integer(max(5, min(n_sites/sites_per_quantile, max_quantiles)))
  if(diagnostics) message(sprintf("Using %d quantiles for %d sites", n_quantiles, n_sites))

  # Create quantiles
  probs <- seq(0, 1, length.out = n_quantiles + 1)
  breaks <- unique(round(quantile(avg_rates, prob = probs, names = FALSE),
                         .bin_digits))

  if(length(breaks) < 3) {
    warning("Too few unique values for quantile binning. Reducing number of quantiles.")
    n_quantiles <- length(unique(avg_rates)) - 1
    if(n_quantiles < 2) stop("Not enough unique values for quantile normalization")
    probs <- seq(0, 1, length.out = n_quantiles + 1)
    breaks <- unique(round(quantile(avg_rates, prob = probs, names = FALSE),
                           .bin_digits))
  }

  rate_quantiles <- cut(avg_rates, breaks = breaks, include.lowest = TRUE)

  # Check that quantiles are valid and match data size
  if (length(rate_quantiles) != n_sites) {
    stop(sprintf("Quantile vector length (%d) doesn't match site count (%d)",
                 length(rate_quantiles), n_sites))
  }

  # Initialize result list and progress bar
  result <- sample_list
  # Messages go to stderr; the progress bar is drawn only in an interactive
  # session, so logs and rendered documents stay clean.
  show_progress <- diagnostics && interactive()
  if(diagnostics) message("Normalizing by quantile...")
  if(show_progress) {
    pb <- txtProgressBar(min = 0, max = length(levels(rate_quantiles)), style = 3)
  }

  # Normalize each quantile
  for(i in seq_along(levels(rate_quantiles))) {
    quant <- levels(rate_quantiles)[i]
    sites_in_quantile <- which(rate_quantiles == quant)

    # Check for empty quantiles
    if(length(sites_in_quantile) == 0) {
      if(diagnostics) message(sprintf("Quantile %s is empty, skipping", quant))
      if(show_progress) setTxtProgressBar(pb, i)
      next
    }

    # Validate site indices
    max_sites <- min(sapply(result, nrow))
    if(any(sites_in_quantile > max_sites)) {
      problem_sites <- sites_in_quantile[sites_in_quantile > max_sites]
      if(diagnostics) message(sprintf("Removing %d out-of-bounds sites",
                                      length(problem_sites)))
      sites_in_quantile <- sites_in_quantile[sites_in_quantile <= max_sites]

      # Skip if no valid sites remain
      if(length(sites_in_quantile) == 0) {
        if(show_progress) setTxtProgressBar(pb, i)
        next
      }
    }

    # Calculate average rate for this quantile
    quant_rates <- sapply(result, function(df) mean(df$rate[sites_in_quantile]))
    avg_quant_rate <- mean(quant_rates)

    # Calculate and apply scaling factors with structure preservation
    result <- map2(result, quant_rates, function(df, sample_rate) {
      dt <- data.table::copy(df)
      data.table::setDT(dt)

      # Safety check for this specific data table
      if(max(sites_in_quantile) > nrow(dt)) {
        valid_sites <- sites_in_quantile[sites_in_quantile <= nrow(dt)]
        if(length(valid_sites) == 0) {
          return(dt)  # Return unchanged if no valid sites
        }
        sites_in_quantile <- valid_sites
      }

      # Get current rates for these sites
      current_rates <- dt$rate[sites_in_quantile]

      # Preserve structure while normalizing to average
      centered <- current_rates - mean(current_rates)
      new_rates <- avg_quant_rate + alpha * centered
      
      # Clamp rates to [0, 1] and ensure mc does not exceed cov
      new_rates <- pmin(1, pmax(0, new_rates))

      # Update methylation counts and rates
      dt[sites_in_quantile, `:=`(
        mc = pmin(cov, as.integer(round(cov * new_rates))),
        rate = new_rates
      )]

      return(dt)
    })

    if(show_progress) setTxtProgressBar(pb, i)
  }

  if(show_progress) close(pb)
  if(diagnostics) message("Normalization complete!")

  return(result)
}


# Generate bedGraph Functions ---------------------------------------------

#' Process a single methylation file into bedGraph format
#' @param file_path Path to methylation file
#' @param type Type of bedGraph to create ("meth" or "cov")
#' @param out Output directory path
#' @param sample_name Name of the sample (optional).
#' @return Path to created bedGraph file
#' @keywords internal
process_single_file <- function(file_path, type, out, sample_name = NULL,
                                .options = list(use_cpp_impl = TRUE)) {

  # Validate file path
  if (!file.exists(file_path)) {
    stop("Input file does not exist: ", file_path)
  }

  # Load and process file
  data <- load_single_file(file_path)

  process_single_sample(data, type, out, sample_name = sample_name, .options = .options)
}


#' Process a single sample data frame into bedGraph format
#' @param sample_df Data frame containing methylation data
#' @param type Type of bedGraph to create ("meth" or "cov")
#' @param out Output directory path
#' @param group_name Name of group (optional).
#' @param sample_name Name of the sample (optional).
#' @param .options Internal options list containing implementation choices
#' @return Path to created bedGraph file
#' @keywords internal
process_single_sample <- function(sample_df, type, out, group_name = NULL,
                                  sample_name = NULL, .options = list(use_cpp_impl = TRUE)) {

  value_col <- if (type == "meth") "rate" else "cov"
  required_cols <- c("chr", "pos", value_col)
  if (!all(required_cols %in% names(sample_df))) {
    stop("Input data frame missing required columns: ",
         paste(setdiff(required_cols, names(sample_df)), collapse = ", "))
  }
  if (is.null(sample_name)) sample_name <- "sample"

  # Only the columns written, with UCSC-style chromosome names. Both writers
  # used to put "chr" in front of every name, so "chr1" became "chrchr1".
  out_df <- data.frame(
    chr = .ucsc_chr_names(sample_df[["chr"]]),
    pos = as.integer(sample_df[["pos"]]),
    stringsAsFactors = FALSE
  )
  out_df[[value_col]] <- as.numeric(sample_df[[value_col]])
  sample_df <- out_df

  # Select implementation based on options
  if (.options$use_cpp_impl) {
    return(process_single_sample_cpp(sample_df, type, out, group_name, sample_name))
  } else {
    return(process_single_sample_r(sample_df, type, out, group_name, sample_name))
  }
}

#' Process a single sample using C++ implementation
#' @keywords internal
process_single_sample_cpp <- function(sample_df, type, out, group_name = NULL, sample_name = NULL) {
  tryCatch({
    out_file <- createBedgraphCpp(sample_df, type, out, group_name, sample_name)
    if (!file.exists(out_file)) {
      stop("C++ implementation failed to create output file: ", out_file)
    }
    return(out_file)
  }, error = function(e) {
    stop("Error in C++ implementation: ", conditionMessage(e))
  })
}

#' Process a single sample using R implementation
#' @importFrom data.table fwrite
#' @keywords internal
process_single_sample_r <- function(sample_df, type, out, group_name = NULL, sample_name = NULL) {
  # Convert to bedGraph format (your original R code)
  bg_data <- data.table::data.table(
    chr = sample_df$chr,
    start = sample_df$pos - 1L,  # BED format is 0-based
    end = sample_df$pos,
    value = if(type == "meth") sample_df$rate else sample_df$cov
  )

  # Create output filename
  file_prefix <- if(is.null(group_name)) sample_name else paste(group_name, sample_name, sep="_")

  out_file <- file.path(out, paste0(file_prefix, "_", type, ".bedGraph"))

  data.table::fwrite(bg_data, out_file,
                     sep = "\t",
                     quote = FALSE,
                     col.names = FALSE,
                     scipen = 999
  )

  return(out_file)
}


#' UCSC-style chromosome names
#'
#' Adds the "chr" prefix when no name has it (and turns "MT" into "chrM");
#' names that already carry the prefix are returned unchanged. The same rule
#' as \code{.harmonize_chr_names()} uses for bigWig export.
#'
#' @param chr Chromosome names.
#' @return Character vector.
#' @keywords internal
#' @noRd
.ucsc_chr_names <- function(chr) {
  chr <- as.character(chr)
  if (any(startsWith(unique(chr), "chr"))) return(chr)
  chr <- paste0("chr", chr)
  chr[chr == "chrMT"] <- "chrM"
  chr
}


#' Process all samples in a data list into bedGraph format
#' @param data_list List of lists containing methylation data
#' @param type Type of bedGraph to create ("rate" or "cov")
#' @param out Output directory path
#' @return Vector of paths to created bedGraph files
#' @keywords internal
process_all_samples <- function(data_list, type, out,
                                .options = list(use_cpp_impl = TRUE)) {
  out_files <- lapply(names(data_list), function(group) {
    lapply(names(data_list[[group]]), function(sample) {
      process_single_sample(
        sample_df = data_list[[group]][[sample]],
        type = type,
        out = out,
        group_name = group,
        sample_name = sample,
        .options = .options
      )
    })
  })

  return(unlist(out_files))
}


# =============================================================================
# Numeric UID encoding for fast site identification
#
# Replaces the slow paste(chr, pos, site, sep="_") with a numeric key that
# is 10-50x faster to create and uses ~8 bytes per site vs ~40-60 bytes for
# the string version. data.table operations (setkey, merge, intersect) are
# also significantly faster on numeric keys.
#
# Encoding: uid = chr_int * 1e10 + pos * 10 + site_int
#   - chr_int: 1-22 for autosomes, 23=X, 24=Y, 25=M/MT. When any other name is
#              present (chrUn_*, *_random, chrI, chr2L, chr23, scaffolds ...)
#              ALL names are numbered instead by their rank among the names
#              in the data being compared: standard names first in their usual
#              order, then the others alphabetically (see .chr_levels).
#   - pos:     genomic position (up to ~250M for human, fits in 1e10 space)
#   - site_int: 1-9 for methylation context (GCH=1, HCG=2, GCG=3, etc.)
#
# R double precision: exact integers up to 2^53 ~ 9e15, so chr_int can reach
# about 900,000 with no loss of precision.
# =============================================================================

#' Fixed codes for the standard chromosome names
#'
#' @param chr Character vector of chromosome names.
#' @return Integer vector: 1-22, X = 23, Y = 24, M/MT = 25 (with or without a
#'   "chr" prefix), and NA for any other name.
#' @keywords internal
#' @noRd
.standard_chr_code <- function(chr) {
  u <- unique(as.character(chr))
  clean <- sub("^chr", "", u, ignore.case = TRUE)
  code <- rep(NA_integer_, length(u))
  autosome <- clean %in% as.character(1:22)
  code[autosome] <- as.integer(clean[autosome])
  code[clean == "X"] <- 23L
  code[clean == "Y"] <- 24L
  code[clean %in% c("M", "MT")] <- 25L
  code[match(as.character(chr), u)]
}

#' Chromosome names needing their own numbering
#'
#' @param x A table with a \code{chr} column, or a list of such tables.
#' @return NULL when every chromosome name is standard (the fixed codes are
#'   used); otherwise the unique names across all tables (standard names first
#'   in their usual order, then the others alphabetically), whose ranks are
#'   used as codes. Identifiers built from the same levels are comparable.
#' @keywords internal
#' @noRd
.chr_levels <- function(x) {
  if (is.data.frame(x)) x <- list(x)
  names_found <- unique(unlist(lapply(x, function(dt) unique(as.character(dt$chr))),
                               use.names = FALSE))
  code <- .standard_chr_code(names_found)
  if (!anyNA(code)) return(NULL)
  # Standard names first, in their usual order, then the rest alphabetically.
  names_found[order(code, names_found, na.last = TRUE, method = "radix")]
}

#' Chromosome name to integer mapping
#'
#' @param chr Character vector of chromosome names.
#' @param chr_levels NULL for the fixed codes of the standard names, or the
#'   names to number by rank (from \code{.chr_levels}). When NULL and \code{chr}
#'   holds a non-standard name, the names in \code{chr} itself are used.
#' @keywords internal
#' @noRd
.chr_to_int <- function(chr, chr_levels = NULL) {
  if (is.null(chr_levels)) {
    code <- .standard_chr_code(chr)
    if (!anyNA(code)) return(code)
    chr_levels <- .chr_levels(list(list(chr = chr)))
  }
  code <- match(as.character(chr), chr_levels)
  if (anyNA(code)) {
    stop("Chromosome name(s) not in chr_levels: ",
         paste(utils::head(unique(chr[is.na(code)]), 5), collapse = ", "),
         call. = FALSE)
  }
  code
}

#' Site context to integer mapping
#' @keywords internal
#' @noRd
.site_to_int <- function(site) {
  # Use a fixed lookup for known contexts; fallback for unknowns
  lookup <- c(
    "GCH" = 1L, "HCG" = 2L, "GCG" = 3L,
    "CCG" = 4L, "CCC" = 5L, "CCH" = 6L,
    "HCH" = 7L, "GCC" = 8L
  )
  result <- lookup[site]
  # Any unrecognized context gets 9
  result[is.na(result)] <- 9L
  return(as.integer(result))
}


#' Create numeric unique ID for methylation sites
#'
#' Generates a numeric key from chr, pos, and site columns that uniquely
#' identifies each methylation site. This is ~10-50x faster than the string-based
#' \code{paste(chr, pos, site, sep = "_")} and uses significantly less memory.
#'
#' @param dt A data.table with columns \code{chr}, \code{pos}, and \code{site}.
#' @param in_place Logical. If TRUE (default), adds the \code{uid} column to
#'   \code{dt} by reference. If FALSE, returns the uid vector without modifying \code{dt}.
#' @return If \code{in_place = TRUE}, returns \code{dt} invisibly (modified by reference).
#'   If \code{in_place = FALSE}, returns a numeric vector of uid values.
#'
#' @param chr_levels Chromosome numbering to use, from \code{.chr_levels} on
#'   every table whose identifiers will be compared. NULL (default) uses the
#'   fixed codes when all names are standard, else the names in \code{dt}.
#'
#' @details The encoding is: \code{uid = chr_int * 1e10 + pos * 10 + site_int}, where
#'   chr_int maps chromosomes to integers (1-22, X=23, Y=24, M=25, or ranks of
#'   \code{chr_levels}) and site_int maps methylation contexts to integers
#'   (GCH=1, HCG=2, etc.). This fits within R's double-precision integer range
#'   (exact up to 2^53) with no loss of precision.
#'
#' @keywords internal
#' @noRd
.create_uid <- function(dt, in_place = TRUE, chr_levels = NULL) {
  uid_vec <- .chr_to_int(dt$chr, chr_levels) * 1e10 + dt$pos * 10 +
    .site_to_int(dt$site)

  if (in_place) {
    dt[, uid := uid_vec]
    return(invisible(dt))
  } else {
    return(uid_vec)
  }
}


#' Create string uniqueID from components (for RSE export only)
#'
#' This is the legacy string-based ID. Use only when a character identifier
#' is required (e.g., matrix rownames, GRanges metadata).
#'
#' @param dt A data.table with columns \code{chr}, \code{pos}, and \code{site}.
#' @param in_place Logical. If TRUE, adds \code{uniqueID} column by reference.
#' @return If \code{in_place = TRUE}, returns \code{dt} invisibly.
#'   If \code{in_place = FALSE}, returns a character vector.
#'
#' @keywords internal
#' @noRd
.create_string_uid <- function(dt, in_place = TRUE) {
  uid_str <- paste(dt$chr, dt$pos, dt$site, sep = "_")
  
  if (in_place) {
    dt[, uniqueID := uid_str]
    return(invisible(dt))
  } else {
    return(uid_str)
  }
}

