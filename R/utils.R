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
#' @param files.list list or character vector of file names
#' @return list or character vector of file names
#' @keywords internal
#' @noRd
.check_files_exist <- function(files.list) {
  existing_files <- file.exists(files.list)
  if (!all(existing_files)) {
    warning(sprintf("%d files from sample sheet not found in directory",
                    sum(!existing_files)))
    message("Missing files: ", paste(files.list[!existing_files], collapse=", "))
    files.list <- files.list[existing_files]
  }

  # Validation step
  if (length(files.list) == 0) {
    stop("No files found in the directory matching the sample sheet.")
  }
  if (length(files.list) < 2) {
    stop("At least two files are required for analysis.")
  }

  message(sprintf("Loading %d files...", length(files.list)))

  return(files.list)
}

# Internal function to load using C++ implementation
.load_data_cpp <- function(files.list, cores) {
  # Call the C++ function
  raw_data <- readMethylationFiles(files.list)
  names(raw_data) <- names(files.list)

  # Convert to data.tables and add calculated columns
  all_samples <- lapply(seq_along(raw_data), function(i) {
    df <- as.data.table(raw_data[[i]])
    df[, rate := ifelse(cov > 0, mc/cov, 0)]
    .create_uid(df)  # fast numeric key
    return(df)
  })

  names(all_samples) <- names(raw_data)
  return(all_samples)
}

# Internal function as fallback using pure R
.load_data_r <- function(files.list, cores) {
  # Import helper functions using ::
  # This ensures the package works even if these packages aren't loaded
  if (cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    # Use parallel processing
    all_samples <- .load_data_parallel(files.list, cores)
  } else {
    # Use sequential processing with progress bar
    all_samples <- .load_data_sequential(files.list)
  }

  return(all_samples)
}

#' Internal function to load a single file using R
#' @return data.table of methylation data
#' @keywords internal
#' @noRd
.load_single_file_r <- function(file) {
  tryCatch({
    # Use base R to read the gzipped file
    con <- gzfile(file, "r")
    data_lines <- readLines(con)
    close(con)

    # Split the tab-delimited lines
    split_lines <- strsplit(data_lines, "\t")

    # Extract the columns we need
    n_rows <- length(split_lines)

    # Preallocate vectors for efficiency
    chr <- character(n_rows)
    pos <- integer(n_rows)
    strand <- character(n_rows)
    site <- character(n_rows)
    mc <- integer(n_rows)
    cov <- integer(n_rows)

    # Extract data column by column
    for (i in 1:n_rows) {
      line <- split_lines[[i]]
      chr[i] <- line[1]
      pos[i] <- as.integer(line[2])
      strand[i] <- line[3]
      site[i] <- line[4]
      mc[i] <- as.integer(line[5])
      cov[i] <- as.integer(line[6])
    }

    # Create the data.table
    df <- data.table::data.table(
      chr = chr,
      pos = pos,
      strand = strand,
      site = site,
      mc = mc,
      cov = cov
    )

    # Safely calculate rate
    df[, rate := ifelse(cov > 0, mc/cov, 0)]
    .create_uid(df)

    return(df)
  }, error = function(e) {
    message("Error processing file: ", file)
    message("Error message: ", e$message)
    return(NULL)
  })
}

#' Internal function to sequentially load files using R implementation
#' @return data.table of methylation data
#' @keywords internal
#' @noRd
.load_data_sequential <- function(files.list) {
  # Create progress bar
  pb <- progress::progress_bar$new(
    format = "  Loading [:bar] :percent in :elapsed",
    total = length(files.list),
    clear = FALSE
  )

  all_samples <- lapply(files.list, function(file) {
    pb$tick()
    .load_single_file_r(file)
  })

  names(all_samples) <- names(files.list)
  return(all_samples)
}

#' Internal function to sequentially load files in parallel using R implementation
#' @return data.table of methylation data
#' @keywords internal
#' @noRd
.load_data_parallel <- function(files.list, cores) {
  all_samples <- NULL

  if (cores > 1) {
    # Parallel processing
    tryCatch({
      message("Cluster created at: ", Sys.time())
      cl <- parallel::makeCluster(min(cores, parallel::detectCores() - 1))
      on.exit(parallel::stopCluster(cl), add = TRUE)

      # Export the loading function
      parallel::clusterExport(cl = cl, varlist = c("load_single_file_r"), envir = environment())

      parallel::clusterEvalQ(cl, {
        library(data.table)
        NULL
      })

      all_samples <- parallel::parLapply(cl, files.list, .load_single_file_r)
      names(all_samples) <- names(files.list)

    }, error = function(e) {
      message("Error in parallel processing: ", e$message)
      message("Falling back to sequential processing...")
      return(NULL)
    })
  }

  # If parallel processing failed, fall back to sequential
  if (is.null(all_samples)) {
    all_samples <- .load_data_sequential(files.list)
  }

  return(all_samples)
}

#' Internal function post-process loaded samples
#' @return data.table, or list of data.tables, of methylation data
#' @keywords internal
#' @noRd
.post_process_samples <- function(all_samples, sample_data, files.list) {
  # Common post-processing for both implementations
  valid_samples <- !sapply(all_samples, is.null)
  all_samples <- all_samples[valid_samples]

  if (length(all_samples) == 0) {
    warning("No files were successfully loaded!")
    attr(all_samples, "sample_metadata") <- sample_data[0,]
    return(all_samples)
  }

  valid_sample_ids <- names(files.list)[valid_samples]
  filtered_sample_data <- sample_data[match(valid_sample_ids, sample_data$sample_id)]
  attr(all_samples, "sample_metadata") <- filtered_sample_data

  return(all_samples)
}


# Coverage and Methylation Normalization Utilities ------------------------

#' Internal function to normalize a group of samples
#' @param replicates List of replicate data frames
#' @param sf Vector of scaling factors
#' @return List of normalized data frames
#' @keywords internal
.normalize_group <- function(replicates, sf) {
  purrr::map2(replicates, sf, function(df, s) {
    if(is.na(s)) stop("Scaling factor is NA")
    dt <- data.table::copy(df)
    data.table::setDT(dt)
    data.table::setDT(dt)[, c("cov", "mc", "rate") := list(cov/s,
                                                           mc/s,
                                                           mc/cov
    )]
    return(dt)
  })
}

#' Internal function to normalize methylation rates within a set of samples
#' @param sample_list List of replicate data frames
#' @param alpha Numeric between 0 and 1 controlling structure preservation (default = 0.3)
#' @param sites_per_quantile Target number of sites per quantile (default = 1000)
#' @param max_quantiles Maximum number of quantiles to use (default = 50)
#' @param diagnostics Logical indicating whether to print diagnostic information
#' @keywords internal
.normalize_methylation_within_set <- function(sample_list,
                                             alpha = 0.3,
                                             sites_per_quantile = 1000,
                                             max_quantiles = 50,
                                             diagnostics = TRUE) {

  if(diagnostics) cat("Calculating average rates...\n")

  # Calculate average methylation rates
  rate_mat <- do.call(cbind, lapply(sample_list, `[[`, "rate"))
  avg_rates <- rowMeans(rate_mat)
  
  if (!is.matrix(rate_mat) || ncol(rate_mat) != length(sample_list)) {
    stop("Samples have inconsistent row counts — check find_shared_sites was applied")
  }

  # Dynamic quantile calculation
  n_sites <- length(avg_rates)
  n_quantiles <- as.integer(max(5, min(n_sites/sites_per_quantile, max_quantiles)))
  if(diagnostics) cat(sprintf("Using %d quantiles for %d sites\n", n_quantiles, n_sites))

  # Create quantiles
  probs <- seq(0, 1, length.out = n_quantiles + 1)
  breaks <- unique(quantile(avg_rates, prob = probs, names = FALSE))

  if(length(breaks) < 3) {
    warning("Too few unique values for quantile binning. Reducing number of quantiles.")
    n_quantiles <- length(unique(avg_rates)) - 1
    if(n_quantiles < 2) stop("Not enough unique values for quantile normalization")
    probs <- seq(0, 1, length.out = n_quantiles + 1)
    breaks <- unique(quantile(avg_rates, prob = probs, names = FALSE))
  }

  rate_quantiles <- cut(avg_rates, breaks = breaks, include.lowest = TRUE)

  # Check that quantiles are valid and match data size
  if (length(rate_quantiles) != n_sites) {
    stop(sprintf("Quantile vector length (%d) doesn't match site count (%d)",
                 length(rate_quantiles), n_sites))
  }

  # Initialize result list and progress bar
  result <- sample_list
  if(diagnostics) {
    cat("Normalizing by quantile...\n")
    pb <- txtProgressBar(min = 0, max = length(levels(rate_quantiles)), style = 3)
  }

  # Normalize each quantile
  for(i in seq_along(levels(rate_quantiles))) {
    quant <- levels(rate_quantiles)[i]
    sites_in_quantile <- which(rate_quantiles == quant)

    # Check for empty quantiles
    if(length(sites_in_quantile) == 0) {
      if(diagnostics) cat(sprintf("Warning: Quantile %s is empty, skipping\n", quant))
      if(diagnostics) setTxtProgressBar(pb, i)
      next
    }

    # Validate site indices
    max_sites <- min(sapply(result, nrow))
    if(any(sites_in_quantile > max_sites)) {
      problem_sites <- sites_in_quantile[sites_in_quantile > max_sites]
      if(diagnostics) cat(sprintf("Warning: Removing %d out-of-bounds sites\n",
                                  length(problem_sites)))
      sites_in_quantile <- sites_in_quantile[sites_in_quantile <= max_sites]

      # Skip if no valid sites remain
      if(length(sites_in_quantile) == 0) {
        if(diagnostics) setTxtProgressBar(pb, i)
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

    if(diagnostics) setTxtProgressBar(pb, i)
  }

  if(diagnostics) {
    close(pb)
    cat("\nNormalization complete!\n")
  }

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
  data <- load_data(single_file = file_path)

  # Ensure required columns
  required_cols <- c("chr", "pos", "rate", "cov")
  if (!all(required_cols %in% names(data))) {
    stop("Loaded data missing required columns: ",
         paste(setdiff(required_cols, names(data)), collapse = ", "))
  }

  # Ensure correct column types
  data$chr <- as.character(data$chr)
  data$pos <- as.integer(data$pos)
  data$rate <- as.numeric(data$rate)
  data$cov <- as.integer(data$cov)

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
  # Validate input data frame
  required_cols <- c("chr", "pos", "rate", "cov")
  if (!all(required_cols %in% names(sample_df))) {
    stop("Input data frame missing required columns: ",
         paste(setdiff(required_cols, names(sample_df)), collapse = ", "))
  }

  # Ensure correct types
  sample_df$chr <- as.character(sample_df$chr)
  sample_df$pos <- as.integer(sample_df$pos)
  sample_df$rate <- as.numeric(sample_df$rate)
  sample_df$cov <- as.integer(sample_df$cov)

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
    chr = paste0("chr", sample_df$chr),
    start = sample_df$pos - 1,  # BED format is 0-based
    end = sample_df$pos,
    value = if(type == "meth") sample_df$rate else sample_df$cov
  )

  # Create output filename
  if(is.null(sample_name)) {
    file_prefix <- "sample"
    sample_name <- basename(tempfile())
  } else {
    file_prefix <- if(is.null(group_name)) sample_name else paste(group_name, sample_name, sep="_")
  }

  out_file <- file.path(out, paste0(file_prefix, "_", type, ".bedGraph"))

  data.table::fwrite(bg_data, out_file,
                     sep = "\t",
                     quote = FALSE,
                     col.names = FALSE,
                     scipen = 999
  )

  return(out_file)
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
#   - chr_int: 1-22 for autosomes, 23=X, 24=Y, 25=M/MT
#   - pos:     genomic position (up to ~250M for human, fits in 1e10 space)
#   - site_int: 1-9 for methylation context (GCH=1, HCG=2, GCG=3, etc.)
#
# Max possible value: 25e10 + 2.5e9 + 9 ≈ 2.75e10
# R double precision: exact integers up to 2^53 ≈ 9e15 → no loss of precision
# =============================================================================

#' Chromosome name to integer mapping
#' @keywords internal
#' @noRd
.chr_to_int <- function(chr) {
  # Strip "chr" prefix if present
  chr_clean <- sub("^chr", "", chr, ignore.case = TRUE)
  
  # Map to integer
  result <- integer(length(chr_clean))
  result[chr_clean == "X"] <- 23L
  result[chr_clean == "Y"] <- 24L
  result[chr_clean %in% c("M", "MT")] <- 25L
  
  # Numeric chromosomes
  numeric_idx <- !chr_clean %in% c("X", "Y", "M", "MT")
  result[numeric_idx] <- as.integer(chr_clean[numeric_idx])
  
  return(result)
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
#' @details The encoding is: \code{uid = chr_int * 1e10 + pos * 10 + site_int}, where
#'   chr_int maps chromosomes to integers (1-22, X=23, Y=24, M=25) and site_int
#'   maps methylation contexts to integers (GCH=1, HCG=2, etc.). This fits within
#'   R's double-precision integer range (exact up to 2^53) with no loss of precision.
#'
#' @keywords internal
#' @noRd
.create_uid <- function(dt, in_place = TRUE) {
  uid_vec <- .chr_to_int(dt$chr) * 1e10 + dt$pos * 10 + .site_to_int(dt$site)
  
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

