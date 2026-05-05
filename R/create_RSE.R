#' Create a RangedSummarizedExperiment from Methylation list of list of dataframes containing methylation data
#'
#' This function takes a nested named list of methylation data frames and converts them
#' into a single \code{RangedSummarizedExperiment} object. It expects all data frames to
#' contain filtered and common sites across all samples.
#'
#' @param lol A **doubly-nested named list**.
#' \itemize{
#'   \item \bold{Tier 1 (Names):} Group identifiers (e.g., "Group1", "Group2").
#'   \item \bold{Tier 2 (Names):} Sample identifiers (e.g., "Rep1", "Rep2").
#'   \item \bold{Elements:} Data frames or data.tables containing at least
#'   \code{chr, pos, strand, site, uniqueID, cov}, and \code{rate}.
#' }
#' @param groups Groups to include in the summarized experiment object. If not specified, all groups are retained.
#' @param exp_metadata Any additional list containing experiment-level metadata to be stored in the
#' \code{metadata} slot of the SummarizedExperiment.
#' @param verbose Logical. If TRUE (default), prints progress messages during processing.
#'
#' @return A \code{RangedSummarizedExperiment} object with two assays: \code{coverage}
#' and \code{rate}. Site-level genomic metadata is preserved in the
#' \code{rowRanges} slot as a \code{GRanges} object.
#'
#' @export
create_RSE <- function(lol, groups = NULL, exp_metadata = list(), verbose = TRUE) {
  
  if (!is.list(lol) || is.null(names(lol))) {
    stop("Input 'lol' must be a named list.")
  }
  if (!all(sapply(lol, function(x) is.list(x) && !is.null(names(x))))) {
    stop("All elements in 'lol' must be named lists.")
  }
  
  if(!is.null(groups)){
    if (verbose) message("Subsetting groups...")
    lol <- lol[groups]
  }
  
  # ── PART 1: Metadata Extraction ──
  key_cols <- c("chr", "pos", "strand", "site")
  
  if (verbose) message("Extracting master metadata and rowRanges...")
  first_s <- names(lol)[1]
  first_r <- names(lol[[first_s]])[1]
  
  required_cols <- c("chr", "pos", "strand", "site", "cov", "rate")
  first_df <- lol[[first_s]][[first_r]]
  if (!all(required_cols %in% names(first_df))) {
    stop("Missing required columns: ", paste(setdiff(required_cols, names(first_df)), collapse = ", "))
  }
  
  meta <- as.data.table(first_df)[, ..key_cols]
  
  # Create uid if not present, for fast alignment
  if (!"uid" %in% names(meta)) {
    .create_uid(meta)
  }
  setkey(meta, uid)
  
  # Generate string uniqueID for GRanges metadata and matrix rownames
  .create_string_uid(meta)
  
  # Construct RowRanges
  row_ranges <- GenomicRanges::GRanges(
    seqnames = as.character(meta$chr),
    ranges   = IRanges::IRanges(start = meta$pos, width = 1),
    strand   = meta$strand,
    site     = meta$site,
    uniqueID = meta$uniqueID
  )
  
  # ── PART 2: Flattening and Matrix Pre-allocation ──
  if (verbose) message("Flattening structure and pre-allocating matrices...")
  
  flat_reps <- unlist(lol, recursive = FALSE)
  rep_names <- names(flat_reps)
  
  num_rows <- nrow(meta)
  num_reps <- length(flat_reps)
  
  cov_matrix  <- matrix(NA_real_, nrow = num_rows, ncol = num_reps)
  rate_matrix <- matrix(NA_real_, nrow = num_rows, ncol = num_reps)
  
  # ── PART 3: Alignment Loop (using numeric uid) ──
  if (verbose) message("Aligning samples to the column metadata...")
  
  iteration_names <- names(flat_reps)
  
  for (i in seq_along(iteration_names)) {
    curr_rep_name <- iteration_names[i]
    if (verbose) message("   -> Processing: ", curr_rep_name)
    
    dt <- as.data.table(flat_reps[[curr_rep_name]])
    
    # Ensure uid exists
    if (!"uid" %in% names(dt)) {
      .create_uid(dt)
    }
    setkey(dt, uid)
    
    # Fast numeric key join
    aligned_dt <- dt[meta[, .(uid)]]
    
    cov_matrix[, i]  <- as.numeric(aligned_dt$cov)
    rate_matrix[, i] <- as.numeric(aligned_dt$rate)
    
    flat_reps[[curr_rep_name]] <- NULL
  }
  
  # ── PART 4: Column Metadata ──
  if (verbose) message("Building colData...")
  
  clean_names <- gsub("\\.", "_", rep_names)
  colnames(cov_matrix) <- colnames(rate_matrix) <- clean_names
  rownames(cov_matrix) <- rownames(rate_matrix) <- meta$uniqueID  # string for rownames
  
  name_split <- stringr::str_split(clean_names, pattern = "_", n = 2, simplify = TRUE)
  
  col_metadata <- data.frame(
    group      = as.factor(name_split[, 1]),
    sample     = name_split[, 2],
    row.names  = clean_names,
    stringsAsFactors = FALSE
  )
  
  # ── PART 5: Assembly ──
  if (verbose) message("Assembling RangedSummarizedExperiment...")
  
  rse <- SummarizedExperiment::SummarizedExperiment(
    assays    = S4Vectors::SimpleList(
      coverage = cov_matrix,
      rate     = rate_matrix
    ),
    rowRanges = row_ranges,
    colData   = col_metadata,
    metadata  = list(experiment = exp_metadata)
  )
  
  return(rse)
}

