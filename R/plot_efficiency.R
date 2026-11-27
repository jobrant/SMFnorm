#' Efficiency Curve: Each Sample Against the Average, by Methylation Rate
#'
#' Shows what SMFnorm corrects. Sites are grouped into bins of similar
#' methylation rate, and each sample's mean rate in a bin is compared with the
#' mean over the samples of its group (the default) or over all samples. A
#' sample labeled less efficiently than the others sits below the average,
#' usually by an amount that changes with the rate, so before normalization
#' the samples fan out along curves. After normalization they fall on the
#' reference line.
#'
#' @param meth_data A nested list of methylation data (group -> sample ->
#'   data.table), as returned by \code{\link{split_by_groups}} or the
#'   normalization functions, or a named list of such objects to compare
#'   stages, e.g. \code{list(Before = split_data, After = normalized_data)}.
#'   A \code{RangedSummarizedExperiment} or \code{BSseq} object is accepted
#'   wherever a nested list is; its samples are grouped by the \code{colData}
#'   column \code{group}.
#'   Within a stage, all samples must hold the same sites in the same order.
#' @param stages Optional character vector of stage names. Defaults to the
#'   names of \code{meth_data} when it holds several stages.
#' @param reference What each sample is compared with:
#'   \describe{
#'     \item{"group"}{(default) the mean of the samples in its own group. This
#'       is what within-group normalization corrects, and it says nothing
#'       about differences between groups. Groups with a single sample are
#'       left out.}
#'     \item{"all"}{the mean of all samples in all groups. This is what the
#'       between-group step corrects; use it to see differences between
#'       groups, or after normalizing with \code{rate_between_groups = TRUE}.}
#'   }
#' @param n_bins Number of rate bins, each holding about the same number of
#'   sites. Default 50, the default number SMFnorm normalizes over. Fewer bins
#'   are used when the data hold too few distinct rates.
#' @param scale What to show on the y-axis:
#'   \describe{
#'     \item{"difference"}{(default) sample mean minus the all-sample mean, so
#'       the reference is a horizontal line at 0 and small differences are
#'       visible.}
#'     \item{"ratio"}{sample mean divided by the all-sample mean: the sample's
#'       efficiency relative to the average. A flat line means the same
#'       proportional loss at every rate.}
#'     \item{"rate"}{the sample mean itself, with the diagonal as reference.}
#'   }
#' @param color_by Color the lines by \code{"group"} (default) or by
#'   \code{"sample"}, to tell individual samples apart.
#' @param points Logical. Mark each bin with a point (default \code{TRUE}).
#' @param theme A ggplot2 theme. Default is \code{\link[ggplot2]{theme_minimal}}.
#'
#' @return A ggplot object, one panel per stage, one line per sample. The plotted values are in its \code{data} element: \code{stage},
#'   \code{group}, \code{sample}, \code{bin}, \code{n_sites},
#'   \code{reference_rate}, \code{sample_rate}, \code{difference} and
#'   \code{ratio}.
#'
#' @details
#' Bins are formed separately for each stage, and with
#' \code{reference = "group"} for each group, from the mean rate of each site
#' over the samples compared, in the way the normalization forms them.
#' All sites are used; nothing is subsampled. Bin means are plain means over
#' sites, not weighted by coverage. With more bins here than the normalization
#' used, the curve also shows the variation left inside each normalization bin.
#'
#' With \code{reference = "group"}, within-group normalization brings every
#' line to the reference. With \code{reference = "all"}, data normalized within
#' groups only show the samples of a group together and the groups apart: the
#' between-group difference, which may be technical or biological and which
#' only the between-group step (\code{rate_between_groups = TRUE}) removes.
#'
#' @examples
#' \dontrun{
#' normalized_data <- normalize_methylation_data(split_data,
#'                                               rate_between_groups = FALSE)
#' plot_efficiency(list(Before = split_data, After = normalized_data))
#'
#' # Relative efficiency of each sample within its group
#' plot_efficiency(split_data, scale = "ratio")
#'
#' # Against all samples, to see differences between groups
#' plot_efficiency(split_data, reference = "all")
#' }
#'
#' @seealso \code{\link{normalize_methylation_rates}} for the correction this
#'   plot displays.
#'
#' @export
plot_efficiency <- function(meth_data,
                            stages = NULL,
                            reference = c("group", "all"),
                            n_bins = 50,
                            scale = c("difference", "ratio", "rate"),
                            color_by = c("group", "sample"),
                            points = TRUE,
                            theme = ggplot2::theme_minimal()) {

  scale <- match.arg(scale)
  reference <- match.arg(reference)
  color_by <- match.arg(color_by)
  if (!is.numeric(n_bins) || length(n_bins) != 1L || is.na(n_bins) || n_bins < 2) {
    stop("'n_bins' must be a single number of at least 2.")
  }

  stage_list <- .detect_stages(meth_data, stages)
  curve <- data.table::rbindlist(lapply(names(stage_list), function(stage_name) {
    .efficiency_table(stage_list[[stage_name]], n_bins, stage_name, reference)
  }))
  curve[, stage := factor(stage, levels = names(stage_list))]
  curve[, group := factor(group, levels = unique(group))]

  ref_name <- if (reference == "group") "group" else "all-sample"
  y_col <- switch(scale, difference = "difference", ratio = "ratio",
                  rate = "sample_rate")
  y_lab <- switch(scale,
                  difference = paste("Sample mean rate -", ref_name, "mean rate"),
                  ratio = paste("Sample mean rate /", ref_name, "mean rate"),
                  rate = "Sample mean rate")
  curve <- as.data.frame(curve)
  curve$line_id <- interaction(curve$stage, curve$group, curve$sample, drop = TRUE)

  p <- ggplot2::ggplot(curve, ggplot2::aes(
    x = .data[["reference_rate"]], y = .data[[y_col]],
    group = .data[["line_id"]], color = .data[[color_by]]
  ))
  p <- p + switch(
    scale,
    difference = ggplot2::geom_hline(yintercept = 0, color = "grey35",
                                     linetype = "dashed"),
    ratio = ggplot2::geom_hline(yintercept = 1, color = "grey35",
                                linetype = "dashed"),
    rate = ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey35",
                                linetype = "dashed")
  )
  p <- p + ggplot2::geom_line(linewidth = 0.6, alpha = 0.85, na.rm = TRUE)
  if (points) p <- p + ggplot2::geom_point(size = 0.9, alpha = 0.85, na.rm = TRUE)

  groups <- levels(curve$group)
  if (color_by == "group" && length(groups) <= length(.group_palette)) {
    p <- p + ggplot2::scale_color_manual(
      values = stats::setNames(.group_palette[seq_along(groups)], groups))
  }

  p <- p +
    ggplot2::facet_wrap(~stage) +
    ggplot2::labs(
      title = "Efficiency curve",
      subtitle = "Each line is a sample; each point a bin of sites with similar methylation rate",
      x = if (reference == "group") {
        "Mean methylation rate of the bin (samples of the group)"
      } else {
        "Mean methylation rate of the bin (all samples)"
      },
      y = y_lab,
      color = if (color_by == "group") "Group" else "Sample"
    ) +
    theme +
    ggplot2::theme(legend.position = "bottom")

  p
}


# Colors for groups: the Okabe-Ito set, distinguishable with the common forms
# of color blindness.
.group_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                    "#E69F00", "#56B4E9", "#000000", "#F0E442")


#' Mean rate of each sample in each rate bin, for one stage
#'
#' @param stage_data Nested group/sample list.
#' @param n_bins Number of bins asked for.
#' @param stage_name Label for the \code{stage} column.
#' @param reference "group" or "all".
#' @return data.table with one row per sample and bin.
#' @keywords internal
#' @noRd
.efficiency_table <- function(stage_data, n_bins, stage_name = "Data",
                              reference = "all") {
  is_group <- vapply(stage_data, function(x) is.list(x) && !is.data.frame(x),
                     logical(1))
  if (!is.list(stage_data) || is.data.frame(stage_data) || !all(is_group) ||
      is.null(names(stage_data))) {
    stop("Each stage must be a named list of groups, each a named list of ",
         "samples, as returned by split_by_groups().")
  }
  sample_groups <- rep(names(stage_data), lengths(stage_data))
  sample_names <- unlist(lapply(stage_data, names), use.names = FALSE)
  samples <- unlist(unname(stage_data), recursive = FALSE)
  if (length(samples) < 2L) {
    stop("An efficiency curve needs at least two samples (stage '",
         stage_name, "' has ", length(samples), ").")
  }
  if (length(sample_names) != length(samples)) {
    stop("Every sample in stage '", stage_name, "' needs a name.")
  }
  names(samples) <- paste(sample_groups, sample_names, sep = "/")

  if (reference == "group") {
    # Each group against its own mean, with its own bins
    n_in_group <- lengths(stage_data)
    single <- names(stage_data)[n_in_group < 2L]
    if (length(single) == length(stage_data)) {
      stop("No group in stage '", stage_name, "' has two or more samples. ",
           "Use reference = \"all\" to compare samples across groups.")
    }
    if (length(single)) {
      message("Stage '", stage_name, "': group(s) with a single sample left ",
              "out of the within-group curve: ", paste(single, collapse = ", "))
    }
    return(data.table::rbindlist(lapply(setdiff(names(stage_data), single), function(g) {
      use <- sample_groups == g
      .efficiency_curve(samples[use], sample_groups[use], sample_names[use],
                        n_bins, stage_name)
    })))
  }
  .efficiency_curve(samples, sample_groups, sample_names, n_bins, stage_name)
}


#' Bin means of a set of samples against the mean of that set
#' @keywords internal
#' @noRd
.efficiency_curve <- function(samples, sample_groups, sample_names, n_bins,
                              stage_name) {
  .check_sites_aligned(samples)

  rates <- lapply(samples, function(df) {
    if ("rate" %in% names(df)) {
      as.numeric(df[["rate"]])
    } else if (all(c("mc", "cov") %in% names(df))) {
      ifelse(df[["cov"]] > 0, df[["mc"]] / df[["cov"]], 0)
    } else {
      stop("Samples need a 'rate' column, or 'mc' and 'cov'.")
    }
  })

  # Bins of the across-sample mean rate, formed as the normalization forms them
  mean_rate <- round(Reduce(`+`, rates) / length(rates), .bin_digits)
  ok <- !is.na(mean_rate)
  if (sum(ok) < 2L) stop("No sites with a rate in stage '", stage_name, "'.")
  probs <- seq(0, 1, length.out = as.integer(n_bins) + 1L)
  breaks <- unique(round(stats::quantile(mean_rate[ok], probs = probs, names = FALSE),
                         .bin_digits))
  if (length(breaks) < 3L) {
    stop("The mean rates in stage '", stage_name, "' take too few distinct ",
         "values to form rate bins.")
  }
  bin <- findInterval(mean_rate, breaks, rightmost.closed = TRUE, left.open = TRUE)
  bin[bin == 0L] <- 1L                 # sites at the lowest boundary
  n_found <- length(breaks) - 1L
  n_sites <- tabulate(bin[ok], nbins = n_found)
  bin_mean <- function(x) {
    out <- rep(NA_real_, n_found)
    use <- ok & !is.na(x)
    s <- rowsum(x[use], bin[use])
    out[as.integer(rownames(s))] <- s[, 1] /
      tabulate(bin[use], nbins = n_found)[as.integer(rownames(s))]
    out
  }
  reference <- bin_mean(mean_rate)

  keep <- which(n_sites > 0L)
  data.table::rbindlist(lapply(seq_along(rates), function(i) {
    sample_rate <- bin_mean(rates[[i]])
    data.table::data.table(
      stage = stage_name,
      group = sample_groups[i],
      sample = sample_names[i],
      bin = keep,
      n_sites = n_sites[keep],
      reference_rate = reference[keep],
      sample_rate = sample_rate[keep],
      difference = sample_rate[keep] - reference[keep],
      ratio = ifelse(reference[keep] > 0, sample_rate[keep] / reference[keep],
                     NA_real_)
    )
  }))
}
