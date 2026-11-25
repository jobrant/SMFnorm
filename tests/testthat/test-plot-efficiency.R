# The efficiency curve: each sample's mean rate per rate bin against the
# all-sample mean. It should show an efficiency difference before
# normalization and none after.

# Two groups; sample j has every true rate multiplied by eff[j].
.eff_data <- function(eff = c(1, 0.8, 0.7, 0.55), n_sites = 4000, seed = 1) {
    set.seed(seed)
    base <- data.table::data.table(
        chr = "1", pos = seq_len(n_sites) * 10L, strand = "+", site = "GCH")
    p <- stats::rbeta(n_sites, 1.2, 3)
    mk <- function(e) {
        dt <- data.table::copy(base)
        dt[, cov := 200L]
        dt[, mc := stats::rbinom(.N, cov, p * e)]
        dt[, rate := mc / cov]
        dt[]
    }
    list(A = list(a1 = mk(eff[1]), a2 = mk(eff[2])),
         B = list(b1 = mk(eff[3]), b2 = mk(eff[4])))
}

test_that("the curve recovers each sample's efficiency relative to the average", {
    eff <- c(1, 0.8, 0.7, 0.55)
    p <- plot_efficiency(.eff_data(eff), n_bins = 20, scale = "ratio",
                         reference = "all")
    expect_s3_class(p, "ggplot")

    d <- data.table::as.data.table(p$data)
    expect_setequal(names(d)[1:9],
                    c("stage", "group", "sample", "bin", "n_sites",
                      "reference_rate", "sample_rate", "difference", "ratio"))
    expect_equal(data.table::uniqueN(d$bin), 20L)
    # Multiplicative efficiency: a flat ratio of eff / mean(eff), away from the
    # lowest bins where counts are few.
    got <- d[bin > 5, .(ratio = mean(ratio)), by = sample]
    expect_equal(got$ratio, eff / mean(eff), tolerance = 0.03)
    expect_equal(d$difference, d$sample_rate - d$reference_rate)
    # Every site is in exactly one bin.
    expect_equal(d[sample == "a1", sum(n_sites)], 4000L)
})

test_that("after full normalization every sample lies on the reference", {
    # 40 normalization bins, so the 20 plotted bins are no finer than the
    # correction. (Plotted finer than it was normalized, the curve shows the
    # variation left inside each normalization bin.)
    raw <- .eff_data()
    norm <- suppressMessages(normalize_methylation_data(
        raw, min_coverage = 0, sites_per_quantile = 100, diagnostics = FALSE, rate_between_groups = TRUE))
    p <- plot_efficiency(list(Before = raw, After = norm), n_bins = 20,
                         reference = "all")
    d <- data.table::as.data.table(p$data)

    expect_identical(levels(d$stage), c("Before", "After"))
    expect_gt(d[stage == "Before", max(abs(difference))], 0.05)
    expect_lt(d[stage == "After", max(abs(difference))], 0.01)
    expect_no_error(ggplot2::ggplot_build(p))
})

test_that("within-group normalization leaves the groups apart", {
    raw <- .eff_data(eff = c(1, 0.9, 0.6, 0.5))
    within <- suppressWarnings(suppressMessages(normalize_methylation_data(
        raw, min_coverage = 0, rate_between_groups = FALSE, diagnostics = FALSE)))
    d <- data.table::as.data.table(
        plot_efficiency(within, n_bins = 10, reference = "all")$data)
    by_group <- d[bin > 3, .(m = mean(difference)), by = group]
    expect_gt(by_group[group == "A", m], 0.02)
    expect_lt(by_group[group == "B", m], -0.02)

    # The default compares each sample with its own group, where within-group
    # normalization leaves nothing; before it, the replicates differ.
    expect_identical(formals(plot_efficiency)$reference[[2]], "group")
    d0 <- data.table::as.data.table(plot_efficiency(raw, n_bins = 10)$data)
    dw <- data.table::as.data.table(plot_efficiency(within, n_bins = 5)$data)
    expect_gt(d0[, max(abs(difference))], 0.02)
    expect_lt(dw[, max(abs(difference))], 0.002)
})

test_that("the group reference recovers efficiency relative to the group", {
    eff <- c(1, 0.8, 0.7, 0.55)
    d <- data.table::as.data.table(
        plot_efficiency(.eff_data(eff), n_bins = 20, scale = "ratio")$data)
    got <- d[bin > 5, .(ratio = mean(ratio)), by = sample]
    expect_equal(got$ratio, c(1 / 0.9, 0.8 / 0.9, 0.7 / 0.625, 0.55 / 0.625),
                 tolerance = 0.03)
})

test_that("single-sample groups are left out of the group reference", {
    raw <- .eff_data()
    one <- list(A = raw$A, B = raw$B["b1"])
    expect_message(p <- plot_efficiency(one), "single sample left out.*B")
    expect_identical(as.character(unique(p$data$group)), "A")
    expect_error(plot_efficiency(list(A = raw$A["a1"], B = raw$B["b1"])),
                 "reference = \"all\"")
    expect_no_error(plot_efficiency(list(A = raw$A["a1"], B = raw$B["b1"]),
                                    reference = "all"))
})

test_that("scale and color_by options build", {
    raw <- .eff_data()
    for (sc in c("difference", "ratio", "rate")) {
        expect_no_error(ggplot2::ggplot_build(plot_efficiency(raw, scale = sc)))
    }
    p <- plot_efficiency(raw, color_by = "sample", points = FALSE)
    expect_identical(p$labels$colour, "Sample")
    expect_no_error(ggplot2::ggplot_build(p))
})

test_that("unusable input is an error that says why", {
    raw <- .eff_data()
    flat <- unlist(unname(raw), recursive = FALSE)
    expect_error(plot_efficiency(list(A = raw$A["a1"]), reference = "all"),
                 "at least two samples")
    expect_error(plot_efficiency(raw, n_bins = 1), "n_bins")
    bad <- raw; bad$B$b2 <- bad$B$b2[-1]
    expect_error(plot_efficiency(bad), "different numbers of sites")
    constant <- lapply(raw, function(g) lapply(g, function(d) {
        d <- data.table::copy(d); d[, rate := 0.5]; d }))
    expect_error(plot_efficiency(constant), "too few distinct")
})

test_that("few distinct rates give fewer bins rather than an error", {
    raw <- lapply(.eff_data(), function(g) lapply(g, function(d) {
        d <- data.table::copy(d); d[, rate := round(rate, 1)]; d }))
    d <- data.table::as.data.table(plot_efficiency(raw, n_bins = 50)$data)
    expect_lt(data.table::uniqueN(d$bin), 50L)
    expect_gt(data.table::uniqueN(d$bin), 3L)
})

test_that("visualize_normalization includes the curve and writes nothing unasked", {
    raw <- .eff_data()
    norm <- suppressMessages(normalize_methylation_data(
        raw, min_coverage = 0, diagnostics = FALSE, rate_between_groups = TRUE))
    wd <- tempfile("wd_"); dir.create(wd)
    old <- setwd(wd); on.exit(setwd(old))

    res <- suppressMessages(visualize_normalization(list(Before = raw, After = norm)))
    expect_named(res$plots, c("efficiency", "density"))
    expect_length(list.files(wd), 0L)
    expect_length(res$report_files, 0L)

    out <- file.path(wd, "report", "norm.pdf")
    res <- suppressMessages(visualize_normalization(
        list(Before = raw, After = norm), plots = c("efficiency", "pca"),
        output_file = out))
    expect_true(file.exists(out))
    expect_true(file.exists(file.path(wd, "report", "efficiency_report.pdf")))
    expect_error(visualize_normalization(list(Before = raw), plots = "qq"),
                 "Invalid plot type")
})
