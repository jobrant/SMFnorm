# RangedSummarizedExperiment and BSseq objects in, the same class out.
#
# The object route must give exactly the numbers the list route gives: it is
# the same normalization, with a conversion on each side.

.obj_lol <- function(groups = c("A", "B"), n_sites = 600, seed = 11) {
    make_test_split_data(n_reps_per_group = 2, n_sites = n_sites,
                         group_names = groups, seed = seed)
}

.rates <- function(lol) do.call(cbind, lapply(unlist(unname(lol), recursive = FALSE),
                                              function(d) d$rate))

# Row order of an object made by create_RSE(), relative to the list's rows.
.row_order <- function(rse, dt) {
    rr <- SummarizedExperiment::rowRanges(rse)
    match(paste(as.character(GenomicRanges::seqnames(rr)), GenomicRanges::start(rr)),
          paste(dt$chr, dt$pos))
}


test_that("an RSE gives the same rates as the list route, and comes back an RSE", {
    lol <- .obj_lol()
    rse <- create_RSE(lol, verbose = FALSE)
    for (between in c(TRUE, FALSE)) {
        ref <- quiet_within_only(suppressMessages(normalize_methylation_data(
            lol, min_coverage = 0, rate_between_groups = between)))
        out <- quiet_within_only(suppressMessages(normalize_methylation_data(
            rse, min_coverage = 0, rate_between_groups = between)))

        expect_s4_class(out, "RangedSummarizedExperiment")
        expect_identical(dim(out), dim(rse))
        expect_identical(colnames(out), colnames(rse))
        expect_identical(SummarizedExperiment::rowRanges(out),
                         SummarizedExperiment::rowRanges(rse))
        ord <- .row_order(rse, lol[[1]][[1]])
        expect_equal(unname(SummarizedExperiment::assay(out, "rate")),
                     unname(.rates(ref))[ord, ], tolerance = 1e-12)
        # Coverage is never changed.
        expect_identical(SummarizedExperiment::assay(out, "coverage"),
                         SummarizedExperiment::assay(rse, "coverage"))
        # The input object is not modified.
        expect_false(isTRUE(all.equal(SummarizedExperiment::assay(out, "rate"),
                                      SummarizedExperiment::assay(rse, "rate"))))
    }
})

test_that("the settings used are recorded in the object's metadata", {
    rse <- create_RSE(.obj_lol(), exp_metadata = list(run = "x"), verbose = FALSE)
    out <- suppressMessages(normalize_methylation_data(
        rse, min_coverage = 0, rate_between_groups = TRUE, max_quantiles = 20))
    md <- S4Vectors::metadata(out)
    expect_identical(md$experiment, list(run = "x"))
    expect_true(md$SMFnorm$rate_between_groups)
    expect_equal(md$SMFnorm$max_quantiles, 20)
    expect_equal(md$SMFnorm$n_sites, nrow(out))
})

test_that("the coverage filter subsets the object's rows, in order", {
    lol <- .obj_lol()
    rse <- create_RSE(lol, verbose = FALSE)
    thr <- 18
    cov <- SummarizedExperiment::assay(rse, "coverage")
    keep <- rowSums(cov >= thr) == ncol(cov)
    expect_true(any(keep) && !all(keep))

    out <- suppressMessages(normalize_methylation_data(
        rse, min_coverage = thr, rate_between_groups = TRUE))
    expect_equal(nrow(out), sum(keep))
    expect_identical(rownames(out), rownames(rse)[keep])

    ref <- suppressMessages(normalize_methylation_data(
        lol, min_coverage = thr, rate_between_groups = TRUE))
    kept_dt <- ref[[1]][[1]]
    ord <- .row_order(out, kept_dt)
    expect_equal(unname(SummarizedExperiment::assay(out, "rate")),
                 unname(.rates(ref))[ord, ], tolerance = 1e-12)
    expect_error(suppressMessages(normalize_methylation_data(
        rse, min_coverage = 1e6, rate_between_groups = TRUE)), "No sites are left")
})

test_that("group_names and group_col select and define the groups", {
    lol <- .obj_lol(groups = c("A", "B", "C"))
    rse <- create_RSE(lol, verbose = FALSE)
    out <- suppressMessages(normalize_methylation_data(
        rse, group_names = c("C", "A"), min_coverage = 0, rate_between_groups = TRUE))
    expect_identical(as.character(SummarizedExperiment::colData(out)$group),
                     c("C", "C", "A", "A"))
    ref <- suppressMessages(normalize_methylation_data(
        lol, group_names = c("C", "A"), min_coverage = 0, rate_between_groups = TRUE))
    ord <- .row_order(out, lol[[1]][[1]])
    expect_equal(unname(SummarizedExperiment::assay(out, "rate")),
                 unname(.rates(ref))[ord, ], tolerance = 1e-12)

    # Another grouping column
    SummarizedExperiment::colData(rse)$batch <- rep(c("x", "y"), 3)
    by_batch <- suppressMessages(normalize_methylation_data(
        rse, group_col = "batch", min_coverage = 0, rate_between_groups = FALSE))
    expect_false(isTRUE(all.equal(SummarizedExperiment::assay(by_batch, "rate"),
                                  SummarizedExperiment::assay(out, "rate")[, 1:2])))

    expect_error(normalize_methylation_data(rse, group_col = "nope",
                                            rate_between_groups = TRUE),
                 "no column 'nope'.*Columns available: group, sample, batch")
    expect_error(normalize_methylation_data(rse, group_names = "Z",
                                            rate_between_groups = TRUE),
                 "not found in colData column 'group': Z")
})

test_that("the between-group choice is required for objects too", {
    rse <- create_RSE(.obj_lol(), verbose = FALSE)
    expect_error(suppressMessages(normalize_methylation_data(rse, min_coverage = 0)),
                 "'rate_between_groups' must be set")
    expect_error(suppressMessages(normalize_methylation_rates(rse, diagnostics = FALSE)),
                 "'between_groups' must be set")
})

test_that("normalize_methylation_rates accepts an object", {
    lol <- .obj_lol()
    rse <- create_RSE(lol, verbose = FALSE)
    out <- suppressMessages(normalize_methylation_rates(
        rse, between_groups = TRUE, diagnostics = FALSE))
    ref <- suppressMessages(normalize_methylation_rates(
        lol, between_groups = TRUE, diagnostics = FALSE))
    ord <- .row_order(rse, lol[[1]][[1]])
    expect_equal(unname(SummarizedExperiment::assay(out, "rate")),
                 unname(.rates(ref))[ord, ], tolerance = 1e-12)
})

test_that("sites with a missing value are dropped, with a message", {
    rse <- create_RSE(.obj_lol(), verbose = FALSE)
    r <- SummarizedExperiment::assay(rse, "rate"); r[c(3, 40), 2] <- NA
    SummarizedExperiment::assay(rse, "rate") <- r
    expect_message(out <- normalize_methylation_data(
        rse, min_coverage = 0, rate_between_groups = TRUE),
        "Dropping 2 of 600 sites")
    expect_equal(nrow(out), 598L)
    expect_false(anyNA(SummarizedExperiment::assay(out, "rate")))
})

test_that("objects without the needed assays or names are refused clearly", {
    rse <- create_RSE(.obj_lol(), verbose = FALSE)
    bad <- rse
    SummarizedExperiment::assays(bad) <- list(
        counts = SummarizedExperiment::assay(rse, "coverage"))
    expect_error(normalize_methylation_data(bad, rate_between_groups = TRUE),
                 "needs assays 'coverage' and 'rate'.*found: counts")
})

test_that("a SummarizedExperiment without ranges is accepted", {
    lol <- .obj_lol()
    rse <- create_RSE(lol, verbose = FALSE)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(coverage = SummarizedExperiment::assay(rse, "coverage"),
                      rate = SummarizedExperiment::assay(rse, "rate")),
        colData = SummarizedExperiment::colData(rse))
    out <- suppressMessages(normalize_methylation_data(
        se, min_coverage = 0, rate_between_groups = TRUE))
    ref <- suppressMessages(normalize_methylation_data(
        rse, min_coverage = 0, rate_between_groups = TRUE))
    expect_equal(SummarizedExperiment::assay(out, "rate"),
                 SummarizedExperiment::assay(ref, "rate"))
})

test_that("a BSseq object gives counts matching the list route, as a BSseq", {
    skip_if_not_installed("bsseq")
    lol <- .obj_lol()
    flat <- unlist(unname(lol), recursive = FALSE)
    first <- flat[[1]]
    M <- do.call(cbind, lapply(flat, function(d) as.integer(d$mc)))
    Cov <- do.call(cbind, lapply(flat, function(d) as.integer(d$cov)))
    o <- order(first$chr, first$pos)
    bs <- suppressWarnings(bsseq::BSseq(
        chr = first$chr[o], pos = first$pos[o], M = M[o, ], Cov = Cov[o, ],
        sampleNames = names(flat)))
    SummarizedExperiment::colData(bs)$condition <- rep(names(lol), lengths(lol))

    expect_error(normalize_methylation_data(bs, rate_between_groups = TRUE),
                 "no column 'group'")
    out <- suppressMessages(normalize_methylation_data(
        bs, group_col = "condition", min_coverage = 0, rate_between_groups = TRUE))
    ref <- suppressMessages(normalize_methylation_data(
        lol, min_coverage = 0, rate_between_groups = TRUE))

    expect_s4_class(out, "BSseq")
    expect_identical(dim(out), dim(bs))
    ref_mc <- do.call(cbind, lapply(unlist(unname(ref), recursive = FALSE),
                                    function(d) as.integer(d$mc)))
    rr <- SummarizedExperiment::rowRanges(out)
    ord <- match(paste(as.character(GenomicRanges::seqnames(rr)), GenomicRanges::start(rr)),
                 paste(first$chr, first$pos))
    got <- as.matrix(bsseq::getCoverage(out, type = "M"))
    expect_identical(unname(got), unname(ref_mc)[ord, ])
    expect_identical(as.matrix(bsseq::getCoverage(out, type = "Cov")),
                     as.matrix(bsseq::getCoverage(bs, type = "Cov")))
    expect_true(all(got <= as.matrix(bsseq::getCoverage(out, type = "Cov"))))
    expect_type(got, "integer")
})

test_that("an HDF5-backed object is read and normalized like an in-memory one", {
    skip_if_not_installed("HDF5Array")
    rse <- create_RSE(.obj_lol(), verbose = FALSE)
    h5 <- rse
    SummarizedExperiment::assays(h5, withDimnames = FALSE) <- list(
        coverage = HDF5Array::writeHDF5Array(SummarizedExperiment::assay(rse, "coverage")),
        rate = HDF5Array::writeHDF5Array(SummarizedExperiment::assay(rse, "rate")))
    expect_s4_class(SummarizedExperiment::assay(h5, "rate", withDimnames = FALSE),
                    "DelayedMatrix")
    out <- suppressMessages(normalize_methylation_data(
        h5, min_coverage = 0, rate_between_groups = TRUE))
    ref <- suppressMessages(normalize_methylation_data(
        rse, min_coverage = 0, rate_between_groups = TRUE))
    expect_equal(unname(as.matrix(SummarizedExperiment::assay(out, "rate"))),
                 unname(SummarizedExperiment::assay(ref, "rate")))
})

test_that("the plots accept objects, alone or as stages", {
    lol <- .obj_lol(n_sites = 2000)
    rse <- create_RSE(lol, verbose = FALSE)
    norm <- suppressWarnings(suppressMessages(normalize_methylation_data(
        rse, min_coverage = 0, rate_between_groups = FALSE)))
    p <- plot_efficiency(list(Before = rse, After = norm), n_bins = 10)
    d <- data.table::as.data.table(p$data)
    expect_identical(levels(d$stage), c("Before", "After"))
    expect_setequal(as.character(unique(d$group)), c("A", "B"))
    ref <- data.table::as.data.table(plot_efficiency(lol, n_bins = 10)$data)
    expect_equal(d[stage == "Before"][order(group, sample, bin), difference],
                 ref[order(group, sample, bin), difference])
    expect_no_error(ggplot2::ggplot_build(plot_density(rse)))
    expect_no_error(ggplot2::ggplot_build(plot_pca(norm)))
})
