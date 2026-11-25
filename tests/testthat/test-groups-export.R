# Group handling in normalization and the RSE round trip, and bedGraph export.
#
# Four faults these pin down:
#  - create_RSE() recovered group and sample by splitting the column name at
#    its first underscore, so group "WT_treated" became group "WT";
#  - create_bedgraph() put "chr" in front of every chromosome name ("chrchr1")
#    and failed on a list of samples, because a data frame is also a list;
#  - `group_names` in the normalization functions was accepted and ignored;
#  - a leftover test function (testCpp) was compiled into the package.

.lol <- function(groups = c("WT_treated", "WT.ctrl"), n_sites = 60, chr = c("1", "2")) {
    x <- make_test_split_data(n_reps_per_group = 2, n_sites = n_sites,
                              group_names = groups, seed = 7)
    if (!identical(chr, c("1", "2"))) {
        for (g in names(x)) for (s in names(x[[g]])) {
            x[[g]][[s]][, chr := ifelse(chr == "1", ..chr[1], ..chr[2])]
        }
    }
    x
}


# --- create_RSE / convertSE_to_list ---------------------------------------

test_that("create_RSE keeps group names that contain underscores or dots", {
    lol <- .lol()
    rse <- create_RSE(lol, verbose = FALSE)
    cd <- SummarizedExperiment::colData(rse)

    expect_identical(as.character(cd$group),
                     rep(c("WT_treated", "WT.ctrl"), each = 2))
    expect_identical(levels(cd$group), c("WT_treated", "WT.ctrl"))
    expect_identical(cd$sample, unlist(lapply(lol, names), use.names = FALSE))
    expect_identical(colnames(rse), paste(cd$group, cd$sample, sep = "_"))
})

test_that("the RSE round trip returns the same groups, samples and values", {
    lol <- .lol()
    back <- convertSE_to_list(create_RSE(lol, verbose = FALSE), verbose = FALSE)

    expect_identical(names(back), names(lol))
    for (g in names(lol)) {
        expect_identical(names(back[[g]]), names(lol[[g]]))
        for (s in names(lol[[g]])) {
            a <- lol[[g]][[s]]; b <- back[[g]][[s]]
            b <- b[match(a$uid, b$uid)]
            expect_equal(b$rate, a$rate)
            expect_equal(b$cov, as.numeric(a$cov))
            expect_equal(b$mc, as.numeric(a$mc))
            expect_type(b$chr, "character")
        }
    }
})

test_that("convertSE_to_list selects a group whose name has an underscore", {
    rse <- create_RSE(.lol(), verbose = FALSE)
    one <- convertSE_to_list(rse, groups = "WT_treated", verbose = FALSE)
    expect_named(one, "WT_treated")
    expect_length(one[[1]], 2L)
    expect_error(convertSE_to_list(rse, groups = "WT", verbose = FALSE),
                 "None of the specified groups")
})

test_that("an RSE without group and sample columns still converts", {
    rse <- create_RSE(.lol(groups = c("A", "B")), verbose = FALSE)
    SummarizedExperiment::colData(rse)$group <- NULL
    SummarizedExperiment::colData(rse)$sample <- NULL
    back <- convertSE_to_list(rse, verbose = FALSE)
    expect_named(back, c("A", "B"))
    expect_named(back$A, c("A_Rep1", "A_Rep2"))
})

test_that("create_RSE rejects an unknown group and does not modify its input", {
    lol <- .lol()
    before <- data.table::copy(lol[[1]][[1]])
    expect_error(create_RSE(lol, groups = "nope", verbose = FALSE), "not found")
    sub <- create_RSE(lol, groups = "WT.ctrl", verbose = FALSE)
    expect_equal(ncol(sub), 2L)
    expect_identical(lol[[1]][[1]], before)
})


# --- group_names ----------------------------------------------------------

test_that("group_names restricts normalization to the named groups", {
    lol <- .lol(groups = c("A", "B", "C"), n_sites = 400)
    two <- suppressMessages(normalize_methylation_data(
        lol, group_names = c("C", "A"), min_coverage = 0, rate_between_groups = TRUE))

    expect_named(two, c("C", "A"))
    expect_setequal(attr(two, "sample_metadata")$group_id, c("A", "C"))

    # The same as normalizing a list that holds only those groups.
    ref <- suppressMessages(normalize_methylation_data(lol[c("C", "A")],
                                                       min_coverage = 0, rate_between_groups = TRUE))
    for (g in names(two)) for (s in names(two[[g]])) {
        expect_equal(two[[g]][[s]]$rate, ref[[g]][[s]]$rate)
    }
    # And not the same as normalizing all three groups together.
    all3 <- suppressMessages(normalize_methylation_data(lol, min_coverage = 0, rate_between_groups = TRUE))
    expect_false(isTRUE(all.equal(all3$A[[1]]$rate, two$A[[1]]$rate)))
})

test_that("group_names is honoured by normalize_methylation_rates too", {
    lol <- .lol(groups = c("A", "B", "C"), n_sites = 400)
    res <- suppressMessages(normalize_methylation_rates(
        lol, group_names = "B", between_groups = FALSE, diagnostics = FALSE))
    expect_named(res, "B")
})

test_that("unknown or malformed group_names is an error", {
    lol <- .lol(groups = c("A", "B"))
    expect_error(normalize_methylation_data(lol, group_names = c("A", "Z"), rate_between_groups = TRUE),
                 "not found in data_list: Z")
    expect_error(normalize_methylation_rates(lol, group_names = 1:2,
                                             diagnostics = FALSE, between_groups = TRUE),
                 "character vector")
    expect_error(normalize_methylation_data(lol, group_names = c("A", "A"), rate_between_groups = TRUE),
                 "distinct")
})

test_that("a flat list of samples is rejected with a pointer to split_by_groups", {
    flat <- unlist(unname(.lol(groups = c("A", "B"))), recursive = FALSE)
    expect_error(normalize_methylation_data(flat, rate_between_groups = TRUE), "split_by_groups")
})

test_that("normalization diagnostics are messages, with nothing on stdout", {
    lol <- .lol(groups = c("A", "B"), n_sites = 400)
    out <- capture.output(
        res <- suppressMessages(normalize_methylation_data(lol, min_coverage = 0,
                                                           diagnostics = TRUE, rate_between_groups = TRUE)))
    expect_length(out, 0L)
    expect_message(normalize_methylation_data(lol, min_coverage = 0,
                                              diagnostics = TRUE, rate_between_groups = TRUE),
                   "within-group")
})


# --- create_bedgraph ------------------------------------------------------

.read_bg <- function(f) data.table::fread(f, header = FALSE,
                                          col.names = c("chr", "start", "end", "value"),
                                          colClasses = list(character = 1))

test_that("chromosome names with a chr prefix are written unchanged", {
    lol <- .lol(groups = c("A", "B"), chr = c("chr1", "chrX"))
    for (cpp in c(TRUE, FALSE)) {
        out <- tempfile("bg_")
        f <- create_bedgraph(lol$A[[1]], out = out, sample_name = "s", use_cpp = cpp)
        expect_setequal(unique(.read_bg(f)$chr), c("chr1", "chrX"))
    }
})

test_that("names without the prefix get it, and MT becomes chrM", {
    lol <- .lol(groups = c("A", "B"), chr = c("1", "MT"))
    for (cpp in c(TRUE, FALSE)) {
        out <- tempfile("bg_")
        f <- create_bedgraph(lol$A[[1]], out = out, sample_name = "s", use_cpp = cpp)
        expect_setequal(unique(.read_bg(f)$chr), c("chr1", "chrM"))
    }
})

test_that("both writers give the same bedGraph, for rates and for coverage", {
    dt <- .lol(groups = c("A", "B"))$A[[1]]
    for (type in c("meth", "cov")) {
        a <- .read_bg(create_bedgraph(dt, type = type, out = tempfile("bg_"),
                                      sample_name = "s", use_cpp = TRUE))
        b <- .read_bg(create_bedgraph(dt, type = type, out = tempfile("bg_"),
                                      sample_name = "s", use_cpp = FALSE))
        expect_identical(a[, .(chr, start, end)], b[, .(chr, start, end)])
        expect_equal(a$value, b$value, tolerance = 1e-6)
        expect_identical(a$end - a$start, rep(1L, nrow(a)))
        expect_equal(a$value, if (type == "meth") dt$rate else as.numeric(dt$cov),
                     tolerance = 1e-6)
    }
})

test_that("a list of groups writes one file per sample, named group_sample_type", {
    lol <- .lol(groups = c("A", "B"))
    out <- tempfile("bg_")
    f <- create_bedgraph(lol, out = out)
    expect_length(f, 4L)
    expect_setequal(basename(f),
                    paste0(rep(c("A", "B"), each = 2), "_",
                           unlist(lapply(lol, names)), "_meth.bedGraph"))
    expect_true(all(file.exists(f)))
})

test_that("a flat list of samples works, with and without a group name", {
    lol <- .lol(groups = c("A", "B"))
    f <- create_bedgraph(lol$A, out = tempfile("bg_"), type = "cov")
    expect_setequal(basename(f), paste0(names(lol$A), "_cov.bedGraph"))
    expect_true(all(file.exists(f)))

    g <- create_bedgraph(lol$A, out = tempfile("bg_"), group_name = "grp")
    expect_setequal(basename(g), paste0("grp_", names(lol$A), "_meth.bedGraph"))
})

test_that("create_bedgraph rejects lists it cannot name or interpret", {
    lol <- .lol(groups = c("A", "B"))
    expect_error(create_bedgraph(unname(lol$A), out = tempfile("bg_")), "unique name")
    expect_error(create_bedgraph(list(a = lol$A[[1]], b = lol$B), out = tempfile("bg_")),
                 "list of samples")
    expect_error(create_bedgraph(lol, all = FALSE, out = tempfile("bg_")), "all=FALSE")
    expect_error(create_bedgraph(lol$A[[1]][, .(chr, pos)], out = tempfile("bg_")),
                 "missing required columns: rate")
})

test_that("a file path input is loaded and named after the file", {
    dt <- .lol(groups = c("A", "B"))$A[[1]]
    dir <- tempfile("in_"); dir.create(dir)
    path <- file.path(dir, "sampleX.tsv.gz")
    data.table::fwrite(dt[, .(chr, pos, strand, site, mc, cov)], path,
                       sep = "\t", col.names = FALSE)
    f <- suppressMessages(create_bedgraph(path, out = tempfile("bg_")))
    expect_identical(basename(f), "sampleX_meth.bedGraph")
    expect_equal(nrow(.read_bg(f)), nrow(dt))
})


# --- flat lists elsewhere -------------------------------------------------

test_that("summarize_meth handles the flat list returned by load_data", {
    lol <- .lol(groups = c("A", "B"))
    flat <- unlist(unname(lol), recursive = FALSE)
    out <- capture.output(res <- summarize_meth(flat))
    expect_true(any(grepl("flat sample list", out)))
    out <- capture.output(res <- summarize_meth(lol))
    expect_true(any(grepl("nested", out)))
})

test_that("create_RSE rejects a flat list of samples", {
    flat <- unlist(unname(.lol(groups = c("A", "B"))), recursive = FALSE)
    expect_error(create_RSE(flat, verbose = FALSE), "split_by_groups")
})


# --- leftovers ------------------------------------------------------------

test_that("the leftover C++ test function is gone", {
    expect_false(exists("testCpp", envir = asNamespace("SMFnorm"), inherits = FALSE))
    expect_false(exists("extract_group_names", envir = asNamespace("SMFnorm"),
                        inherits = FALSE))
})


# --- the between-group step has no default --------------------------------

test_that("the between-group choice must be made when there are several groups", {
    lol <- .lol(groups = c("A", "B"), n_sites = 400)
    expect_error(normalize_methylation_data(lol, min_coverage = 0),
                 "'rate_between_groups' must be set to TRUE or FALSE")
    expect_error(normalize_methylation_rates(lol, diagnostics = FALSE),
                 "'between_groups' must be set to TRUE or FALSE")
    expect_error(normalize_methylation_data(lol, rate_between_groups = NA),
                 "must be TRUE or FALSE")
    expect_error(normalize_methylation_data(lol, rate_between_groups = "yes"),
                 "must be TRUE or FALSE")
})

test_that("the choice is not asked for when it cannot matter", {
    lol <- .lol(groups = c("A", "B"), n_sites = 400)
    # One group
    one <- suppressMessages(normalize_methylation_data(lol, group_names = "A",
                                                       min_coverage = 0))
    expect_named(one, "A")
    expect_no_error(suppressMessages(
        normalize_methylation_rates(lol["B"], diagnostics = FALSE)))
    # No rate normalization at all
    expect_no_error(normalize_methylation_data(lol, normalize_rates = FALSE,
                                               min_coverage = 0))
})

test_that("TRUE and FALSE give the two different normalizations", {
    lol <- .lol(groups = c("A", "B"), n_sites = 400)
    both <- suppressMessages(normalize_methylation_data(
        lol, min_coverage = 0, rate_between_groups = TRUE))
    within <- suppressWarnings(suppressMessages(normalize_methylation_data(
        lol, min_coverage = 0, rate_between_groups = FALSE)))
    expect_false(isTRUE(all.equal(both$A[[1]]$rate, within$A[[1]]$rate)))
})
