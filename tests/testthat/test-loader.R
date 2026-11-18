# Loading files: both readers, damaged input, chromosome names, shared sites.
#
# The C++ and R readers used to disagree (a header line was skipped by one and
# kept as a row of NAs by the other), a damaged or unreadable file loaded as a
# smaller or missing sample without an error, chromosome names other than
# 1-22/X/Y/M got no site identifier, and find_shared_sites() reordered the
# caller's tables.

.records <- function(n = 40, chrs = c("chr1", "chr2", "chrX"), seed = 1) {
    set.seed(seed)
    dt <- data.table::data.table(
        chr = rep(chrs, length.out = n),
        pos = sort(sample(1e6, n)),
        strand = sample(c("+", "-"), n, replace = TRUE),
        site = "GCH",
        cov = sample(1:40, n, replace = TRUE))
    dt[, mc := stats::rbinom(n, cov, 0.3)]
    dt[, .(chr, pos, strand, site, mc, cov)]
}

.write_allc <- function(dt, path = tempfile(fileext = ".tsv.gz"), header = FALSE,
                        extra_col = TRUE) {
    lines <- do.call(paste, c(as.list(dt), sep = "\t"))
    if (extra_col) lines <- paste(lines, 1, sep = "\t")
    if (header) lines <- c("chr\tpos\tstrand\tsite\tmc\tcov\tflag", lines)
    con <- gzfile(path, "w")
    writeLines(lines, con)
    close(con)
    path
}

.both <- function(path, ...) list(
    cpp = suppressMessages(load_single_file(path, use_cpp = TRUE, ...)),
    r   = suppressMessages(load_single_file(path, use_cpp = FALSE, ...)))

.sheet <- function(paths, groups = rep("A", length(paths))) data.table::data.table(
    group_id = groups, replicate = seq_along(paths),
    sample_id = paste0("s", seq_along(paths)), file_name = basename(paths))


# --- The two readers -------------------------------------------------------

test_that("the C++ and R readers return the same table", {
    rec <- .records()
    out <- .both(.write_allc(rec))

    expect_equal(as.data.frame(out$cpp), as.data.frame(out$r))
    expect_equal(nrow(out$cpp), nrow(rec))
    expect_equal(out$cpp$pos, rec$pos)
    expect_equal(out$cpp$mc, rec$mc)
    expect_equal(out$cpp$rate, rec$mc / rec$cov)
    expect_false(anyNA(out$cpp$uid))
})


test_that("a header line is skipped by both readers", {
    rec <- .records()
    out <- .both(.write_allc(rec, header = TRUE))

    expect_equal(nrow(out$cpp), nrow(rec))
    expect_equal(as.data.frame(out$cpp), as.data.frame(out$r))
    expect_false(anyNA(out$r$pos))
})


test_that("files with exactly six columns and plain-text files load", {
    rec <- .records()
    six <- .both(.write_allc(rec, extra_col = FALSE))
    expect_equal(nrow(six$cpp), nrow(rec))
    expect_equal(as.data.frame(six$cpp), as.data.frame(six$r))

    plain <- tempfile(fileext = ".tsv")
    writeLines(do.call(paste, c(as.list(rec), sep = "\t")), plain)
    txt <- .both(plain)
    expect_equal(nrow(txt$cpp), nrow(rec))
    expect_equal(as.data.frame(txt$cpp), as.data.frame(txt$r))
})


test_that("min_coverage filters the same rows in both readers", {
    rec <- .records(n = 200)
    out <- .both(.write_allc(rec), min_coverage = 20L)

    expect_equal(nrow(out$cpp), sum(rec$cov >= 20))
    expect_true(all(out$cpp$cov >= 20))
    expect_equal(as.data.frame(out$cpp), as.data.frame(out$r))
})


# --- Damaged input ---------------------------------------------------------

test_that("a malformed line is an error naming the file and line", {
    rec <- .records()
    path <- .write_allc(rec)
    lines <- readLines(gzfile(path))
    lines[7] <- "chr1\tnot_a_position\t+\tGCH\t1\t5"
    con <- gzfile(path, "w"); writeLines(lines, con); close(con)

    for (cpp in c(TRUE, FALSE)) {
        expect_error(suppressMessages(load_single_file(path, use_cpp = cpp)),
                     "line 7 is not a valid record")
    }
})


test_that("a line with too few columns is an error", {
    path <- .write_allc(.records())
    lines <- readLines(gzfile(path))
    lines[3] <- "chr1\t100\t+"
    con <- gzfile(path, "w"); writeLines(lines, con); close(con)

    for (cpp in c(TRUE, FALSE)) {
        expect_error(suppressMessages(load_single_file(path, use_cpp = cpp)),
                     "line 3 is not a valid record")
    }
})


test_that("more methylated reads than coverage is an error", {
    path <- .write_allc(.records())
    lines <- readLines(gzfile(path))
    lines[5] <- "chr1\t100\t+\tGCH\t9\t4\t1"
    con <- gzfile(path, "w"); writeLines(lines, con); close(con)

    for (cpp in c(TRUE, FALSE)) {
        expect_error(suppressMessages(load_single_file(path, use_cpp = cpp)),
                     "line 5 has impossible counts")
    }
})


test_that("an empty file and a header-only file are errors", {
    empty <- tempfile(fileext = ".tsv.gz")
    con <- gzfile(empty, "w"); close(con)
    header_only <- tempfile(fileext = ".tsv.gz")
    con <- gzfile(header_only, "w")
    writeLines("chr\tpos\tstrand\tsite\tmc\tcov", con); close(con)

    for (path in c(empty, header_only)) for (cpp in c(TRUE, FALSE)) {
        expect_error(suppressMessages(load_single_file(path, use_cpp = cpp)),
                     "No records found")
    }
})


test_that("a truncated gzip file is an error, not a shorter sample", {
    path <- .write_allc(.records(n = 5000, seed = 3))
    bytes <- readBin(path, "raw", file.size(path))
    cut <- tempfile(fileext = ".tsv.gz")
    writeBin(bytes[seq_len(length(bytes) %/% 2)], cut)

    for (cpp in c(TRUE, FALSE)) {
        expect_error(suppressMessages(suppressWarnings(
            load_single_file(cut, use_cpp = cpp))))
    }
})


test_that("load_data stops when a file in the sample sheet is missing or damaged", {
    dir <- tempfile("allc_"); dir.create(dir)
    p1 <- .write_allc(.records(seed = 1), file.path(dir, "a.tsv.gz"))
    p2 <- .write_allc(.records(seed = 2), file.path(dir, "b.tsv.gz"))
    sheet <- .sheet(c(p1, p2, file.path(dir, "missing.tsv.gz")))

    expect_error(suppressMessages(load_data(dir, sheet)), "missing.tsv.gz")

    # A damaged third file must not leave a two-sample result.
    bad <- file.path(dir, "missing.tsv.gz")
    con <- gzfile(bad, "w"); writeLines(c("chr1\t10\t+\tGCH\t1\t5", "garbage"), con)
    close(con)
    for (cpp in c(TRUE, FALSE)) {
        expect_error(suppressMessages(load_data(dir, sheet, use_cpp = cpp)),
                     "line 2 is not a valid record")
    }
})


test_that("load_data returns the same samples from both readers", {
    dir <- tempfile("allc_"); dir.create(dir)
    paths <- vapply(1:3, function(i) .write_allc(
        .records(seed = i), file.path(dir, paste0("s", i, ".tsv.gz"))), character(1))
    sheet <- .sheet(paths, c("A", "A", "B"))

    a <- suppressMessages(load_data(dir, sheet, use_cpp = TRUE))
    b <- suppressMessages(load_data(dir, sheet, use_cpp = FALSE))

    expect_identical(names(a), c("s1", "s2", "s3"))
    expect_equal(lapply(a, as.data.frame), lapply(b, as.data.frame))
    expect_identical(attr(a, "sample_metadata")$sample_id, c("s1", "s2", "s3"))
    expect_null(attr(a, "chr_levels"))
})


test_that("the parallel R loader runs on workers and matches the sequential one", {
    skip_on_cran()
    skip_if_not(nzchar(system.file(package = "SMFnorm")), "SMFnorm not installed")
    skip_if(parallel::detectCores() < 3, "needs at least 3 cores")
    dir <- tempfile("allc_"); dir.create(dir)
    paths <- vapply(1:3, function(i) .write_allc(
        .records(seed = i), file.path(dir, paste0("s", i, ".tsv.gz"))), character(1))
    files <- stats::setNames(paths, paste0("s", 1:3))

    msgs <- testthat::capture_messages(
        par <- .load_data_parallel(files, cores = 2))
    seq <- suppressMessages(.load_data_sequential(files))

    expect_true(any(grepl("on 2 workers", msgs)))
    expect_false(any(grepl("sequentially", msgs)))
    expect_equal(lapply(par, as.data.frame), lapply(seq, as.data.frame))
})


# --- Chromosome names ------------------------------------------------------

test_that("standard chromosome names keep their fixed codes", {
    expect_equal(.chr_to_int(c("chr1", "22", "chrX", "Y", "chrM", "MT")),
                 c(1L, 22L, 23L, 24L, 25L, 25L))
    expect_null(.chr_levels(data.table::data.table(chr = c("chr1", "X", "chrM"))))
})


test_that("any other chromosome name gets an identifier, with no NA and no warning", {
    chrs <- c("chrUn_GL000220v1", "chr1_KI270706v1_random", "chrI", "chr2L",
              "chr23", "scaffold_12", "chr1", "chrX")
    dt <- data.table::data.table(chr = chrs, pos = 100L, site = "GCH")

    expect_no_warning(uid <- .create_uid(dt, in_place = FALSE))
    expect_false(anyNA(uid))
    expect_equal(length(unique(uid)), length(chrs))
    # Standard names first in their usual order, then the rest alphabetically.
    expect_identical(.chr_levels(dt),
                     c("chr1", "chrX", "chr1_KI270706v1_random", "chr23", "chr2L",
                       "chrI", "chrUn_GL000220v1", "scaffold_12"))
})


test_that("chromosome 23 and X do not share an identifier", {
    dt <- data.table::data.table(chr = c("chr23", "chrX"), pos = 100L, site = "GCH")
    uid <- .create_uid(dt, in_place = FALSE)
    expect_false(uid[1] == uid[2])
})


test_that("samples with different sets of contigs are matched correctly", {
    dir <- tempfile("allc_"); dir.create(dir)
    common <- data.table::data.table(
        chr = c("chrI", "chrI", "chr2L", "chrUn_1", "chr1"),
        pos = c(10L, 20L, 30L, 40L, 50L), strand = "+", site = "GCH",
        mc = 1L, cov = 10L)
    only1 <- data.table::data.table(chr = c("aaa_contig", "chrI"), pos = c(5L, 99L),
                                    strand = "+", site = "GCH", mc = 2L, cov = 10L)
    only2 <- data.table::data.table(chr = c("zzz_contig", "chr2L"), pos = c(5L, 77L),
                                    strand = "+", site = "GCH", mc = 3L, cov = 10L)
    p1 <- .write_allc(rbind(only1, common), file.path(dir, "a.tsv.gz"))
    p2 <- .write_allc(rbind(common[5:1], only2), file.path(dir, "b.tsv.gz"))

    for (cpp in c(TRUE, FALSE)) {
        expect_no_warning(
            x <- suppressMessages(load_data(dir, .sheet(c(p1, p2)), use_cpp = cpp)))
        expect_false(anyNA(x$s1$uid))
        expect_false(is.null(attr(x, "chr_levels")))

        shared <- find_shared_sites(x, quiet = TRUE)
        expect_equal(nrow(shared$s1), 5L)
        expect_identical(shared$s1[, .(chr, pos)], shared$s2[, .(chr, pos)])
        expect_setequal(paste(shared$s1$chr, shared$s1$pos),
                        paste(common$chr, common$pos))
        expect_identical(shared$s1$uid, shared$s2$uid)
    }
})


test_that("separately loaded samples with different contigs are still matched", {
    a <- data.table::data.table(chr = c("ctgA", "ctgB", "ctgC"), pos = c(1L, 2L, 3L),
                                strand = "+", site = "GCH", mc = 1L, cov = 5L)
    b <- data.table::data.table(chr = c("ctgB", "ctgC", "ctgD"), pos = c(2L, 3L, 4L),
                                strand = "+", site = "GCH", mc = 1L, cov = 5L)
    x <- list(s1 = suppressMessages(load_single_file(.write_allc(a))),
              s2 = suppressMessages(load_single_file(.write_allc(b))))
    # Numbered within each file, so the stored identifiers disagree...
    expect_false(x$s1$uid[x$s1$chr == "ctgB"] == x$s2$uid[x$s2$chr == "ctgB"])

    # ...and find_shared_sites() must not rely on them.
    shared <- find_shared_sites(x, quiet = TRUE)
    expect_identical(shared$s1$chr, c("ctgB", "ctgC"))
    expect_identical(shared$s2$chr, c("ctgB", "ctgC"))
    expect_identical(shared$s1$uid, shared$s2$uid)
})


# --- find_shared_sites -----------------------------------------------------

test_that("find_shared_sites leaves the caller's tables untouched", {
    samples <- make_test_dataset(n_samples = 3, n_sites = 100, shared_frac = 0.8,
                                 group_ids = c("A", "A", "B"), seed = 7)
    # Put each sample in its own row order.
    set.seed(1)
    for (i in seq_along(samples)) samples[[i]] <- samples[[i]][sample(.N)]
    before <- lapply(samples, data.table::copy)

    out <- find_shared_sites(samples, quiet = TRUE)

    for (i in seq_along(samples)) {
        expect_identical(samples[[i]], before[[i]])
        expect_null(data.table::key(samples[[i]]))
    }
    # The result is aligned and sorted whatever the input order.
    expect_identical(out[[1]]$uid, out[[2]]$uid)
    expect_identical(out[[1]]$uid, out[[3]]$uid)
    expect_false(is.unsorted(out[[1]]$uid))
    expect_identical(out[[1]][, .(chr, pos)], out[[3]][, .(chr, pos)])
})


test_that("each sample's values stay with their own sites", {
    samples <- make_test_dataset(n_samples = 2, n_sites = 60, shared_frac = 0.5,
                                 group_ids = c("A", "B"), seed = 9)
    out <- find_shared_sites(samples, quiet = TRUE)

    for (i in 1:2) {
        orig <- samples[[i]][match(out[[i]]$uid, uid)]
        expect_identical(out[[i]]$mc, orig$mc)
        expect_identical(out[[i]]$cov, orig$cov)
    }
})


test_that("a sample containing the same site twice is an error", {
    samples <- make_test_dataset(n_samples = 2, n_sites = 30, shared_frac = 1,
                                 group_ids = c("A", "B"), seed = 3)
    samples[[2]] <- rbind(samples[[2]], samples[[2]][4])

    expect_error(find_shared_sites(samples, quiet = TRUE),
                 "contains the same site more than once")
})


test_that("tables with only a legacy uniqueID column are handled", {
    a <- data.table::data.table(uniqueID = c("1_10_GCH", "1_20_GCH", "2_5_GCH"), cov = 1:3)
    b <- data.table::data.table(uniqueID = c("2_5_GCH", "1_20_GCH", "3_1_GCH"), cov = 4:6)

    ids <- find_shared_sites(list(a = a, b = b), filter = FALSE, quiet = TRUE)
    expect_setequal(ids, c("1_20_GCH", "2_5_GCH"))

    out <- find_shared_sites(list(a = a, b = b), quiet = TRUE)
    expect_identical(out$a$uniqueID, out$b$uniqueID)
    expect_identical(out$a$cov, c(2L, 3L))
    expect_identical(out$b$cov, c(5L, 4L))
})


test_that("no shared sites gives a warning", {
    a <- data.table::data.table(chr = "chr1", pos = 1:3, site = "GCH")
    b <- data.table::data.table(chr = "chr2", pos = 1:3, site = "GCH")
    expect_warning(find_shared_sites(list(a = a, b = b), quiet = TRUE),
                   "No sites are shared")
})
