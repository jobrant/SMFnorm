# tests/testthat/test-uid-functions.R

test_that(".chr_to_int maps standard chromosomes correctly", {
  # Bare numbers
  expect_equal(.chr_to_int("1"), 1L)
  expect_equal(.chr_to_int("22"), 22L)

  # With chr prefix
  expect_equal(.chr_to_int("chr1"), 1L)
  expect_equal(.chr_to_int("chr22"), 22L)

  # Sex and mitochondrial
  expect_equal(.chr_to_int("X"), 23L)
  expect_equal(.chr_to_int("Y"), 24L)
  expect_equal(.chr_to_int("M"), 25L)
  expect_equal(.chr_to_int("MT"), 25L)
  expect_equal(.chr_to_int("chrX"), 23L)

  # Vectorized
  result <- .chr_to_int(c("1", "X", "22", "M"))
  expect_equal(result, c(1L, 23L, 22L, 25L))
})


test_that(".site_to_int maps known contexts", {
  expect_equal(.site_to_int("GCH"), 1L)
  expect_equal(.site_to_int("HCG"), 2L)
  expect_equal(.site_to_int("GCG"), 3L)

  # Unknown context gets 9
  expect_equal(.site_to_int("ZZZ"), 9L)

  # Vectorized
  result <- .site_to_int(c("GCH", "HCG", "GCG"))
  expect_equal(result, c(1L, 2L, 3L))
})


test_that(".create_uid produces correct encoding", {
  dt <- data.table(
    chr = c("1", "2", "X"),
    pos = c(100, 200, 300),
    site = c("GCH", "HCG", "GCH")
  )

  .create_uid(dt)

  expect_true("uid" %in% names(dt))
  expect_equal(dt$uid[1], 1 * 1e10 + 100 * 10 + 1)   # chr1, pos 100, GCH
  expect_equal(dt$uid[2], 2 * 1e10 + 200 * 10 + 2)   # chr2, pos 200, HCG
  expect_equal(dt$uid[3], 23 * 1e10 + 300 * 10 + 1)  # chrX, pos 300, GCH
})


test_that(".create_uid with in_place=FALSE returns vector without modifying dt", {
  dt <- data.table(
    chr = c("1", "2"),
    pos = c(100, 200),
    site = c("GCH", "GCH")
  )

  result <- .create_uid(dt, in_place = FALSE)

  expect_true(is.numeric(result))
  expect_equal(length(result), 2)
  expect_false("uid" %in% names(dt))
})


test_that(".create_uid values are unique for different sites", {
  dt <- data.table(
    chr  = c("1", "1", "1", "2"),
    pos  = c(100, 100, 200, 100),
    site = c("GCH", "HCG", "GCH", "GCH")
  )

  .create_uid(dt)

  # All should be unique
  expect_equal(length(unique(dt$uid)), 4)
})


test_that(".create_string_uid matches legacy paste approach", {
  dt <- data.table(
    chr = c("1", "2", "X"),
    pos = c(100, 200, 300),
    site = c("GCH", "HCG", "GCH")
  )

  .create_string_uid(dt)

  expect_equal(dt$uniqueID[1], "1_100_GCH")
  expect_equal(dt$uniqueID[2], "2_200_HCG")
  expect_equal(dt$uniqueID[3], "X_300_GCH")
})
