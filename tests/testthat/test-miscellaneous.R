## Miscellaneous tests
## these are all relatively low priority and/or for files with only a few things
## some tests are just for internal functions

test_that("utils functions work as they should", {
    expect_true(Biostrings:::isNumericOrNAs(NA_character_))
    expect_true(Biostrings:::isNumericOrNAs(NA_real_))
    expect_true(Biostrings:::isNumericOrNAs(1))
    expect_true(Biostrings:::isNumericOrNAs(1L))
    expect_true(Biostrings:::isNumericOrNAs(1:2))
    expect_false(Biostrings:::isNumericOrNAs(NULL))
    expect_false(Biostrings:::isNumericOrNAs("a"))

    expect_error(Biostrings:::pow.int("a", 1), "must be a numeric vector")
    expect_identical(Biostrings:::pow.int(3,5), as.integer(3**5))

    expect_error(Biostrings:::normargUseNames(NA), "must be TRUE or FALSE")
    expect_true(Biostrings:::normargUseNames(NULL))
    expect_true(Biostrings:::normargUseNames(TRUE))
    expect_false(Biostrings:::normargUseNames(FALSE))
})

test_that("longestConsecutive still functions", {
    ## adapted from the examples in the man page
    v <- c("AAACTGTGFG", "GGGAATT", "CCAAAAAAAAAATT")
    expect_equal(longestConsecutive(v, "A"), c(3L, 2L, 10L))
    expect_equal(longestConsecutive(v, "C"), c(1L, 0L, 2L))
    expect_equal(longestConsecutive(v, "C"), c(1L, 0L, 2L))
    expect_error(longestConsecutive(v, NA),
                 "'letter' must be a character variable")
    expect_error(longestConsecutive(NA, "A"), "'x' must be a string")
})

