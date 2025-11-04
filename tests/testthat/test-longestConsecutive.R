
test_that("longestConsecutive() works as expected", {
    ## adapted from the examples in the man page
    v <- c("AAACTGTGFG", "GGGAATT", "CCAAAAAAAAAATT")
    expect_equal(longestConsecutive(v, "A"), c(3L, 2L, 10L))
    expect_equal(longestConsecutive(v, "C"), c(1L, 0L, 2L))
    expect_equal(longestConsecutive(v, "C"), c(1L, 0L, 2L))
    expect_error2(longestConsecutive(v, NA),
                  "'letter' must be a character variable")
    expect_error2(longestConsecutive(NA, "A"), "'x' must be a string")
})

