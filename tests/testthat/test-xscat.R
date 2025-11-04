
test_that("xscat() works as expected", {
    s1 <- "ATGATG"
    s2 <- "CCACCA"

    ## character
    expect_true(xscat(s1, s2) == BString(paste0(s1, s2)))
    expect_error2(xscat(s1, NA),
                  "arguments must be character vectors \\(with no NAs\\)")

    ## XString
    d1 <- DNAString(s1)
    d2 <- DNAString(s2)
    d3 <- DNAString(paste0(s1, s2))
    d4 <- DNAString(paste0(s2, s1))
    expect_true(xscat(d1, d2) == d3)

    ## XStringSet
    dss1 <- DNAStringSet(list(d1, d2))
    dss2 <- DNAStringSet(list(d2, d1))
    dss3 <- DNAStringSet(list(d3, d4))
    expect_true(all(xscat(dss1, dss2) == dss3))

    v1 <- Views(d3, start=c(1,7), width=3)
    v2 <- Views(d3, start=c(4,10), width=3)
    expect_true(all(xscat(v1, v2) == DNAStringSet(c(s1, s2))))
})

