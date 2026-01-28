
test_that("strsplit() and unstrsplit() methods work as expected", {
    ds <- "ATGCATGC"
    d1 <- DNAStringSet(ds)
    d2 <- DNAStringSet("ATCATC")
    d_spl <- DNAStringSet(c("AT", "CAT", "C"))

    expect_true(all(strsplit(d1, "G") == d_spl))

    expect_true(unstrsplit(strsplit(d1, "G")) == d2)

    ## round trip identity
    expect_true(unstrsplit(strsplit(d1, "G"), "G") == d1)

    ## unstrsplit on a stringset should be a noop
    expect_identical(unstrsplit(d1), d1)
})

