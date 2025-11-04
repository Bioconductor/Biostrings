
test_that("toComplex() works as expected", {
    seq <- DNAString("agtcagtcagtcagtc")
    baseValues1 <- c(A=1+0i, G=0+1i, T=-1+0i, C=0-1i)
    expect_equal(toComplex(seq, baseValues1), unname(rep(baseValues1, 4)))

    ## GC content:
    baseValues2 <- c(A=0, C=1, G=1, T=0)
    expect_equal(sum(as.integer(toComplex(seq, baseValues2))), 8L)

    expect_error2(toComplex(seq, 1:4), "'baseValues' must have names")
    expect_error2(toComplex(seq, c(W=1,X=1,Y=1,Z=1)),
                  "'baseValues' names must be valid DNA letters")
})

