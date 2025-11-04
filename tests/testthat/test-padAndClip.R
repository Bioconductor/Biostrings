
test_that("padAndClip, stackStrings have correct functionality", {
    ## just replicating the stuff in the examples
    ## can add more tests later if necessary, I don't see these used super often
    x <- BStringSet(c(seq1="ABCD", seq2="abcdefghijk", seq3="", seq4="XYZ"))

    p1 <- padAndClip(x, IRanges(3, 8:5),
                     Lpadding.letter=">", Rpadding.letter="<")
    expect_equal(unname(as.character(p1)),
                 c("CD<<<<", "cdefg", "<<<<", "Z<<"))

    p2 <- padAndClip(x, IRanges(1:-2, 7),
                     Lpadding.letter=">", Rpadding.letter="<")
    expect_equal(unname(as.character(p2)),
                 c("ABCD<<<", ">abcdefg", ">><<<<<<<", ">>>XYZ<<<<"))

    expect_equal(width(stackStrings(x, 2, 8)), rep(7L, 4L))

    exp_out <- c("###ABCD....", "ijk........", "#########..", "##########X")
    s1 <- stackStrings(x, -2, 8, shift=c(0, -11, 6, 7),
                       Lpadding.letter="#", Rpadding.letter=".")
    expect_equal(width(s1), rep(11, 4L))
    expect_equal(unname(as.character(s1)), exp_out)

    s2 <- stackStrings(x, -2, 8, shift=c(0, -14, 6, 7),
                       Lpadding.letter="#", Rpadding.letter=".",
                       remove.out.of.view.strings=TRUE)
    expect_equal(names(s2), paste0("seq", c(1,3,4)))
    expect_equal(unname(as.character(s2)), exp_out[c(1,3,4)])
})

