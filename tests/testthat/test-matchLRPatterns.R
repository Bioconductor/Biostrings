## matchLRPatterns.R exports matchLRPatterns for the following:
## - XString
## - XStringViews
## - MaskedXString
##
## note that Lfixed, Rfixed are only defined for DNAString, RNAString

test_that("matchLRPatterns works for all inputs", {
	d <- DNAString("AAATTAACCCTT")

	## different values of max.mismatch
	expect_equal(matchLRPatterns("AA", "TT", 0, d), Views(d, start=2, end=5))
	expect_equal(matchLRPatterns("AA", "TT", 1, d), Views(d, start=c(1,2), end=c(5,5)))
	expect_equal(matchLRPatterns("AA", "TT", 3, d), Views(d, start=c(1,2,6), end=c(5,5,12)))
	expect_equal(matchLRPatterns("AA", "TT", 7, d), Views(d, start=c(1,2,2,6), end=c(5,5,12,12)))

	## L,Rmismatch
	expect_equal(matchLRPatterns("AA", "TT", 0, subject=d, max.Lmismatch=2), Views(d, start=c(2,9), end=c(5,12)))
	expect_equal(matchLRPatterns("AA", "TT", 0, subject=d, max.Rmismatch=2), Views(d, start=c(1,2,6), end=c(4,5,9)))

	## L,Rindels
	d2 <- DNAString("ATATTA")
	expect_equal(matchLRPatterns("AA", "TT", 0, d2, max.Lmismatch=1, with.Lindels=FALSE), Views(d2, start=2, end=5))
	expect_equal(matchLRPatterns("AA", "TT", 0, d2, max.Lmismatch=1, with.Lindels=TRUE), Views(d2, start=3, end=5))
	d3 <- DNAString("AATAT")
	expect_equal(matchLRPatterns("AA", "TT", 0, subject=d3, max.Rmismatch=1, with.Rindels=FALSE), Views(d3,start=1,end=4))
	expect_equal(matchLRPatterns("AA", "TT", 0, subject=d3, max.Rmismatch=1, with.Rindels=TRUE), Views(d3,start=1,end=3))

	## Lfixed, Rfixed
	expect_equal(matchLRPatterns("NN", "TT", 0, subject=d), Views(d))
	expect_equal(matchLRPatterns("MM", "TT", 0, subject=d, Lfixed="subject"), Views(d, start=c(2,9), end=c(5,12)))
	expect_equal(matchLRPatterns("AA", "WY", 0, subject=d, Rfixed="subject"), Views(d, start=c(1,2), end=c(4,5)))

	## with a mask
	dm <- d
	masks(dm) <- Mask(length(dm), start=8, end=10)
	expect_equal(matchLRPatterns("AA", "TT", 0, dm), Views(unmasked(dm), start=2, end=5))
	expect_equal(matchLRPatterns("AA", "TT", 1, dm), Views(unmasked(dm), start=c(1,2), end=c(5,5)))
	expect_equal(matchLRPatterns("AA", "TT", 3, dm), Views(unmasked(dm), start=c(1,2), end=c(5,5)))

	## on a Views object
	v <- Views(d, start=c(1,6), end=c(5,12))
	expect_equal(matchLRPatterns("AA", "TT", 0, v), Views(d, start=2, end=5))
	expect_equal(matchLRPatterns("AA", "TT", 1, v), Views(d, start=c(1,2), end=c(5,5)))
	expect_equal(matchLRPatterns("AA", "TT", 3, v), Views(d, start=c(1,2,6), end=c(5,5,12)))
	expect_equal(matchLRPatterns("AA", "TT", 7, v), Views(d, start=c(1,2,6), end=c(5,5,12)))
})