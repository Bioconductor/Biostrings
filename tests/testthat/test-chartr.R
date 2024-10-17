## chartr.R exports the following:
## - chartr
## - replaceAmbiguities
##
## chartr() is defined for XString, XStringSet,
## XStringViews, MaskedXString


test_that("chartr has correct behavior on Biostrings objects", {
	s1 <- BString("apple AppLE potato")
	s2 <- BString("peels AeeLE eotpto")
	expect_true(chartr("apple", "peels", s1) == s2)
	expect_true(all(chartr("apple", "peels", BStringSet(list(s1,s1))) == BStringSet(list(s2,s2))))

	v1 <- Views(s1, start=c(1,7,13), width=c(5,5,6))
	v2 <- Views(s2, start=c(1,7,13), width=c(5,5,6))
	expect_equal(as.character(chartr("apple", "peels", v1)), as.character(v2))

	## will break when we implement this to remind me to write tests
	ms1 <- s1
	ms2 <- s2
	m1 <- Mask(length(s1), start=7, width=5)
	masks(ms1) <- m1
	expect_error(chartr("apple", "peels", ms1), "Please complain!")
})

test_that('replaceAmbiguities works as expected', {
	dna1 <- DNAString(paste(DNA_ALPHABET, collapse=''))
	dna2 <- DNAString(paste(c(DNA_BASES, rep("N", length(dna1)-7L), "-+."), collapse=''))

	expect_equal(replaceAmbiguities(dna1), dna2)
	expect_equal(replaceAmbiguities(as(dna1, "RNAString")), as(dna2, "RNAString"))

	aa <- AAString(paste(AA_ALPHABET, collapse=''))
	bb <- BString(paste(LETTERS, collapse=''))
	expect_error(replaceAmbiguities(aa), "only supported for DNA and RNA")
	expect_error(replaceAmbiguities(bb), "only supported for DNA and RNA")
	expect_error(replaceAmbiguities("test"), "only supported for DNA and RNA")

	expect_true(all(replaceAmbiguities(DNAStringSet(list(dna1, dna1))) == DNAStringSet(list(dna2, dna2))))
	expect_true(all(replaceAmbiguities(Views(dna1, start=c(1,4,7), width=3)) == Views(dna2, start=c(1,4,7), width=3)))
})