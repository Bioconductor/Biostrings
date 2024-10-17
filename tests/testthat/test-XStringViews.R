## XStringViews-class.R exports the following:
## - `==` for XStringViews against XStringViews, XString, character
## - `!=` for the same as above
## - as.character
## - as.matrix
## - toString


slen <- 25L
dnas <- sample(DNA_BASES, slen, replace=TRUE)
d <- DNAString(paste(dnas, collapse=''))

test_that('XStringViews equality works correctly', {
	n_ranges <- 5L
	starts <- sample(slen-5L, n_ranges)
	ends <- starts + sample(5L, n_ranges)
	strs <- vapply(seq_len(n_ranges), \(i) paste(dnas[seq(starts[i], ends[i])], collapse=''), character(1L))
	v <- Views(d, start=starts, end=ends)
	expect_equal(seqtype(v), seqtype(d))

	## only the diagonal should be true
	exp_opt <- matrix(FALSE, nrow=n_ranges, ncol=n_ranges)
	diag(exp_opt) <- TRUE
	test1 <- matrix(unlist(v == DNAStringSet(strs)), nrow=n_ranges, byrow=TRUE)
	expect_identical(test1, exp_opt)

	## only the diagonal should be false
	exp_opt <- !exp_opt
	test2 <- matrix(unlist(v != DNAStringSet(strs)), nrow=n_ranges, byrow=TRUE)
	expect_identical(test2, exp_opt)

	## comparison against Views objects should be supported
	expect_identical(v==v, rep(TRUE, n_ranges))

	## comparison with non-BString is not supported
	expect_error(v == strs, "is not supported")
	expect_error(strs == v, "is not supported")

	## comparing BString against character
	v2 <- v
	expect_identical((seqtype(v2) <- "B"), "B")
	expect_identical(v2 == strs, rep(TRUE, n_ranges))
	expect_identical(strs == v2, rep(TRUE, n_ranges))
	expect_identical(v2 != strs, rep(FALSE, n_ranges))
	expect_identical(strs != v2, rep(FALSE, n_ranges))

	expect_identical(toString(v2), paste(strs, collapse=', '))
	expect_identical(nchar(v2), nchar(strs))
})
