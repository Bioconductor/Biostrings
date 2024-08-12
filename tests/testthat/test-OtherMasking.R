## Additional tests for the following files:
## - maskMotif.R (contains `maskMotif`, `mask`)
## - injectHardMask.R (contains `injectHardMask`)

test_that("maskMotif works for XStrings", {
	s1 <- BString("aabbaacc")
	s2 <- maskMotif(s1, "bb")
	expect_equal(as.character(s2), "aa##aacc")
	s3 <- maskMotif(s2, "baac")
	expect_equal(as.character(s3), "aa##aacc")

	## order of operations matters, masking should skip masked positions
	s4 <- maskMotif(s1, "baac")
	expect_equal(as.character(s4), "aab####c")
	expect_equal(as.character(gaps(s4)), "###baac#")
	expect_equal(as.character(as(s4, "Views")), c("aab", "c"))

	## mask DNAString by multiple motifs
	## example from ?maskMotif.R, just using DNAString objects
	x <- DNAString("ACACAACTAGATAGNACTNNGAGAGACGC")
	x1 <- maskMotif(x, DNAString("N"))
	x2 <- maskMotif(x1, DNAString("AC"))
	x3 <- maskMotif(x2, DNAString("GA"), min.block.width=5)
	expect_equal(width(as(x3, "Views")), c(1,7,1,5,2))
	exp_opt <- c(A=5L, C=5L, N=3L, `#`=16L)
	expect_true(all(table(strsplit(as.character(gaps(x3)), '')[[1]]) == as.table(exp_opt[order(names(exp_opt))])))
})

test_that("injectHardMask works for all supported inputs", {
	set.seed(50L)
	ds <- sample(DNA_BASES, 50L, replace=TRUE)

	d1 <- DNAString(paste(ds, collapse=''))
	d2 <- d1

	mask_starts <- c(4,10,28)
	mask_ends <- c(7,17,38)
	m1 <- Mask(length(ds), start=mask_starts, end=mask_ends)
	masks(d2) <- m1

	all_masked_indices <- unlist(lapply(seq_along(mask_starts), \(i) seq(mask_starts[i], mask_ends[i])))
	exp_opt <- ds
	exp_opt[all_masked_indices] <- "+"
	expect_equal(as.character(injectHardMask(d2)), paste(exp_opt, collapse=''))

	exp_opt[all_masked_indices] <- "M"
	expect_equal(as.character(injectHardMask(d2, "M")), paste(exp_opt, collapse=''))

	v <- Views(d1, start=mask_starts, end=mask_ends)
	## injecting on a Views does the opposite
	all_vmasked <- seq_len(length(ds))[-all_masked_indices]
	exp_vopt <- ds
	exp_vopt[all_vmasked] <- "+"
	expect_equal(as.character(injectHardMask(v)), paste(exp_vopt, collapse=''))

	exp_vopt[all_vmasked] <- "M"
	expect_equal(as.character(injectHardMask(v, "M")), paste(exp_vopt, collapse=''))
})