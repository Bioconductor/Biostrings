## Miscellaneous tests
## these are all relatively low priority and/or for files with only a few things
## some tests are just for internal functions

test_that("coloring works for DNA, RNA, and AA", {
	## not a super important test
	dna_rna_expected <- c(DNA_BASES, "U", DNA_ALPHABET[-c(1:4,16:18)])
	expect_true(!any(duplicated(Biostrings:::make_DNA_AND_RNA_COLORED_LETTERS())))
	expect_equal(sort(names(Biostrings:::make_DNA_AND_RNA_COLORED_LETTERS())), sort(dna_rna_expected))

	aa_expected <- AA_ALPHABET[-c(27:30)]
	expect_true(!any(duplicated(Biostrings:::make_AA_COLORED_LETTERS())))
	expect_equal(sort(names(Biostrings:::make_AA_COLORED_LETTERS())), sort(aa_expected))
})

test_that("utils functions work as they should", {
	expect_true(Biostrings:::isNumericOrNAs(NA_character_))
	expect_true(Biostrings:::isNumericOrNAs(NA_real_))
	expect_true(Biostrings:::isNumericOrNAs(1))
	expect_true(Biostrings:::isNumericOrNAs(1L))
	expect_true(Biostrings:::isNumericOrNAs(1:2))
	expect_false(Biostrings:::isNumericOrNAs(NULL))
	expect_false(Biostrings:::isNumericOrNAs('a'))

	expect_error(Biostrings:::pow.int('a', 1), "must be a numeric vector")
	expect_identical(Biostrings:::pow.int(3,5), as.integer(3**5))

	expect_error(Biostrings:::normargUseNames(NA), "must be TRUE or FALSE")
	expect_true(Biostrings:::normargUseNames(NULL))
	expect_true(Biostrings:::normargUseNames(TRUE))
	expect_false(Biostrings:::normargUseNames(FALSE))
})
