## Miscellaneous tests
## these are all relatively low priority and/or for files with only a few things




test_that("coloring works for DNA, RNA, and AA", {
	## not a super important test
	dna_rna_expected <- c(DNA_BASES, "U", DNA_ALPHABET[-c(1:4,16:18)])
	expect_true(!any(duplicated(Biostrings:::make_DNA_AND_RNA_COLORED_LETTERS())))
	expect_equal(sort(names(Biostrings:::make_DNA_AND_RNA_COLORED_LETTERS())), sort(dna_rna_expected))

	aa_expected <- AA_ALPHABET[-c(27:30)]
	expect_true(!any(duplicated(Biostrings:::make_AA_COLORED_LETTERS())))
	expect_equal(sort(names(Biostrings:::make_AA_COLORED_LETTERS())), sort(aa_expected))
})

