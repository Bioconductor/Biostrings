## translate.R exports the following:
## - translate (DNAString, RNAString, DNAStringSet, RNAStringSet,
##								MaskedDNAString, MaskedRNAString, Default)
## - codons (DNAString, RNAString, MaskedDNAString, MaskedRNAString, Default)

d <- DNAString(paste(names(GENETIC_CODE), collapse=''))
r <- as(d, "RNAString")
a <- AAString(paste(GENETIC_CODE, collapse=''))

test_that("translate works properly for all dna,rna", {
	## basic conversion
	expect_identical(translate(d), a)
	expect_identical(translate(r), a)

	dd <- DNAStringSet(list(d,d))
	rr <- RNAStringSet(list(r,r))
	aa <- AAStringSet(list(a,a))
	expect_identical(translate(dd) == aa, c(TRUE, TRUE))
	expect_identical(translate(rr) == aa, c(TRUE, TRUE))

	## errors
	expect_error(translate(BString()), "unable to find an inherited method")
	expect_error(translate(AAString()), "unable to find an inherited method")
	expect_error(translate("ATG"), "unable to find an inherited method")

	alt_gen_code <- GENETIC_CODE
	alt_gen_code[] <- rep(c("A", "B"), each=32L)
	alt_a <- AAString(paste(rep(c("A", "B"), each=32L), collapse=''))
	alt_aa <- AAStringSet(list(alt_a, alt_a))
	expect_identical(translate(d, genetic.code=alt_gen_code), alt_a)
	expect_identical(translate(r, genetic.code=alt_gen_code), alt_a)
	expect_identical(translate(dd, genetic.code=alt_gen_code) == alt_aa, c(TRUE, TRUE))
	expect_identical(translate(rr, genetic.code=alt_gen_code) == alt_aa, c(TRUE, TRUE))

	## fault paths
	expect_error(translate(d, genetic.code=character(0L)),
		"'genetic.code' must have the same names as predefined constant GENETIC_CODE", fixed=TRUE)
	expect_error(translate(d, genetic.code=GENETIC_CODE[seq_len(10L)]),
		"'genetic.code' must have the same names as predefined constant GENETIC_CODE", fixed=TRUE)
	alt_gen_code[] <- ''
	expect_error(translate(d, genetic.code=alt_gen_code),
		"'genetic.code' must contain 1-letter strings", fixed=TRUE)

	## this is promoted to an error with this test case
	alt_gen_code[] <- '/'
	expect_error(translate(d, genetic.code=alt_gen_code),
		"some codons in 'genetic.code' are mapped to letters not in the Amino Acid\n  alphabet", fixed=TRUE)

	## alt_init_codon args
	## error message could be more informative here
	error_msg <- "'genetic.code' must have an \"alt_init_codons\" attribute"
	alt_gen_code <- GENETIC_CODE
	attr(alt_gen_code, "alt_init_codons")[] <- c("BLAH", "NOTANAMINOACID")
	expect_error(translate(d, genetic.code=alt_gen_code), error_msg, fixed=TRUE)
	attr(alt_gen_code, "alt_init_codons")[] <- c("AAA", "AAA")
	expect_error(translate(d, genetic.code=alt_gen_code), error_msg, fixed=TRUE)
	attr(alt_gen_code, "alt_init_codons")[] <- c(NA_character_, "AAA")
	expect_error(translate(d, genetic.code=alt_gen_code), error_msg, fixed=TRUE)
	attr(alt_gen_code, "alt_init_codons")[] <- c(1,2)
	expect_error(translate(d, genetic.code=alt_gen_code), error_msg, fixed=TRUE)
	attr(alt_gen_code, "alt_init_codons") <- NULL
	expect_error(translate(d, genetic.code=alt_gen_code), error_msg, fixed=TRUE)

	## alternate init codons
	dna1 <- DNAString("TTGATATGGCCCTTATAA")
	expect_identical(translate(dna1), AAString("MIWPL*"))
	expect_identical(translate(dna1, no.init.codon=TRUE), AAString("LIWPL*"))
	rna1 <- as(dna1, "RNAString")
	expect_identical(translate(rna1), AAString("MIWPL*"))
	expect_identical(translate(rna1, no.init.codon=TRUE), AAString("LIWPL*"))

	## fuzzy codons
	dna1 <- DNAString("TTYTTTTTC")
	expect_error(translate(dna1), "not a base at pos 3")
	expect_identical(translate(dna1, if.fuzzy.codon='solve'), AAString("FFF"))
	expect_identical(translate(dna1, if.fuzzy.codon='X'), AAString("XFF"))
	expect_identical(translate(dna1, if.fuzzy.codon='error.if.X'), AAString("XFF"))

	dna2 <- DNAString("TTNTTYTTC")
	expect_identical(translate(dna2, if.fuzzy.codon='solve'), AAString("XFF"))
	expect_identical(translate(dna2, if.fuzzy.codon='X'), AAString("XXF"))
	expect_error(translate(dna2, if.fuzzy.codon='error.if.X'), "ambiguous fuzzy codon starting at pos 1")

	expect_error(translate(dna2, if.fuzzy.codon='blahblahblaherror'), '\'arg\' should be one of ', fixed=TRUE)
})

test_that("translate, codons work properly for masked strings and string sets", {
	mask0 <- Mask(mask.width=nchar(d), start=c(4,10), width=3)
	mask1 <- Mask(mask.width=nchar(d), start=c(4,10), width=4)
	mask2 <- Mask(mask.width=nchar(d), start=1, width=5)

	dm <- d
	rm <- r
	masks(dm) <- masks(rm) <- mask0

	m_a <- AAString(paste(GENETIC_CODE[-c(2,4)], collapse=''))
	expect_identical(translate(dm), m_a)
	expect_identical(translate(rm), m_a)

	codon_views_d <- codons(dm)
	codon_views_r <- codons(rm)
	svec <- c(1L, 7L, seq(13L,nchar(d),3L))
	expect_identical(start(codon_views_d), svec)
	expect_identical(start(codon_views_r), svec)
	expect_identical(end(codon_views_d), svec+2L)
	expect_identical(end(codon_views_r), svec+2L)

	masks(d) <- mask1
	expect_warning(translate(d), "last base was ignored")

})

test_that("codon function works properly", {
	expect_identical(codons(d), Views(d, start=seq(1,nchar(d),3), width=3))
	expect_identical(codons(r), Views(r, start=seq(1,nchar(d),3), width=3))

	expect_warning(codons(DNAString("AT")), "the number of nucleotides in 'x' is not a multiple of 3")
	expect_warning(codons(DNAString("ATATA")), "the number of nucleotides in 'x' is not a multiple of 3")
	expect_warning(codons(RNAString("AU")), "the number of nucleotides in 'x' is not a multiple of 3")
	expect_warning(codons(RNAString("AUAUA")), "the number of nucleotides in 'x' is not a multiple of 3")

	expect_error(codons(AAString()), "unable to find an inherited method for function")
	expect_error(codons(BString()), "unable to find an inherited method for function")
	expect_error(codons("ATG"), "unable to find an inherited method for function")
})