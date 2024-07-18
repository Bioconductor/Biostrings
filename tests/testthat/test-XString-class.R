dnastr <- paste(DNA_ALPHABET, collapse='')
rnastr <- paste(RNA_ALPHABET, collapse='')
aastr <- paste(AA_ALPHABET, collapse='')
bstr <- rawToChar(as.raw(32:126))

test_allstrings_for_error <- function(input, exp_error, ignore=''){
	if(!"DNA"%in%ignore) expect_error(DNAString(input), exp_error)
	if(!"RNA"%in%ignore) expect_error(RNAString(input), exp_error)
	if(!"AA"%in%ignore)  expect_error(AAString(input), exp_error)
	if(!"B"%in%ignore)   expect_error(BString(input), exp_error)
}

test_that("seqtype correctly infers types", {
	expect_equal(seqtype(BString("ABC")), "B")
	expect_equal(seqtype(RNAString("AUG")), "RNA")
	expect_equal(seqtype(DNAString("ATG")), "DNA")
	expect_equal(seqtype(AAString("ARN")), "AA")
})


test_that("character, vector conversion works properly", {
	expect_equal(as.character(DNAString(dnastr)), dnastr)
	expect_equal(as.character(RNAString(rnastr)), rnastr)
	expect_equal(as.character(AAString(aastr)), aastr)
	expect_equal(as.character(BString(bstr)), bstr)
})

test_that("encode/decode tables work correctly", {
	expect_equal(DNAString(tolower(dnastr)), DNAString(dnastr))
	expect_equal(RNAString(tolower(rnastr)), RNAString(rnastr))
	expect_equal(AAString(tolower(aastr)), AAString(aastr))

	expect_equal(DNAString(as.factor(dnastr)), DNAString(dnastr))
	expect_equal(RNAString(as.factor(rnastr)), RNAString(rnastr))
	expect_equal(AAString(as.factor(aastr)), AAString(aastr))
	expect_equal(BString(as.factor(bstr)), BString(bstr))

	test_allstrings_for_error(bstr, "not in lookup table", ignore="B")
})

test_that("'letter' correctly extracts elements", {
	expect_equal(letter(DNAString("ATGCATGC"), c(1,3,6)), "AGT")
	expect_equal(letter(RNAString("AUGCAUGC"), c(1,3,6)), "AGU")
	expect_equal( letter(AAString("ARNDARND"), c(1,3,6)), "ANR")
	expect_equal(  letter(BString("ABCDEFGH"), c(1,3,6)), "ACF")
	expect_equal(letter(bstr, seq(nchar(bstr), 1)), rawToChar(rev(charToRaw(bstr))))

	## TODO: Better error messages
	expect_error(letter(DNAString(""), 10), "> length\\(x\\)")
	expect_error(letter(RNAString(""), 10), "> length\\(x\\)")
	expect_error(letter(AAString(""), 10), "> length\\(x\\)")
	expect_error(letter(BString(""), 10), "> length\\(x\\)")
	expect_error(letter("", 10), "out of bounds")
})

test_that("constructors handle invalid and non-char input correctly", {
	# note that this only displays for character input
	# meaning DNAString(1) returns a non-informative error
	test_allstrings_for_error(NA_character_, "input must be a single non-NA string")
	test_allstrings_for_error(as.factor(1:4), "input must be a single non-NA string")
	test_allstrings_for_error(c("A", "A"), "input must be a single non-NA string")
})

test_that("conversion between XStrings works properly", {
	# DNA <-> RNA
	expect_equal(RNAString(DNAString(dnastr)), RNAString(rnastr))
	expect_equal(DNAString(RNAString(rnastr)), DNAString(dnastr))

	# X -> B
	expect_equal(BString(DNAString(dnastr)), BString(dnastr))
	expect_equal(BString(RNAString(rnastr)), BString(rnastr))
	expect_equal(BString(AAString(aastr)), BString(aastr))

	# valid B -> X
	expect_equal(DNAString(BString(dnastr)), DNAString(dnastr))
	expect_equal(RNAString(BString(rnastr)), RNAString(rnastr))
	expect_equal(AAString(BString(aastr)), AAString(aastr))

	# invalid DNA,RNA <-> AA
	expect_error(AAString(DNAString("ATGC")), "incompatible sequence types")
	expect_error(AAString(RNAString("AUGC")), "incompatible sequence types")
	expect_error(DNAString(AAString("ARND")), "incompatible sequence types")
	expect_error(RNAString(AAString("ARND")), "incompatible sequence types")

	# invalid B -> X
	bbad <- BString(";")
	test_allstrings_for_error(bbad, 'not in lookup table', ignore='B')
})

test_that("as.vector methods work correctly", {
	## TODO: sometimes returns factor, sometimes character
	## 				note the BString case as an example
	dnafac <- factor(seq_len(nchar(dnastr)))
	attr(dnafac, 'levels') <- strsplit(dnastr, '')[[1]]

	rnafac <- factor(seq_len(nchar(rnastr)))
	attr(rnafac, 'levels') <- strsplit(rnastr, '')[[1]]

	aafac <- factor(seq_len(nchar(aastr)))
	attr(aafac, 'levels') <- strsplit(aastr, '')[[1]]

	bfac <- factor(seq_len(nchar(bstr)))
	attr(bfac, 'levels') <- strsplit(bstr, '')[[1]]

	expect_equal(as.vector(DNAString(dnastr)), dnafac)
	expect_equal(as.vector(RNAString(rnastr)), rnafac)
	expect_equal(as.vector(AAString(aastr)), aafac)
	#expect_equal(as.vector(BString(bstr)), bfac)
})

test_that("equality methods work as advertised", {
	expect_equal(DNAString(dnastr) == DNAString(dnastr), TRUE)
	expect_equal(RNAString(rnastr) == RNAString(rnastr), TRUE)
	expect_equal(AAString(aastr) == AAString(aastr), TRUE)
	expect_equal(BString(bstr) == BString(bstr), TRUE)

	expect_equal(DNAString(dnastr) == DNAString(substr(dnastr, 1, 7)), FALSE)
	expect_equal(RNAString(rnastr) == RNAString(substr(rnastr, 1, 7)), FALSE)
	expect_equal(AAString(aastr) == AAString(substr(aastr, 1, 7)), FALSE)
	expect_equal(BString(bstr) == BString(substr(bstr, 1, 7)), FALSE)

	expect_equal(DNAString(dnastr) != DNAString(substr(dnastr, 1, 7)), TRUE)
	expect_equal(RNAString(rnastr) != RNAString(substr(rnastr, 1, 7)), TRUE)
	expect_equal(AAString(aastr) != AAString(substr(aastr, 1, 7)), TRUE)
	expect_equal(BString(bstr) != BString(substr(bstr, 1, 7)), TRUE)


	## DNA <-> RNA comparison
	expect_equal(DNAString(dnastr) == RNAString(rnastr), TRUE)

	## other B comparisons
	expect_equal(AAString(aastr) == BString(aastr), TRUE)
	expect_equal(AAString(aastr) != BString(bstr), TRUE)
	expect_equal(bstr == BString(bstr), TRUE)

	## invalid comparisons
	expect_error(DNAString(dnastr) == AAString(aastr), 'comparison between a "DNAString" instance and a "AAString" instance is not supported')
	expect_error(DNAString(dnastr) == BString(bstr), 'comparison between a "DNAString" instance and a "BString" instance is not supported')
	expect_error(RNAString(rnastr) == AAString(aastr), 'comparison between a "RNAString" instance and a "AAString" instance is not supported')
	expect_error(RNAString(rnastr) == BString(bstr), 'comparison between a "RNAString" instance and a "BString" instance is not supported')
})

test_that("output works correctly", {
	s <- paste(rep("A", 100L), collapse='')
	expect_output(show(DNAString(s)), "100-letter DNAString object\\nseq: .+\\.\\.\\..+$", width=80)
	expect_output(show(RNAString(s)), "100-letter RNAString object\\nseq: .+\\.\\.\\..+$", width=80)
	expect_output(show(AAString(s)), "100-letter AAString object\\nseq: .+\\.\\.\\..+$", width=80)
	expect_output(show(BString(s)), "100-letter BString object\\nseq: A{36}\\.\\.\\.A{36}$", width=80)

	# width of sequence is 5 less than that of 'width'
	expect_output(show(DNAString(s)), "100-letter DNAString object\\nseq: .+\\.\\.\\..+$", width=10)
	expect_output(show(RNAString(s)), "100-letter RNAString object\\nseq: .+\\.\\.\\..+$", width=10)
	expect_output(show(AAString(s)), "100-letter AAString object\\nseq: .+\\.\\.\\..+$", width=10)
	expect_output(show(BString(s)), "100-letter BString object\\nseq: AA\\.\\.\\.AA$", width=10)
})

test_that("substr, substring methods work correctly", {
	d <- DNAString(dnastr)
	r <- RNAString(rnastr)
	a <- AAString(aastr)
	b <- BString(bstr)

	expect_equal(as.character(substr(d, 1, 10)), substr(dnastr, 1, 10))
	expect_equal(as.character(substr(r, 1, 10)), substr(rnastr, 1, 10))
	expect_equal(as.character(substr(a, 1, 10)), substr(aastr, 1, 10))
	expect_equal(as.character(substr(b, 1, 10)), substr(bstr, 1, 10))

	expect_equal(as.character(substring(d, 5, 10)), substring(dnastr, 5, 10))
	expect_equal(as.character(substring(r, 5, 10)), substring(rnastr, 5, 10))
	expect_equal(as.character(substring(a, 5, 10)), substring(aastr, 5, 10))
	expect_equal(as.character(substring(b, 5, 10)), substring(bstr, 5, 10))

	expect_error(substring(d, 10, 5), "Invalid sequence coordinates")
	expect_error(substring(r, 10, 5), "Invalid sequence coordinates")
	expect_error(substring(a, 10, 5), "Invalid sequence coordinates")
	expect_error(substring(b, 10, 5), "Invalid sequence coordinates")

	# `[` dispatch
	expect_equal(as.character(d[1:10]), substr(dnastr, 1, 10))
	expect_equal(as.character(d[-1]), substr(dnastr, 2, nchar(dnastr)))
})

## Porting RUnit tests
test_that("alphabet finds the correct values", {
	expect_equal(DNAString(dnastr), strsplit(dnastr, '')[[1]])
	expect_equal(RNAString(rnastr), strsplit(rnastr, '')[[1]])
	expect_equal(AAString(aastr), strsplit(aastr, '')[[1]])
	expect_equal(BString(bstr), NULL)
})