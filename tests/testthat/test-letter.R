## letter.R exports `letter` for the following classes:
## - character
## - XString
## - XStringViews
## - MaskedXString


test_that("letter generic works properly for all classes", {
	s <- paste(letters, collapse='')
	sdna <- sample(DNA_BASES, 20L, replace=TRUE)

	dna <- DNAString(paste(sdna, collapse=''))
	viewdna <- Views(dna, start=rep(1L,3L), width=10)

	m <- Mask(mask.width=20L, start=1, end=5L)
	maskdna <- dna
	masks(maskdna) <- m

	## happy path testing
	test_cases <- list(1.0, 1:3, 1:3*2L, integer(0L))
	for(i in seq_along(test_cases)){
		## character
		expect_equal(letter(s, test_cases[[i]]), paste(letters[test_cases[[i]]], collapse=''))

		## XString
		expect_equal(letter(dna, test_cases[[i]]), paste(sdna[test_cases[[i]]], collapse=''))

		## XStringViews
		expect_equal(letter(viewdna, test_cases[[i]]), rep(paste(sdna[test_cases[[i]]], collapse=''), 3L))

		## MaskedXString -- note that this ignores the mask
		expect_equal(letter(maskdna, test_cases[[i]]), paste(sdna[test_cases[[i]]], collapse=''))
	}

	## sad path testing
	sad_list <- list('a', NA_real_, 10000)
	error_msgs <- c("must be an NA-free numeric vector",
									"must be an NA-free numeric vector",
									"out of bounds")
	for(i in seq_along(sad_list)){
		expect_error(letter(s, sad_list[[i]]), error_msgs[i])
		expect_error(letter(dna, sad_list[[i]]), error_msgs[i])
		expect_error(letter(viewdna, sad_list[[i]]), error_msgs[i])
		expect_error(letter(maskdna, sad_list[[i]]), error_msgs[i])
	}
})