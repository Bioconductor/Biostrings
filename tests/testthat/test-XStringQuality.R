## XStringQuality-class.R exports the following:
## - PhredQuality
## - SolexaQuality
## - IlluminaQuality
## - offset
## - alphabet
## - encoding
## - coercion methods

## Phred: 0-99 quality as (ASCII)-33 (! = 0, note that this can't support >=95)
## Solexa: -5-99 quality as (ASCII)-64 (; = 0, no valid characters above 63)
## Illumina: 0-99 quality as (ASCII)-64 (@ = 0, no valid characters above 63)

test_that("Quality constructors work correctly", {
	expect_s4_class(PhredQuality(seq(0,99)), "PhredQuality")
	expect_s4_class(SolexaQuality(seq(-5,99)), "SolexaQuality")
	expect_s4_class(IlluminaQuality(0:99), "IlluminaQuality")
	p <- PhredQuality(seq(0,99))
	s <- SolexaQuality(seq(-5,99))
	i <- IlluminaQuality(0:99)

	## note that PhredQuality cannot represent scores greater than 93
	## (and it auto-rounds all these values down)
	## all other qualities can support the values in their range
	complete_ascii <- strsplit(rawToChar(as.raw(33:126)), '')[[1]]

	## check alphabets
	expect_equal(alphabet(p), complete_ascii[seq_len(94)])
	expect_equal(alphabet(s)[seq_len(68)], complete_ascii[seq(27,94)])
	expect_equal(alphabet(i)[seq_len(63)], complete_ascii[seq(32,94)])

	## check encodings
	mapping <- seq(0,131)
	names(mapping) <- complete_ascii
	expect_equal(encoding(p), mapping[seq_len(94)])
	expect_equal(encoding(s), mapping[seq(27,131)] - 31L)
	expect_equal(encoding(i), mapping[seq(32,131)] - 31L)

	## check probability conversion
	prob_conversions <- c("numeric", "NumericList")

	## probability definitions taken from the man pages, should match
	phred_scores <- c(0:93, rep(93,6L))
	solexa_scores <- seq(-5L, 99L)
	illumina_scores <- seq(0L, 99L)

	phred_probs <- 10^(-0.1 * phred_scores)
	solexa_probs <- 10^(-0.1*solexa_scores) / (1 + 10^(-0.1*solexa_scores))
	illumina_probs <- 10^(-0.1 * illumina_scores)
	expect_equal(as.numeric(p), phred_probs)
	expect_equal(as.numeric(s), solexa_probs)
	expect_equal(as.numeric(i), illumina_probs)
	expect_equal(as.numeric(as(p, "NumericList")[[1]]), phred_probs)
	expect_equal(as.numeric(as(s, "NumericList")[[1]]), solexa_probs)
	expect_equal(as.numeric(as(i, "NumericList")[[1]]), illumina_probs)

	## score conversions
	score_conversions <- c("integer", "matrix", "IntegerList")
	expect_equal(as.integer(p), phred_scores)
	expect_equal(as.integer(s), solexa_scores)
	expect_equal(as.integer(i), illumina_scores)
	expect_equal(as(p,'matrix')[1,], phred_scores)
	expect_equal(as(s,'matrix')[1,], solexa_scores)
	expect_equal(as(i,'matrix')[1,], illumina_scores)
	expect_equal(as.integer(as(p,'IntegerList')[[1]]), phred_scores)
	expect_equal(as.integer(as(s,'IntegerList')[[1]]), solexa_scores)
	expect_equal(as.integer(as(i,'IntegerList')[[1]]), illumina_scores)

	## BString, character conversions
	b <- BString(paste(LETTERS, collapse=''))
	expect_equal(as.character(PhredQuality(b)), as.character(b))
	expect_equal(as.character(SolexaQuality(b)), as.character(b))
	expect_equal(as.character(IlluminaQuality(b)), as.character(b))
	expect_equal(as.character(PhredQuality(as.character(b))), as.character(b))
	expect_equal(as.character(SolexaQuality(as.character(b))), as.character(b))
	expect_equal(as.character(IlluminaQuality(as.character(b))), as.character(b))
})