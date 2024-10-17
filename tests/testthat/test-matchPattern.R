## matchPattern.R exports the following:
## - gregexpr2
## - matchPattern generic
## - countPattern generic
## - vmatchPattern generic
## - vcountPattern generic
##
## classes for generics: character, XString,
## 	XStringSet, XStringViews, MaskedXString

sss <- "ATGCATGCATGC"
dna <- DNAString(sss)
dss <- DNAStringSet(list(dna,dna))
maskpos <- c(2L,9L)
m1 <- Mask(length(dna), start=maskpos, end=maskpos+2L)
mdna <- dna
masks(mdna) <- m1
vdna <- Views(dna, start=maskpos, end=maskpos+2L)
# note that gregexpr and gregexpr2 are supported but not documented
all_algos <- c("naive-exact", "naive-inexact", "boyer-moore", "shift-or", "indels")

.test_match_output <- function(subject, pattern, algo, exp_match, ...){
	m <- matchPattern(subject=subject, pattern=pattern, algorithm=algo, ...)
	mc <- countPattern(subject=subject, pattern=pattern, algorithm=algo, ...)
	res <- list(start=start(m), end=end(m), width=width(m), stype=seqtype(m))
	expect_s4_class(m, "XStringViews")
	expect_equal(res, exp_match)
	expect_equal(mc, length(exp_match$start))
}

.test_vmatch_output <- function(subject, pattern, algo, exp_match, ...){
	m <- vmatchPattern(subject=subject, pattern=pattern, algorithm=algo, ...)
	mc <- vcountPattern(subject=subject, pattern=pattern, algorithm=algo, ...)
	res <- lapply(m, \(x) list(start=start(x), end=end(x), width=width(x)))
	expect_s4_class(m, "MIndex")
	expect_equal(res, exp_match)
	expect_equal(mc, vapply(m, \(x) length(start(x)), integer(1L)))
}

test_that("gregexpr2 correctly returns overlapping matches", {
	## should find overlapping matches
	expect_equal(unlist(gregexpr2("aa", c("XaaaYaa", "a"))), c(2L, 3L, 6L, -1L))

	## only supports character vectors of length 1 for 'pattern'
	expect_error(gregexpr2(c("aa", "aaa"), "a"), "invalid pattern")
	expect_error(gregexpr2(NA_character_, "a"), "invalid pattern")
	expect_error(gregexpr2(1, "a"), "invalid pattern")
	expect_error(gregexpr2(character(0L), "a"), "invalid pattern")
})

test_that("matchPattern, countPattern work as expected", {
	## exact matching
	for(algo in all_algos[c(1:4)]){
		.test_match_output(subject=sss, pattern="ATG", algo=algo,
												exp_match=list(start=c(1,5,9), end=c(3,7,11), width=rep(3,3), stype="B"))
		.test_match_output(subject=dna, pattern="ATG", algo=algo,
												exp_match=list(start=c(1,5,9), end=c(3,7,11), width=rep(3,3), stype="DNA"))
		.test_match_output(subject=mdna, pattern="ATG", algo=algo,
												exp_match=list(start=5, end=7, width=3, stype="DNA"))
		.test_match_output(subject=vdna, pattern="ATG", algo=algo,
												exp_match=list(start=9, end=11, width=3, stype="DNA"))
	}

	## inexact matching
	.test_match_output(subject=dna, pattern="ATC", algo="auto",
		exp_match=list(start=c(1,2,5,6,9,10), end=c(3,4,7,8,11,12), width=rep(3,6), stype="DNA"),
		max.mismatch=2)
	.test_match_output(subject=dna, pattern="ATC", algo="auto",
		exp_match=list(start=c(2,6,10), end=c(4,8,12), width=rep(3,3), stype="DNA"),
		max.mismatch=2, min.mismatch=2)
	.test_match_output(subject=dna, pattern="ATC", algo="auto",
		exp_match=list(start=c(1,5,9), end=c(2,6,10), width=rep(2,3), stype="DNA"),
		max.mismatch=1, with.indels=TRUE)

	## fixed args
	dna2 <- DNAString("ACGTNMRWSYKVHDB")
	.test_match_output(subject=dna2, pattern="GTN", algo="auto",
		exp_match=list(start=c(3), end=c(5), width=rep(3,1), stype="DNA"), fixed=TRUE)
	.test_match_output(subject=dna2, pattern="GTN", algo="auto",
		exp_match=list(start=c(3,7,9,12), end=c(5,9,11,14), width=rep(3,4), stype="DNA"), fixed=FALSE)

	## sad paths
	expect_error(matchPattern("ATG", dna, algorithm="indels"), "valid algos for your string matching problem")
	expect_error(matchPattern("ATG", dss), "please use vmatchPattern")
	expect_error(countPattern("ATG", dna, algorithm="indels"), "valid algos for your string matching problem")
	expect_error(countPattern("ATG", dss), "please use vcountPattern")
	expect_error(matchPattern("", dna), "empty patterns are not supported")
	expect_error(matchPattern(1, dna), "'pattern' must be a single string or an XString object")
	expect_error(matchPattern("test", '', algorithm="gregexpr"), "'subject' must be a single (non-empty) string", fixed=TRUE)
	expect_error(matchPattern("test", '', algorithm="gregexpr2"), "'subject' must be a single (non-empty) string", fixed=TRUE)

	## edge cases from matchPattern.R
	## unsure if these are intended functionality, but at least we'll track changes to it
	.test_match_output(subject=DNAString("ACGTGCA"), pattern="---", algo="auto",
		exp_match=list(start=seq(-1,7), end=seq(1,9), width=rep(3,9), stype='DNA'), max.mismatch=3)
	.test_match_output(subject=DNAString("A"), pattern="---", algo="auto",
		exp_match=list(start=integer(0L), end=integer(0L), width=integer(0L), stype='DNA'))
})

test_that("vmatchPattern, vcountPattern work correctly", {
		for(algo in all_algos[c(1:4)]){
			.test_vmatch_output(subject=dss, pattern="ATG", algo=algo,
													exp_match=lapply(seq_len(2), \(x) list(start=c(1,5,9), end=c(3,7,11), width=rep(3,3))))
			.test_vmatch_output(subject=rep(sss,2L), pattern="ATG", algo=algo,
													exp_match=lapply(seq_len(2), \(x) list(start=c(1,5,9), end=c(3,7,11), width=rep(3,3))))
			# XStringViews objects aren't supported yet
			# .test_vmatch_output(subject=vdna, pattern="ATG", algo=algo,
			# 										exp_match=lapply(seq_len(2), \(x) list(start=c(1,5,9), end=c(3,7,11), width=rep(3,3))))
		}

		## sad path
		expect_error(vmatchPattern("ATG", dna), "please use matchPattern")
		expect_error(vmatchPattern("ATG", mdna), "please use matchPattern")
		expect_error(vcountPattern("ATG", dna), "please use countPattern")
		expect_error(vcountPattern("ATG", mdna), "please use countPattern")

		## TODO: more sad path options
		expect_equal(vcountPattern("TG", vdna), c(1L,1L))
})

