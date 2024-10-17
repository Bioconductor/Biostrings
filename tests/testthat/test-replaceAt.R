## testing for both replaceAt.R and replaceLetterAt.R
##
## these export the following:
## - extractAt
## - replaceAt
## - replaceLetterAt
## - (.inplaceReplaceLetterAt)
## These functions work like subseq(<-) but without copying data
## replaceLetterAt will make a copy of the sequence.
## .inplaceReplaceLetterAt will not, but its use is discouraged.
## note that replaceLetterAt is only defined for character, DNAString, DNAStringSet

test_that("replaceLetterAt has expected behavior", {
	str <- "AGAGAGAGAGAG"
	d <- DNAString(str)
	dss <- DNAStringSet(list(d,d))
	pos_replace <- c(1,3,5,7,8)
	mat_replace <- matrix(rep(FALSE, length(d)*2), nrow=2L, byrow=TRUE)
	mat_replace[,pos_replace] <- TRUE
	r_str <- "CGCGCGCCAGAG"
	replace_dss <- DNAStringSet(rep(paste(rep("C", length(pos_replace)), collapse=''), 2))
	#expect_equal(replaceLetterAt(str, pos_replace, "C"), replaced_str)
	expect_equal(replaceLetterAt(d, pos_replace, rep("C", length(pos_replace))), DNAString(r_str))
	expect_equal(replaceLetterAt(d, mat_replace[1,], rep("C", length(pos_replace))), DNAString(r_str))
	expect_equal(replaceLetterAt(dss, mat_replace, replace_dss), DNAStringSet(c(r_str,r_str)))

	d2 <- DNAString("RMWSYK")
	expect_equal(replaceLetterAt(d2, 1:6, rep(c("A","M"), each=3)), DNAString("AAAMMM"))
	expect_equal(replaceLetterAt(d2, 1:6, rep(c("A","M"), each=3), if.not.extending="merge"), DNAString("RMWVHN"))
	expect_equal(replaceLetterAt(d2, 1:6, rep(c("A","M"), each=3), if.not.extending="skip"), d2)
	expect_error(replaceLetterAt(d2, 1:6, rep(c("A","M"), each=3), if.not.extending="error"), "does not extend old letter")
})

## old tests

### From an XStringSet object
### =========================

### Only compare the classes, the shapes (i.e. lengths + elementNROWS +
### names), the inner names, and the sequence contents. Doesn't look at
### the metadata or the metadata columns (outer or inner).
identical_XStringSetList <- function(target, current){
    ok1 <- identical(class(target), class(current))
    ok2 <- identical(elementNROWS(target), elementNROWS(current))
    unlisted_target <- unlist(target, use.names=FALSE)
    unlisted_current <- unlist(current, use.names=FALSE)
    ok3 <- identical(names(unlisted_target), names(unlisted_current))
    ok4 <- all(unlisted_target == unlisted_current)
    ok1 && ok2 && ok3 && ok4
}

### Only compare the classes, lengths, names, and sequence contents.
### Doesn't look at the metadata or the metadata columns.
identical_XStringSet <- function(target, current){
    ok1 <- identical(class(target), class(current))
    ok2 <- identical(names(target), names(current))
    ok3 <- all(target == current)
    ok1 && ok2 && ok3
}

# randomDisjointRanges <- function(min_start, max_end, min_width, max_width, n){
#     offset <- min_start - 1L
#     n1 <- n * (1.1 + max_width/(max_end-offset+1L))  # very approximative

#     some_starts <- sample(max_end-offset+1L, n1, replace=TRUE) + offset
#     some_widths <- sample(min_width:max_width, n1, replace=TRUE)
#     some_ranges <- IRanges(some_starts, width=some_widths)
#     some_ranges <- some_ranges[end(some_ranges) <= max_end]
#     ans <- disjoin(some_ranges)
#     if (min_width == 0L) {
#         is_empty <- width(some_ranges) == 0L
#         some_empty_ranges <- some_ranges[is_empty]
#         ans <- sample(c(ans, some_empty_ranges))
#     }
#     if (length(ans) < n)
#         stop("internal error, sorry")
#     head(ans, n=n)
# }


test_that("old extractAt tests still work", {
	x <- BStringSet(c(seq1="ABCD", seq2="abcdefghijk", seq3="XYZ"))
	at <- IRangesList(IRanges(c(2, 1), c(3, 0)),
	                  IRanges(c(7, 2, 12, 7), c(6, 5, 11, 8)),
	                  IRanges(2, 2))
	### Set inner names on 'at'.
	unlisted_at <- unlist(at)
	names(unlisted_at) <- paste0("rg", sprintf("%02d", seq_along(unlisted_at)))
	at <- relist(unlisted_at, at)

	res1a <- extractAt(x, at)
	res1b <- BStringSetList(mapply(extractAt, x, at))
	expect_true(identical_XStringSetList(res1a, res1b))

	res2a <- extractAt(x, at[3])
	res2b <- BStringSetList(mapply(extractAt, x, at[3]))
	expect_true(identical_XStringSetList(res2a, res2b))
	res2c <- extractAt(x, at[[3]])
	expect_true(identical_XStringSetList(res2a, res2c))

	## performance testing is cut, this is too time consuming
	### Testing performance with 2 millions small extractions at random
	### locations in Fly upstream sequences:
	# dm3_upstream_filepath <- system.file("extdata",
	#                                      "dm3_upstream2000.fa.gz",
	#                                      package="Biostrings")
	# dm3_upstream <- readDNAStringSet(dm3_upstream_filepath)
	# dm3_upstream <- dm3_upstream[width(dm3_upstream) == 2000L]

	# set.seed(33)
	# ## cut this down significantly from 2,000,000
	# ## performance is not important for unit testing
	# NE <- 2000000L  # total number of extractions

	# sample_size <- NE * 1.1
	# some_ranges <- IRanges(sample(2001L, sample_size, replace=TRUE),
	#                        width=sample(0:75, sample_size, replace=TRUE))
	# some_ranges <- head(some_ranges[end(some_ranges) <= 2000L], n=NE)
	# split_factor <- factor(sample(length(dm3_upstream), NE, replace=TRUE),
	#                        levels=seq_along(dm3_upstream))
	# at <- unname(split(some_ranges, split_factor))

	# ### Takes < 1s on my machine.
	# res3a <- extractAt(dm3_upstream, at)

	### Equivalent but about 250x slower than the above on my machine.
	#system.time(res3b <- DNAStringSetList(mapply(extractAt, dm3_upstream, at)))
	#expect_true(identical_XStringSetList(res3a, res3b))
})

test_that("old replaceAt tests still work", {
	### On an XString object
	### ====================

	### Testing performance with half-million small substitutions at random
	### locations in Celegans chrI:
	# library(BSgenome.Celegans.UCSC.ce2)
	# genome <- BSgenome.Celegans.UCSC.ce2
	# chrI <- genome$chrI
	# at4 <- randomDisjointRanges(1L, nchar(chrI), 0L, 20L, 500000L)
	# ### Takes < 1s on my machine (64-bit Ubuntu).
	# system.time(current <- replaceAt(chrI, at4, Views(chrI, at4)))
	# stopifnot(current == chrI)

	# ### Testing performance with half-million single-letter insertions at random
	# ### locations in Celegans chrI:
	# at5 <- randomDisjointRanges(1L, nchar(chrI), 0L, 0L, 500000L)
	# ### Takes < 1s on my machine (64-bit Ubuntu).
	# system.time(current <- replaceAt(chrI, at5, value="-"))
	# m <- matchPattern("-", current)
	# stopifnot(identical(sort(start(at5)), start(m) - seq_along(at5) + 1L))

	# system.time(current2 <- replaceAt(chrI, start(at5), value="-"))
	# stopifnot(identical(current, current2))

	# matchPattern("---", current)

	### On an XStringSet object
	### =======================

	x <- BStringSet(c(seq1="ABCD", seq2="abcdefghijk", seq3="XYZ"))
	at <- IRangesList(IRanges(c(2, 1), c(3, 0)),
	                  IRanges(c(7, 2, 12, 7), c(6, 5, 11, 8)),
	                  IRanges(2, 2))
	### Set inner names on 'at'.
	unlisted_at <- unlist(at)
	names(unlisted_at) <- paste0("rg", sprintf("%02d", seq_along(unlisted_at)))
	at <- relist(unlisted_at, at)

	current <- replaceAt(x, at, value=extractAt(x, at))  # no-op
	expect_true(identical_XStringSet(x, current))

	res1a <- replaceAt(x, at)  # deletions
	res1b <- mendoapply(replaceAt, x, at)
	expect_true(identical_XStringSet(res1a, res1b))

	## cutting this as well
	### Testing performance with half-million single-letter insertions at random
	### locations in Fly upstream sequences:
	##
	# set.seed(33)
	# split_factor <- factor(sample(length(dm3_upstream), 500000L, replace=TRUE),
	#                        levels=seq_along(dm3_upstream))
	# at5 <- unname(split(sample(2001L, 500000L, replace=TRUE),
	#                     split_factor))
	# ### Takes < 1s on my machine.
	# system.time(res5a <- replaceAt(dm3_upstream, at5, value="-"))
	# ### Equivalent but about 1400x slower than the above on my machine.
	# system.time(res5b <- mendoapply(replaceAt,
	#                                 dm3_upstream, as(at5, "List"), as("-", "List")))
	# stopifnot(identical_XStringSet(res5a, res5b))
})