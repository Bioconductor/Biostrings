## PDict objects have the following accessors:
## - length
## - width
## - names
## - head
## - tb
## - tb.width
## - tail
## - (duplicated)
## - (patternFrequency)
## The `[[` operator is also defined

test_that("PDicts can be initialized happy path", {
	codons <- mkAllStrings(DNA_BASES, 3)
	namedcodons <- codons
	names(namedcodons) <- codons

	for(algo in c("ACtree2", "Twobit")){
		expect_s4_class(PDict(codons, algorithm=algo), "TB_PDict")

		p <- PDict(codons, algorithm=algo)
		expect_output(show(p), "TB_PDict object of length 64 and width 3")

		## checking threeparts internals
		expect_s4_class(p@threeparts@pptb, algo)
		expect_equal(lengths(p@threeparts@tail), rep(0L, 64L))
		expect_equal(lengths(p@threeparts@head), rep(0L, 64L))

		## accessors
		expect_equal(length(p), 64L)
		expect_equal(width(p), rep(3L, 64L))

		expect_null(names(p))
		expect_equal(names(PDict(namedcodons, algorithm=algo)), codons)

		expect_equal(duplicated(p), rep(FALSE, 64L))
		expect_equal(patternFrequency(p), rep(1L, 64L))

		expect_equal(duplicated(PDict(c("A","A"), algorithm=algo)), c(FALSE, TRUE))
		expect_equal(patternFrequency(PDict(c("A","A"), algorithm=algo)), c(2L, 2L))


		## tb
		expect_identical(tb(p), p@threeparts@pptb@tb)
		expect_s4_class(tb(p), "DNAStringSet")
		expect_equal(as.character(tb(p)), codons)
		expect_equal(tb.width(p), 3L)

		expect_s4_class(p[[1L]], "DNAString")
		expect_equal(as.character(p[[1L]]), "AAA")

		## initializing with ambiguity codes
		strs <- paste0(DNA_ALPHABET, "AT", DNA_ALPHABET, "A")
		l <- length(DNA_ALPHABET)
		expect_s4_class(PDict(strs, tb.start=2L, tb.width=2L, algorithm=algo), "TB_PDict")
		p <- PDict(strs, tb.start=2L, tb.width=2L, algorithm=algo)
		expect_equal(as.character(tb(p)), rep("AT", l))
		expect_equal(tb.width(p), 2L)
		expect_equal(lengths(head(p)), rep(1L, l))
		expect_identical(as.character(head(p)), DNA_ALPHABET)
		expect_equal(lengths(tail(p)), rep(2L, l))
		expect_identical(as.character(tail(p)), paste0(DNA_ALPHABET, "A"))
		expect_equal(width(p@threeparts), rep(5L, l))

		## variable width tail
		strs[seq(1L,length(strs),2L)] <- paste0(strs[seq(1L,length(strs),2L)], "G")
		p <- PDict(strs, tb.start=2L, tb.width=2L, algorithm=algo)
		expect_equal(lengths(head(p@threeparts)), rep(1L, l))
		expect_identical(as.character(head(p@threeparts)), DNA_ALPHABET)
		expect_equal(lengths(tail(p@threeparts)), rep(c(3L,2L), length.out=l))

		## allowing small number of mismatches
		pats <- apply(expand.grid(mkAllStrings(c("A","T"), 3L), mkAllStrings(c("G","C"), 3L)), 1, paste0, collapse='')
		expect_s4_class(PDict(pats, max.mismatch=1L, algorithm=algo), "MTB_PDict")
		p <- PDict(pats, max.mismatch=1L, algorithm=algo)
		expect_output(show(p), "MTB_PDict object of length 64 and width 6")
		expect_error(p@threeparts, 'no slot of name "threeparts"')
		expect_equal(length(p@threeparts_list), 2L)
		p1 <- p@threeparts_list[[1]]
		p2 <- p@threeparts_list[[2]]
		expect_equal(lengths(p1@head), rep(0L, 64L))
		expect_equal(lengths(p1@tail), rep(3L, 64L))
		expect_equal(as.character(p1@tail), substr(pats,4,6))
		expect_equal(lengths(p2@head), rep(3L, 64L))
		expect_equal(lengths(p2@tail), rep(0L, 64L))
		expect_equal(as.character(p2@head), substr(pats,1,3))

		expect_equal(length(as.list(p)), 2L)
	}
})

test_that("PDict sad path errors correctly", {
	expect_error(PDict(DNA_ALPHABET), "non base DNA letter found in Trusted Band")
	expect_error(PDict(AA_ALPHABET), "not in lookup table")
	expect_error(PDict(1:3), "unable to find an inherited method")
	expect_error(PDict(character(0L)), "must contain at least one pattern")
	p <- PDict(DNA_BASES)
	expect_error({names(p) <- DNA_BASES}, "attempt to modify the names of a TB_PDict instance")
	expect_error(PDict(c("AA", "AAA")), "Trusted Band has a different length")

	codons <- mkAllStrings(DNA_BASES, 3)
	expect_error(PDict(codons, max.mismatch=1), "supported only if the width of dictionary 'x' is >= 6")

	longstrings <- mkAllStrings(c("A","T"), 6)
	expect_error(PDict(longstrings, max.mismatch=7), "must be <= 1 given the width of dictionary 'x'")
	expect_error(PDict(paste0(rep("A", 100), collapse=''), max.mismatch=100), "must be <= 32 given the width of dictionary 'x'")

	expect_error(PDict(codons, max.mismatch=NULL), "must be a single integer or 'NA'")
	expect_error(PDict(codons, max.mismatch='a'), "must be a single integer or 'NA'")
	expect_error(PDict(codons, max.mismatch=c(1,2)), "must be a single integer or 'NA'")

	expect_error(PDict("A", algorithm='erroralgo'), "is not a defined class")

	names(longstrings) <- rep('', length(longstrings))
	expect_error(PDict(longstrings), "'x' has invalid names")
	names(longstrings) <- rep('test', length(longstrings))
	expect_error(PDict(longstrings), "'x' has duplicated names")

	expect_error(PDict("AAAA", tb.start='a'), "must be a single integer or 'NA'")
	expect_error(PDict("AAAA", tb.end='a'), "must be a single integer or 'NA'")
	expect_error(PDict("AAAA", tb.width='a'), "must be a single integer or 'NA'")

	expect_warning(PDict("AAAA", algorithm='ACtree'), "support for ACtree preprocessing algo has been dropped")
})

### Old Tests

## Helper functions ported from old tests
randomDNASequences <- function(n, w)
{
  alphabet <- DNA_BASES
  w <- rep(w, length=n)
  sequences <- sapply(seq(1, n, length=n),
                      function(x) {
                        s <- sample(alphabet, w[x], replace=TRUE)
                        s <- paste(s, collapse="")
                        return(s)
                      })
  return(Biostrings::DNAStringSet(sequences))
}

msubseq <- function(x, ir)
{
  ## differs from subseq in the sense that several subsequences
  ## from the same sequence are extracted
  ## x:  XString
  ## ir: IRanges
  res <- vector("character", length = length(ir))
  for (i in seq(along=res)) {
    res[i] <- as.character(subseq(x, start=ir@start[i], width=width(ir)[i]))
    ## forced cast: chain of tools for DNAString seems interupted for
    ##              some use cases (or I missed something)
  }
  res <- DNAStringSet(res)
  return(res)
}

test_that("PDict works with constant width initialization", {
  set.seed(1)
  l <- 150
  target <- randomDNASequences(1, l)[[1]]
  W <- 20
  L <- 6
  ir <- successiveIRanges(rep(W, L), gapwidth = 1)
  short_sequences <- msubseq(target, ir)
  # shuffle the sequences (they are not in consecutive order)
  o <- sample(seq(along=short_sequences))

  dna_short <- DNAStringSet(short_sequences[o])
  pdict <- PDict(dna_short)
  expect_equal(L, length(pdict))
  expect_equal(rep(W, L), width(pdict))
  expect_equal(NULL, head(pdict))
  expect_equal(W, tb.width(pdict))
  expect_equal(NULL, tail(pdict))
})

test_that("PDict works for variable width lookup", {
  set.seed(1)
  l <- 150
  target <- randomDNASequences(1, l)[[1]]
  W <- 20
  L <- 6
  n_cut <- sample(0:5, L, replace=TRUE)
  ir <- successiveIRanges(rep(W, L) - n_cut, gapwidth = 1 + n_cut[-length(n_cut)])
  short_sequences <- msubseq(target, ir)
  # shuffle the sequences (they are not in consecutive order)
  o <- sample(seq(along=short_sequences))

  dna_var_short <- DNAStringSet(short_sequences[o])

  ## Previous comment: shouldn't 1:min(width) be the default?
  pdict <- PDict(dna_var_short,
                 tb.start=1,
                 tb.width=min(width(short_sequences))
                 )
  expect_equal(L, length(pdict))
  expect_equal((rep(W, L) - n_cut)[o], width(pdict))
  expect_equal(NULL, head(pdict))
  shortest_seq_width <- min(width(dna_var_short))
  expect_equal(shortest_seq_width,
              tb.width(pdict))           # mostly a sanity check
  expect_equal(substring(short_sequences, shortest_seq_width+1)[o],
              as.character(tail(pdict)))
})