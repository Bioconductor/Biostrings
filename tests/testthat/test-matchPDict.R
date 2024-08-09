## matchPDict.R exports the following:
##		- matchPDict
##		- countPDict
##		- whichPDict
##		- vmatchPDict
##		- vcountPDict
##		- vwhichPDict
## These are all S4s defined for "subject", XString,
##	XStringSet, XStringViews, MaskedXString

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

test_that("inexact PDict matching works", {
	pdict <- PDict(c("acgt", "gt", "cgt", "ac"), tb.end=2)
	d <- DNAString("acggaccg")
  expect_equal(lengths(endIndex(matchPDict(pdict, d, max.mismatch=0))), c(0L,0L,0L,2L))
  expect_equal(lengths(endIndex(matchPDict(pdict, d, max.mismatch=1))), c(1L,0L,2L,2L))
  expect_equal(lengths(endIndex(matchPDict(pdict, d, max.mismatch=2))), c(2L,0L,2L,2L))
})

### All old tests
test_that("PDict works with constant width lookup", {
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
  expect_equal( (rep(W, L) - n_cut)[o], width(pdict))
  expect_equal(NULL, head(pdict))
  shortest_seq_width <- min(width(dna_var_short))
  expect_equal(shortest_seq_width,
              tb.width(pdict))           # mostly a sanity check
  expect_equal(substring(short_sequences, shortest_seq_width+1)[o],
              as.character(tail(pdict)))
})


test_that("matchPDict works for constant width lookup", {
  set.seed(1)
  l <- 150
  dna_target <- randomDNASequences(1, l)[[1]]
  W <- 20
  L <- 6
  ir <- successiveIRanges(rep(W, L), gapwidth = 1)
  short_sequences <- msubseq(dna_target, ir)
  # shuffle the sequences (so they are not in consecutive order)
  o <- sample(seq(along=short_sequences))

  dna_short <- DNAStringSet(short_sequences[o])
  pdict <- PDict(dna_short)

  res <- matchPDict(pdict, dna_target)

  # mostly a sanity check
  expect_equal(L, length(res))

  for (i in seq(along=res)) {
    m_start <- ir[o][i]@start
    expect_equal(m_start, start(res[[i]]))
    expect_equal(W, width(res[[i]]))
    expect_equal(m_start + W - 1, end(res[[i]]))  # mostly a sanity check
  }
})

test_that("matchPDict works for variable width lookup", {
  set.seed(1)
  l <- 150
  dna_target <- randomDNASequences(1, l)[[1]]
  W <- 20
  L <- 6
  n_cut <- sample(0:5, L, replace=TRUE)
  ir <- successiveIRanges(rep(W, L) - n_cut, gapwidth = 1 + n_cut[-length(n_cut)])
  short_sequences <- msubseq(dna_target, ir)
  # shuffle the sequences (they are not in consecutive order)
  o <- sample(seq(along=short_sequences))

  dna_var_short <- DNAStringSet(short_sequences[o])

  pdict <- PDict(dna_var_short,
                 tb.start=1,                        # can't this be
                 tb.width=min(width(dna_var_short)) # the default for
                                                    # variable width ?
                 )

  res <- matchPDict(pdict, dna_target)

  # mostly a sanity check
  expect_equal(L, length(res))

  iro <- ir[o]
  for (i in seq(along=res)) {
    expect_equal(start(iro)[i], start(res[[i]]))
    expect_equal(width(iro)[i], width(res[[i]]))
    expect_equal(end(iro)[i], end(res[[i]]))  # mostly a sanity check
  }
})

