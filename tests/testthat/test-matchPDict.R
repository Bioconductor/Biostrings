## matchPDict.R exports the following:
##		- matchPDict
##		- countPDict
##		- whichPDict
##		- vmatchPDict
##		- vcountPDict
##		- vwhichPDict
## These are all S4s defined for "subject", XString,
##	XStringSet, XStringViews, MaskedXString

test_that("general match/count/whichPDict behavior works as expected", {
  d1 <- DNAString("ATGGATACACAA")
  all_codons <- mkAllStrings(DNA_BASES, 3L)

  ## non-pdict with matchPDict
  p1 <- matchPDict(DNAStringSet("ATA"), d1)
  expect_equal(as.integer(c(start(p1), end(p1))), c(5,7))

  pd <- PDict(all_codons)

  ## pd against a DNAStringSet
  p2 <- lengths(matchPDict(pd, d1))
  pos_good <- which(all_codons %in% as.character(Views(d1, start=seq_len(length(d1)-2L), width=3)))
  expect_true(all(p2[pos_good] > 0) && all(p2[-pos_good] == 0))
  expect_equal(whichPDict(pd, d1), pos_good)
  expect_equal(countPDict(pd, d1), p2)

  ## pd against an XStringViews object
  p3 <- lengths(matchPDict(pd, codons(d1)))
  pos_good <- which(all_codons %in% as.character(codons(d1)))
  expect_true(all(p3[pos_good] > 0) && all(p3[-pos_good] == 0))
  expect_equal(whichPDict(pd, codons(d1)), pos_good)
  expect_equal(countPDict(pd, codons(d1)), p3)

  ## pd against a MaskedXString object
  m1 <- Mask(length(d1), start=4, end=8)
  d2 <- d1
  masks(d2) <- m1
  p4 <- lengths(matchPDict(pd, d2))
  pos_good <- which(all_codons %in% as.character(Views(unmasked(d2), start=c(1,9,10), width=3)))
  expect_true(all(p4[pos_good] > 0) && all(p4[-pos_good] == 0))
  expect_equal(whichPDict(pd, d2), pos_good)
  expect_equal(countPDict(pd, d2), p4)


  ## matchPDict() with AA, RNA, B
  aa1 <- AAStringSet(c("DARC", "EGH"))
  aa2 <- AAString("KMFPRNDEGHSTTWTEE")
  paa <- matchPDict(aa1, aa2)
  expect_equal(c(unlist(startIndex(paa)), unlist(endIndex(paa))), c(8,10))
  expect_equal(countPDict(aa1, aa2), c(0,1))
  expect_equal(whichPDict(aa1, aa2), 2)

  b1 <- BStringSet(c("DARC", "EGH"))
  b2 <- BString("KMFPRNDEGHSTTWTEE")
  pb <- matchPDict(b1, b2)
  expect_equal(c(unlist(startIndex(pb)), unlist(endIndex(pb))), c(8,10))
  expect_equal(countPDict(b1, b2), c(0,1))
  expect_equal(whichPDict(b1, b2), 2)
})

test_that("inexact PDict matching works", {
	pdict <- PDict(c("acgt", "gt", "cgt", "ac"), tb.end=2)
	d <- DNAString("acggaccg")
  expect_equal(lengths(endIndex(matchPDict(pdict, d, max.mismatch=0))), c(0L,0L,0L,2L))
  expect_equal(lengths(endIndex(matchPDict(pdict, d, max.mismatch=1))), c(1L,0L,2L,2L))
  expect_equal(lengths(endIndex(matchPDict(pdict, d, max.mismatch=2))), c(2L,0L,2L,2L))

  expect_equal(countPDict(pdict, d, max.mismatch=0), c(0,0,0,2))
  expect_equal(countPDict(pdict, d, max.mismatch=1), c(1,0,2,2))
  expect_equal(countPDict(pdict, d, max.mismatch=2), c(2,0,2,2))
})

test_that("vectorized PDict lookup works", {
  ## canary tests for future changes
  pdict <- PDict(c("acgt", "gt", "cgt", "ac"), tb.end=2)
  d <- DNAString("acggaccg")
  dss <- DNAStringSet(list(d,d))
  dvv <- Views(d, start=rep(1,2), width=rep(length(d), 2))

  # Disallowed cases
  expect_error(matchPDict(pdict, dss), "please use vmatchPDict")
  expect_error(vmatchPDict(pdict, d), "please use matchPDict")
  expect_error(vmatchPDict(pdict, as(d, "MaskedDNAString")), "please use matchPDict")
  expect_error(vwhichPDict(pdict, d), "please use whichPDict")
  expect_error(vwhichPDict(pdict, as(d, "MaskedDNAString")), "please use whichPDict")
  expect_error(vcountPDict(pdict, d), "please use countPDict")
  expect_error(vcountPDict(pdict, as(d, "MaskedDNAString")), "please use countPDict")

  exp_out <- matrix(rep(c(0,2), times=c(6,2)), nrow=4L, byrow=TRUE)
  expect_equal(vcountPDict(pdict, dss), exp_out)
  expect_equal(vwhichPDict(pdict, dss), list(4,4))
  expect_equal(vcountPDict(pdict, dvv), exp_out)
  expect_equal(vwhichPDict(pdict, dvv), list(4,4))

  exp_out[1,] <- 1
  exp_out[3,] <- 2
  expect_equal(vcountPDict(pdict, dss, max.mismatch=1), exp_out)
  expect_equal(vwhichPDict(pdict, dss, max.mismatch=1), list(c(1,3,4),c(1,3,4)))
  expect_equal(vcountPDict(pdict, dvv, max.mismatch=1), exp_out)
  expect_equal(vwhichPDict(pdict, dvv, max.mismatch=1), list(c(1,3,4),c(1,3,4)))

  exp_out[1,] <- 2
  expect_equal(vcountPDict(pdict, dss, max.mismatch=2), exp_out)
  expect_equal(vwhichPDict(pdict, dss, max.mismatch=2), list(c(1,3,4),c(1,3,4)))
  expect_equal(vcountPDict(pdict, dvv, max.mismatch=2), exp_out)
  expect_equal(vwhichPDict(pdict, dvv, max.mismatch=2), list(c(1,3,4),c(1,3,4)))

  ## this will fail when we implement it to remind me to add tests for it
  expect_error(vmatchPDict(pdict, dss), "vmatchPDict() is not ready yet", fixed=TRUE)
})

test_that("MIndex functionality works", {
  ## set up an MIndex object
  s <- mkAllStrings(c("A", "T"), 2L)
  names(s) <- s
  pd <- PDict(s)
  d <- DNAString("ATTATTAATTAA")
  m <- matchPDict(pd, d)
  exp_match_starts <- list(c(7,11), c(1,4,8), c(3,6,10), c(2,5,9))
  expect_equal(length(m), 2^2)
  expect_equal(names(m), names(s))
  expect_equal(startIndex(m), exp_match_starts)
  expect_equal(endIndex(m), lapply(exp_match_starts, \(x) x+1L))
  expect_equal(elementNROWS(m), c(2,3,3,3))
  expect_equal(m[[1L]], IRanges(start=c(7,11), end=c(8,12)))
  expect_s4_class(as(m, "CompressedIRangesList"), "CompressedIRangesList")
  expect_equal(start(unlist(m)), unlist(exp_match_starts))
  ux <- extractAllMatches(d, m)
  exp_ux <- rep(s, times=c(2,3,3,3))
  expect_equal(start(ux), unlist(exp_match_starts))
  expect_equal(as.character(ux), exp_ux)
})

### All old tests
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
