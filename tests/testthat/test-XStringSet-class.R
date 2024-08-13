dnastr <- paste(DNA_ALPHABET, collapse='')
rnastr <- paste(RNA_ALPHABET, collapse='')
aastr <- paste(AA_ALPHABET, collapse='')
bstr <- rawToChar(as.raw(32:126))

d <- DNAStringSet(dnastr)
r <- RNAStringSet(rnastr)
a <- AAStringSet(aastr)
b <- BStringSet(bstr)

## testing functions from old testing files
### '.eltAddresses(x)' collects the addresses of the elements in 'x' (in
### practice 'x' will be a list of external pointers or environments).
.eltAddresses <- function(x) sapply(x, XVector:::address)

### 'x' and 'y' must be XVectorList vectors.
.haveIdenticalPools <- function(x, y)
    identical(.eltAddresses(x@pool@xp_list), .eltAddresses(y@pool@xp_list))

### 'x' must be an XVectorList vector.
.poolEltLengths <- function(x)
{
    pool_len <- length(x@pool)
    if (pool_len == 0L)
        return(integer(0))
    sapply(seq_len(pool_len), function(i) length(x@pool[[i]]))
}

test_that("seqtype correctly infers types", {
  expect_equal(seqtype(d), "DNA")
  expect_equal(seqtype(r), "RNA")
  expect_equal(seqtype(a), "AA")
  expect_equal(seqtype(b), "B")

  expect_equal(seqtype(d[[1]]), "DNA")
  expect_equal(seqtype(r[[1]]), "RNA")
  expect_equal(seqtype(a[[1]]), "AA")
  expect_equal(seqtype(b[[1]]), "B")
})

test_that("seqtype conversion works as expected", {
  ## conversion between RNA and DNA
  expect_equal({x <- d; seqtype(x) <- "RNA"; x}, r)
  expect_equal({x <- r; seqtype(x) <- "DNA"; x}, d)

  ## conversion to BStringSets
  expect_equal({x <- d; seqtype(x) <- "B"; x}, BStringSet(dnastr))
  expect_equal({x <- r; seqtype(x) <- "B"; x}, BStringSet(rnastr))
  expect_equal({x <- a; seqtype(x) <- "B"; x}, BStringSet(aastr))

  expect_equal({x <- BStringSet(dnastr); seqtype(x) <- "DNA"; x}, d)
  expect_equal({x <- BStringSet(rnastr); seqtype(x) <- "RNA"; x}, r)
  expect_equal({x <- BStringSet(aastr); seqtype(x) <- "AA"; x}, a)

  ## invalid conversions
  expect_error(seqtype(d) <- "AA", "incompatible sequence types")
  expect_error(seqtype(r) <- "AA", "incompatible sequence types")
  expect_error(seqtype(a) <- "DNA", "incompatible sequence types")
  expect_error(seqtype(a) <- "RNA", "incompatible sequence types")
  expect_error(seqtype(b) <- "AA", "not in lookup table")
  expect_error(seqtype(b) <- "DNA", "not in lookup table")
  expect_error(seqtype(b) <- "RNA", "not in lookup table")
})

test_that("unlisting works", {
  expect_equal(unlist(d), d[[1]])
  expect_equal(unlist(r), r[[1]])
  expect_equal(unlist(a), a[[1]])
  expect_equal(unlist(b), b[[1]])
})

test_that("width and nchar are correct", {
  expect_equal(width(d), nchar(dnastr))
  expect_equal(width(r), nchar(rnastr))
  expect_equal(width(a), nchar(aastr))
  expect_equal(width(b), nchar(bstr))

  expect_equal(nchar(d), nchar(dnastr))
  expect_equal(nchar(r), nchar(rnastr))
  expect_equal(nchar(a), nchar(aastr))
  expect_equal(nchar(b), nchar(bstr))

  expect_error(width(NA_character_), "NAs in 'x' are not supported")
})

test_that("concatenation and character/vector conversion are correct", {
  expect_s4_class(c(d,d), "DNAStringSet")
  expect_s4_class(c(r,r), "RNAStringSet")
  expect_s4_class(c(a,a), "AAStringSet")
  expect_s4_class(c(b,b), "BStringSet")

  dd <- c(d,d)
  rr <- c(r,r)
  aa <- c(a,a)
  bb <- c(b,b)

  expect_equal(dd == DNAStringSet(c(dnastr, dnastr)), c(TRUE, TRUE))
  expect_equal(rr == RNAStringSet(c(rnastr, rnastr)), c(TRUE, TRUE))
  expect_equal(aa == AAStringSet(c(aastr, aastr)), c(TRUE, TRUE))
  expect_equal(bb == BStringSet(c(bstr, bstr)), c(TRUE, TRUE))

  expect_equal(as.character(dd), c(dnastr, dnastr))
  expect_equal(as.character(rr), c(rnastr, rnastr))
  expect_equal(as.character(aa), c(aastr, aastr))
  expect_equal(as.character(bb), c(bstr, bstr))

  expect_equal(as.vector(dd), c(dnastr, dnastr))
  expect_equal(as.vector(rr), c(rnastr, rnastr))
  expect_equal(as.vector(aa), c(aastr, aastr))
  expect_equal(as.vector(bb), c(bstr, bstr))
})

test_that("matrix conversion is correct", {
  dd <- c(d,d)
  rr <- c(r,r)
  aa <- c(a,a)
  bb <- c(b,b)

  md <- do.call(rbind, rep(strsplit(dnastr, ''), 2))
  mr <- do.call(rbind, rep(strsplit(rnastr, ''), 2))
  ma <- do.call(rbind, rep(strsplit(aastr, ''), 2))
  mb <- do.call(rbind, rep(strsplit(bstr, ''), 2))


  expect_equal(as.matrix(dd), md)
  expect_equal(as.matrix(rr), mr)
  expect_equal(as.matrix(aa), ma)
  expect_equal(as.matrix(bb), mb)

  md <- as.matrix(DNAStringSet(''))
  mr <- as.matrix(RNAStringSet(''))
  ma <- as.matrix(AAStringSet(''))
  mb <- as.matrix(BStringSet(''))

  m_base <- matrix('', nrow=1, ncol=0)
  expect_equal(md, m_base)
  expect_equal(mr, m_base)
  expect_equal(ma, m_base)
  expect_equal(mb, m_base)
})

test_that("factor conversion is correct", {
  expect_equal(as.factor(d), as.factor(dnastr))
  expect_equal(as.factor(r), as.factor(rnastr))
  expect_equal(as.factor(a), as.factor(aastr))
  expect_equal(as.factor(b), as.factor(bstr))
})

test_that("data.frame conversion is correct", {
  expect_equal(as.data.frame(c(d,d)), data.frame(x=c(dnastr, dnastr)))
  expect_equal(as.data.frame(c(r,r)), data.frame(x=c(rnastr, rnastr)))
  expect_equal(as.data.frame(c(a,a)), data.frame(x=c(aastr, aastr)))
  expect_equal(as.data.frame(c(b,b)), data.frame(x=c(bstr, bstr)))
})

test_that("toString conversion is correct", {
  expect_equal(toString(c(d,d)), paste(dnastr, dnastr, sep=', '))
  expect_equal(toString(c(r,r)), paste(rnastr, rnastr, sep=', '))
  expect_equal(toString(c(a,a)), paste(aastr, aastr, sep=', '))
  expect_equal(toString(c(b,b)), paste(bstr, bstr, sep=', '))
})

test_that("show methods output correctly", {
  expect_output(show(c(d,d)), "^DNAStringSet object of length 2:\\n")
  expect_output(show(c(r,r)), "^RNAStringSet object of length 2:\\n")
  expect_output(show(c(a,a)), "^AAStringSet object of length 2:\\n")
  expect_output(show(c(b,b)), "^BStringSet object of length 2:\\n")
})

test_that("showAsCell works correctly", {
  d <- DNAStringSet(paste(rep(dnastr, 10), collapse=''))
  r <- RNAStringSet(paste(rep(rnastr, 10), collapse=''))
  a <- AAStringSet(paste(rep(aastr, 10), collapse=''))

  expect_equal(nchar(showAsCell(d)), 23L)
  expect_equal(nchar(showAsCell(r)), 23L)
  expect_equal(nchar(showAsCell(a)), 23L)
  expect_equal(nchar(showAsCell(b)), 23L)
})

test_that("comparison operators are correct", {
  dna <- DNAStringSet(DNA_ALPHABET)
  rna <- RNAStringSet(RNA_ALPHABET)
  aaa <- AAStringSet(AA_ALPHABET)
  bbb <- BStringSet(LETTERS)

  expect_true(!any(is.na(dna)))
  expect_true(!anyNA(dna))
  expect_equal(match(dna, dna), seq_along(dna))
  expect_equal(aaa[seq_len(26)] < bbb, AA_ALPHABET[seq_len(26L)] < LETTERS)

  expect_equal(match(sort(aaa), bbb, nomatch=0), c(rep(0L, 4L), seq_len(26L)))
  expect_true(all(dna == as.character(dna)))

  expect_error(aaa == dna, "is not supported")
  expect_error(aaa == rna, "is not supported")
  expect_true(all(dna == BStringSet(DNA_ALPHABET)))
  expect_true(all(rna == BStringSet(RNA_ALPHABET)))
  expect_true(all(aaa == BStringSet(AA_ALPHABET)))
  expect_equal(dna == NULL, logical(0L))
})

## Porting RUnit tests
test_that("short RUnit tests continue to pass", {
  ## test_width_character
  x <- safeExplode(rawToChar(as.raw(1:255)))
  expect_equal(width(x), rep.int(1L, 255))

  ## DNAStringSet internal elements ##
  dna <- DNAStringSet(DNA_ALPHABET)

  ## DNAStringSet_constructor
  expect_equal(.poolEltLengths(dna), length(DNA_ALPHABET))

  ## DNAStringSet_width
  expect_equal(width(dna), width(DNA_ALPHABET))

  ## DNAStringSet_unlist
  expect_equal(as.character(unlist(dna)), dnastr)

  ## DNAStringSet_showAsCell
  expect_equal(showAsCell(DNAStringSet()), character(0L))
  expect_equal(showAsCell(dna), DNA_ALPHABET)
})

test_that("RUnit test_DNAStringSet_subsetting", {
  dna <- DNAStringSet(DNA_ALPHABET)
  elementMetadata(dna) <- DataFrame(C1=dna)

  dna0 <- dna[FALSE]
  expect_equal(length(dna0), 0L)
  ## Checking internal representation.
  expect_equal(.poolEltLengths(dna0), integer(0))
  expect_true(.haveIdenticalPools(elementMetadata(dna0)$C1, dna0))
  expect_equal(elementMetadata(dna0)$C1@ranges, dna0@ranges)

  idx <- rep.int((8:6)*2L, 100L)
  dna300 <- dna[idx]
  expect_equal(length(dna300), length(idx))
  ## Checking internal representation.
  expect_true(.haveIdenticalPools(dna300, dna))
  expect_true(.haveIdenticalPools(elementMetadata(dna300)$C1, dna300))
  expect_equal(elementMetadata(dna300)$C1@ranges, dna300@ranges)
})

test_that("RUnit test_DNAStringSet_combining", {
  dna <- DNAStringSet(DNA_ALPHABET)
  elementMetadata(dna) <- DataFrame(C1=dna)

  dna2a <- c(dna, dna)
  dna2b <- rep(dna, 2L)
  expect_equal(dna2a, dna2b)

  ## Checking internal representation.
  expect_true(.haveIdenticalPools(dna2a, dna))
  expect_true(.haveIdenticalPools(dna2a, dna2b))
  expect_true(.haveIdenticalPools(elementMetadata(dna2a)$C1, dna2a))
  expect_equal(elementMetadata(dna2a)$C1@ranges, dna2a@ranges)
})

test_that("RUnit test_DNAStringSet_compaction", {
  dna <- DNAStringSet(DNA_ALPHABET)
  elementMetadata(dna) <- DataFrame(C1=dna)

  idx <- rep.int((8:6)*2L, 100L)
  dna300 <- dna[idx]
  compact_dna300 <- compact(dna300)
  expect_equal(as.character(compact_dna300), as.character(dna300))
  ## Checking internal representation.
  expect_equal(.poolEltLengths(compact_dna300), 3L)
  expect_equal(.poolEltLengths(elementMetadata(compact_dna300)$C1),
                .poolEltLengths(compact_dna300))
  expect_equal(elementMetadata(compact_dna300)$C1@ranges,
                compact_dna300@ranges)
})