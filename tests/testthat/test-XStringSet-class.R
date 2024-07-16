dnastr <- paste(DNA_ALPHABET, collapse='')
rnastr <- paste(RNA_ALPHABET, collapse='')
aastr <- paste(AA_ALPHABET, collapse='')
bstr <- rawToChar(as.raw(32:126))

d <- DNAStringSet(dnastr)
r <- RNAStringSet(rnastr)
a <- AAStringSet(aastr)
b <- BStringSet(bstr)

# d <- DNAString(dnastr)
# r <- RNAString(rnastr)
# a <- AAString(aastr)
# b <- BString(bstr)

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

test_that("concatenation and character conversion are correct", {
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
})