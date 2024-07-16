dnastr <- paste(DNA_ALPHABET, collapse='')
rnastr <- paste(RNA_ALPHABET, collapse='')
aastr <- paste(AA_ALPHABET, collapse='')
bstr <- rawToChar(as.raw(32:126))

# d <- DNAString(dnastr)
# r <- RNAString(rnastr)
# a <- AAString(aastr)
# b <- BString(bstr)

test_that("seqtype correctly infers types", {
  expect_equal(seqtype(DNAStringSet(c(dnastr, dnastr))), "DNA")
  expect_equal(seqtype(RNAStringSet(c(rnastr, rnastr))), "RNA")
  expect_equal(seqtype(AAStringSet(c(aastr, aastr))), "AA")
  expect_equal(seqtype(BStringSet(c(bstr, bstr))), "B")
})