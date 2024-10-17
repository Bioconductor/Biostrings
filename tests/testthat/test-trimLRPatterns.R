## trimLRPatterns.R exports trimLRPatterns for the following classes:
## - character
## - XString
## - XStringSet

test_that("trimLRPatterns works correctly on supported input types", {
	## tests based on the examples in the docs
	L <- "TTCTGCTTG"
  R <- "GATCGGAAG"
  DChr <- "TTCTGCTTGACGTGATCGGA"
  DStr <- DNAString(DChr)
  DSet <- DNAStringSet(c(DChr, DChr))

  exp1 <- "ACGTGATCGGA"
  exp2 <- "TTCTGCTTGACGT"
  exp3 <- "ACGT"

  DSet2 <- DNAStringSet(c("TGCTTGACGGCAGATCGG", "TTCTGCTTGGATCGGAAG"))

  ## Perfect matches on the flanks
  ### left match
  expect_true(trimLRPatterns(Lpattern = L, subject = DChr) == exp1)
  expect_true(trimLRPatterns(Lpattern = L, subject = DStr) == DNAString(exp1))
  expect_true(all(trimLRPatterns(Lpattern = L, subject = DSet) == DNAStringSet(rep(exp1, 2L))))

  ### right match
  expect_true(trimLRPatterns(Rpattern = R, subject = DChr) == exp2)
  expect_true(trimLRPatterns(Rpattern = R, subject = DStr) == DNAString(exp2))
  expect_true(all(trimLRPatterns(Rpattern = R, subject = DSet) == DNAStringSet(rep(exp2, 2L))))

  ### both match, perfect matches on flanking overlaps
  expect_true(trimLRPatterns(Lpattern=L, Rpattern = R, subject = DChr) == exp3)
  expect_true(trimLRPatterns(Lpattern=L, Rpattern = R, subject = DStr) == DNAString(exp3))
  expect_true(all(trimLRPatterns(Lpattern=L, Rpattern = R, subject = DSet) == DNAStringSet(rep(exp3, 2L))))

  ## Mismatches on the flanks

  ## Allow for mismatches on the flanks
  t1 <- trimLRPatterns(Lpattern = L, Rpattern = R, subject = DSet2,
  								max.Lmismatch = 0.2, max.Rmismatch = 0.2)
  expect_true(all(t1 == DNAStringSet(c("ACGGCA", ""))))

  maxMismatches <- as.integer(0.2 * 1:9)
  maxMismatches
  t2 <- trimLRPatterns(Lpattern = L, Rpattern = R, subject = DSet2,
        					max.Lmismatch = maxMismatches, max.Rmismatch = maxMismatches)
  ## This should have a better example
  expect_true(all(t1 == t2))

  t3 <- trimLRPatterns(Lpattern = L, Rpattern = R, subject = DSet2[1],
                 max.Lmismatch = maxMismatches, max.Rmismatch = maxMismatches,
                 with.Rindels=TRUE)
  expect_true(t3 == DNAStringSet("ACGGC"))
  t4 <- trimLRPatterns(Lpattern = L, Rpattern = R, subject = DSet2[1],
                 max.Lmismatch = maxMismatches, max.Rmismatch = maxMismatches,
                 with.Rindels=TRUE, with.Lindels=TRUE)
  expect_true(t4 == DNAStringSet("CGGC"))


  ## Produce ranges that can be an input into other functions
  expect_s4_class(trimLRPatterns(Lpattern=L, Rpattern=R, subject=DSet2, ranges=TRUE), "IRanges")
})