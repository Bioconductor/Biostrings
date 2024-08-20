## findPalindromes.R exports the following:
## - findPalindromes()
## - palindromeArmLength()
## - palindromeLeftArm()
## - palindromeRightArm()
##
## TODOs: missing MaskedXString, XStringViews tests


## new tests
test_that("palindromeLeftArm and *RightArm identify the right arms", {
  text4 <- BString("i45hgfe7d321c3b4a56789uvwWVU98765A4B3C123D7EFGH54I")
  current <- findPalindromes(text4, min.armlength = 5, max.looplength = length(text4), max.mismatch = 3)
  leftarm <- palindromeLeftArm(current)
  expect_equal(start(leftarm), start(current))
  expect_equal(width(leftarm), c(5,1,1,3,1,0))

  rightarm <- palindromeRightArm(current)
  expect_equal(end(rightarm), end(current))
  expect_equal(width(rightarm), width(leftarm))
})

## For now I'm just going to port the old tests
test_that("findPalindromes works on BString objects", {
    text1 <- BString("ABCDCBA")
    current <- findPalindromes(text1, min.armlength=2)
    expect_identical(IRanges(1, 7), ranges(current))
    current <- findPalindromes(text1, min.armlength=2, max.looplength=0)
    expect_identical(IRanges(), ranges(current))

    text2 <- BString("xesopeReste etseReposeygohangasalamiImalasagnahogz")
    current <- findPalindromes(text2, max.looplength=2)
    expect_identical(IRanges(c(2, 24), c(22, 49)), ranges(current))
    expect_identical(c(10L, 12L), palindromeArmLength(current))

    text3 <- BString("AATAAACTNTCAAATYCCY")
    current <- findPalindromes(text3, min.armlength=2)
    expect_identical(IRanges(c(1, 3, 16), c(5, 15, 19)), ranges(current))

    text4 <- BString("i45hgfe7d321c3b4a56789uvwWVU98765A4B3C123D7EFGH54I")
    current <- findPalindromes(text4, min.armlength = 5, max.looplength = length(text4), max.mismatch = 3)
    expect_identical(IRanges(c(18, 16, 14, 10, 8, 7), c(33, 35, 37, 41, 43, 44)), ranges(current))
})

test_that("findPalindromes works on nucleotide sequences", {
  x1 <- DNAString("AATAAACTNTCAAATYCCY")  # same sequence as 'text3'
  for (class in c("DNAString", "RNAString")) {
      x <- as(x1, class)
      current <- findPalindromes(x, min.armlength=2)
      expect_identical(IRanges(), ranges(current))
      current <- findPalindromes(x, min.armlength=2, max.looplength=10)
      expect_identical(IRanges(2, 15), ranges(current))
  }

  x2 <- DNAString("AATTT")
  for (class in c("DNAString", "RNAString")) {
      x <- as(x2, class)
      current <- findPalindromes(x, min.armlength=2)
      target <- IRanges(1, 4:5)
      expect_identical(target, ranges(current))
      current <- findPalindromes(x, min.armlength=2, max.looplength=0)
      expect_identical(target[1], ranges(current))
  }

  x3 <- DNAString("TTTAA")
  for (class in c("DNAString", "RNAString")) {
      x <- as(x3, class)
      current <- findPalindromes(x, min.armlength=2)
      target <- IRanges(1:2, 5)
      expect_identical(target, ranges(current))
      current <- findPalindromes(x, min.armlength=2, max.looplength=0)
      expect_identical(target[2], ranges(current))
  }

  x4 <- DNAString("AAATTT")
  for (class in c("DNAString", "RNAString")) {
      x <- as(x4, class)
      current <- findPalindromes(x, min.armlength=2)
      target <- IRanges(c(1, 1, 2), c(5, 6, 6))
      expect_identical(target, ranges(current))
      current <- findPalindromes(x, min.armlength=3)
      expect_identical(target[2], ranges(current))
      current <- findPalindromes(x, min.armlength=4)
      expect_identical(target[0], ranges(current))
  }

  ## With nested palindromes.
  x5 <- DNAString("CTTGAAATTTGAAG")
  for (class in c("DNAString", "RNAString")) {
      x <- as(x5, class)
      target0 <- IRanges(c(2, 2, 5,  5,  1,  6,  8,  9),
                         c(6, 7, 9, 10, 14, 10, 13, 13))
      current <- findPalindromes(x, min.armlength=2, max.looplength=0)
      expect_identical(target0[4], ranges(current))
      current <- findPalindromes(x, min.armlength=2, max.looplength=1)
      expect_identical(target0[-c(2, 5, 7)], ranges(current))
      for (i in 2:7) {
          current <- findPalindromes(x, min.armlength=2, max.looplength=i)
          expect_identical(target0[-5], ranges(current))
      }
      current <- findPalindromes(x, min.armlength=2, max.looplength=8)
      expect_identical(target0, ranges(current))
  }

  x6 <- DNAString("AAGAATTGTT")
  for (class in c("DNAString", "RNAString")) {
      x <- as(x6, class)
      target0 <- IRanges(c(1, 4, 1, 4), c(7, 7, 10, 10))
      for (i in 0:2) {
          current <- findPalindromes(x, min.armlength=2, max.looplength=i)
          expect_identical(target0[2], ranges(current))
      }
      for (i in 3:5) {
          current <- findPalindromes(x, min.armlength=2, max.looplength=i)
          expect_identical(target0[-3], ranges(current))
      }
      current <- findPalindromes(x, min.armlength=2, max.looplength=6)
      expect_identical(target0, ranges(current))
  }

  # With IUPAC ambiguity codes. Ambiguity codes are not allowed in the arms
  # of a palindrome!
  x7 <- DNAString("NNNNNN")
  for (class in c("DNAString", "RNAString")) {
      x <- as(x7, class)
      current <- findPalindromes(x, min.armlength=2)
      expect_identical(IRanges(), ranges(current))
  }

  x8 <- DNAString("ACCGNAAATTTNCGGT")
  for (class in c("DNAString", "RNAString")) {
      x <- as(x8, class)
      current <- findPalindromes(x, min.armlength=2)
      expect_identical(IRanges(c(6, 6, 7), c(10, 11, 11)), ranges(current))
      current <- findPalindromes(x, min.armlength=4)
      expect_identical(IRanges(), ranges(current))
      current <- findPalindromes(x, min.armlength=3, max.looplength=8)
      expect_identical(IRanges(c(6, 1), c(11, 16)), ranges(current))
  }
})

test_that("palindromeArmLength functions correctly for nucleotide sequences", {
    x2 <- DNAString("AATTT")
    for (class in c("DNAString", "RNAString")) {
        x <- as(x2, class)
        expect_identical(2L, palindromeArmLength(x))
        pals <- findPalindromes(x, min.armlength=2)
        expect_identical(c(2L, 2L), palindromeArmLength(pals))
    }

    x3 <- DNAString("TTTAA")
    for (class in c("DNAString", "RNAString")) {
        x <- as(x3, class)
        expect_identical(2L, palindromeArmLength(x))
        pals <- findPalindromes(x, min.armlength=2)
        expect_identical(c(2L, 2L), palindromeArmLength(pals))
    }

    x4 <- DNAString("AAATTT")
    for (class in c("DNAString", "RNAString")) {
        x <- as(x4, class)
        expect_identical(3L, palindromeArmLength(x))
        pals <- findPalindromes(x, min.armlength=2)
        expect_identical(c(2L, 3L, 2L), palindromeArmLength(pals))
    }

    x5 <- DNAString("CTTGAAATTTGAAG")
    for (class in c("DNAString", "RNAString")) {
        x <- as(x5, class)
        expect_identical(3L, palindromeArmLength(x))
        pals <- findPalindromes(x, min.armlength=2, max.looplength=8)
        expect_identical(c(2L, 2L, 2L, 3L, 3L, 2L, 2L, 2L),
                       palindromeArmLength(pals))
    }
})
