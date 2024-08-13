## Miscellaneous XString Methods
## Currently includes:
## - toComplex
## - xscat

test_that("toComplex conversion works correctly", {
	seq <- DNAString("agtcagtcagtcagtc")
  baseValues1 <- c(A=1+0i, G=0+1i, T=-1+0i, C=0-1i)
  expect_equal(toComplex(seq, baseValues1), unname(rep(baseValues1, 4)))

  ## GC content:
  baseValues2 <- c(A=0, C=1, G=1, T=0)
  expect_equal(sum(as.integer(toComplex(seq, baseValues2))), 8L)
})

test_that("xscat concatenates all types correctly", {
  s1 <- "ATGATG"
  s2 <- "CCACCA"

  ## character
  expect_true(xscat(s1, s2) == BString(paste0(s1, s2)))
  expect_error(xscat(s1, NA), "arguments must be character vectors (with no NAs)", fixed=TRUE)

  ## XString
  d1 <- DNAString(s1)
  d2 <- DNAString(s2)
  d3 <- DNAString(paste0(s1, s2))
  d4 <- DNAString(paste0(s2, s1))
  expect_true(xscat(d1, d2) == d3)

  ## XStringSet
  dss1 <- DNAStringSet(list(d1, d2))
  dss2 <- DNAStringSet(list(d2, d1))
  dss3 <- DNAStringSet(list(d3, d4))
  expect_true(all(xscat(dss1, dss2) == dss3))

  v1 <- Views(d3, start=c(1,7), width=3)
  v2 <- Views(d3, start=c(4,10), width=3)
  expect_true(all(xscat(v1, v2) == DNAStringSet(c(s1, s2))))
})

test_that("padAndClip, stackStrings have correct functionality", {
  ## just replicating the stuff in the examples
  ## can add more tests later if necessary, I don't see these used super often
  x <- BStringSet(c(seq1="ABCD", seq2="abcdefghijk", seq3="", seq4="XYZ"))

  p1 <- padAndClip(x, IRanges(3, 8:5), Lpadding.letter=">", Rpadding.letter="<")
  expect_equal(unname(as.character(p1)), c("CD<<<<", "cdefg", "<<<<", "Z<<"))

  p2 <- padAndClip(x, IRanges(1:-2, 7), Lpadding.letter=">", Rpadding.letter="<")
  expect_equal(unname(as.character(p2)), c("ABCD<<<", ">abcdefg", ">><<<<<<<", ">>>XYZ<<<<"))

  expect_equal(width(stackStrings(x, 2, 8)), rep(7L, 4L))

  exp_out <- c("###ABCD....", "ijk........", "#########..", "##########X")
  s1 <- stackStrings(x, -2, 8, shift=c(0, -11, 6, 7), Lpadding.letter="#", Rpadding.letter=".")
  expect_equal(width(s1), rep(11, 4L))
  expect_equal(unname(as.character(s1)), exp_out)

  s2 <- stackStrings(x, -2, 8, shift=c(0, -14, 6, 7),
               Lpadding.letter="#", Rpadding.letter=".",
               remove.out.of.view.strings=TRUE)
  expect_equal(names(s2), paste0("seq", c(1,3,4)))
  expect_equal(unname(as.character(s2)), exp_out[c(1,3,4)])
})