## Miscellaneous XString Methods
## Currently includes:
## - toComplex
## - xscat
## - strsplit, unstrsplit (for XStringSet, XStringSetList)
## - N50

test_that("toComplex conversion works correctly", {
	seq <- DNAString("agtcagtcagtcagtc")
  baseValues1 <- c(A=1+0i, G=0+1i, T=-1+0i, C=0-1i)
  expect_equal(toComplex(seq, baseValues1), unname(rep(baseValues1, 4)))

  ## GC content:
  baseValues2 <- c(A=0, C=1, G=1, T=0)
  expect_equal(sum(as.integer(toComplex(seq, baseValues2))), 8L)

  expect_error(toComplex(seq, 1:4), "'baseValues' must have names")
  expect_error(toComplex(seq, c(W=1,X=1,Y=1,Z=1)), "'baseValues' names must be valid DNA letters")
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

test_that("Biostrings strsplit, unstrsplit methods work correctly", {
  ds <- "ATGCATGC"
  d1 <- DNAStringSet(ds)
  d2 <- DNAStringSet("ATCATC")
  d_spl <- DNAStringSet(c("AT", "CAT", "C"))

  expect_true(all(strsplit(d1, "G") == d_spl))

  expect_true(unstrsplit(strsplit(d1, "G")) == d2)

  ## round trip identity
  expect_true(unstrsplit(strsplit(d1, "G"), "G") == d1)

  ## unstrsplit on a stringset should be a noop
  expect_identical(unstrsplit(d1), d1)
})

test_that("N50 has correct behavior", {
  set.seed(6L)
  n_seq <- 20L
  possible_len <- seq(100,10000)
  lens <- sample(possible_len, n_seq)
  many_seqs <- DNAStringSet(
    vapply(seq_len(n_seq), \(i){
      paste(sample(DNA_BASES, lens[i], replace=TRUE), collapse='')
    }, character(1L))
  )

  ## N50 is calculated by:
  ##  1. add sizes of contigs until half the total size is reached
  all_widths <- sort(width(many_seqs), decreasing=TRUE)
  total_width <- sum(all_widths)
  cs <- cumsum(all_widths)
  pos_first <- which.max(cs > total_width/2)

  ##  2. take the size of the contig that was added last
  N50_exp <- all_widths[pos_first]

  expect_equal(N50(width(many_seqs)), N50_exp)
})