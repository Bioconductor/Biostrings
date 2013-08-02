test_twoWayAlphabetFrequency <- function() {
  ## Equivalence classes:
  ##
  ## x, y:
  ## - length [Set]: (different length), zero, length one, length > one
  ## - string length: (inconsistent lengths), zero, same, [Set] varying
  ## - content: some matches, every cell
  ## - types: XString, DNAString, XStringSet, DNAStringSet, (different)
  ##
  ## as.prob: TRUE, FALSE, (NA), (not logical)
  ## baseOnly [DNA]: TRUE, FALSE, (NA), (not logical)
  ## collapse [Set]: TRUE, FALSE, (NA), (not logical)

  library(Biostrings)
  library(RUnit)
  
  testMatrix <-
    function(x, y, v, as.prob = FALSE, baseOnly = FALSE, collapse = FALSE) {
      x.alphabet <- Biostrings:::xscodes(x, baseOnly=baseOnly)
      y.alphabet <- Biostrings:::xscodes(y, baseOnly=baseOnly)
      dim <- c(length(x.alphabet), length(y.alphabet))
      dimnames <- list(names(x.alphabet), names(y.alphabet))
      if (baseOnly) {
        if (!is.null(dimnames[[1]]))
          dimnames[[1]] <- c(dimnames[[1]], "other")
        if (!is.null(dimnames[[2]]))
          dimnames[[2]] <- c(dimnames[[2]], "other")
        dim <- dim + 1L
      }
      mode(v) <- "integer"
      if (as.prob) {
        if (collapse)
          v <- v / sum(v)
        else v <- v / rep(nchar(x), each = prod(dim(v)))
      }
      if (is(x, "XStringSet") || is(y, "XStringSet")) {
        dimnames <- c(dimnames, list(NULL))
        dim <- c(dim, nchar(x))
      }
      a <- array(0L, dim = dim, dimnames = dimnames)
      if (is.null(rownames(a))) {
        x.ind <- as.integer(sapply(rownames(v), charToRaw)) + 1L
      } else {
        x.ind <- rownames(v)
      }
      if (is.null(colnames(a))) {
        y.ind <- as.integer(sapply(rownames(v), charToRaw)) + 1L
      } else {
        y.ind <- rownames(v)
      }
      if (length(dim(a)) == 3)
        a[x.ind,y.ind,] <- v
      else a[x.ind,y.ind] <- v
      a
    }
  
  some.a <- "ACGTACGTACGT"
  some.b <- "ACGTTGCAGGGA"
  some.mat <- rbind(A = c(1, 0, 1, 1),
                    C = c(0, 1, 2, 0),
                    G = c(0, 1, 2, 0),
                    T = c(2, 0, 0, 1))
  
  ## CASE: different nchar, DNAString
  checkException(twoWayAlphabetFrequency(DNAString("ACGT"), DNAString("ACG")),
                 silent = TRUE)
  
  ## CASE: zero, DNAString
  x <- y <- DNAString("")
  checkIdentical(testMatrix(x, y, 0L), twoWayAlphabetFrequency(x, y))
  
  ## CASE: some, XString, FALSE
  x <- BString(some.a)
  y <- BString(some.b)
  truth <- testMatrix(x, y, some.mat)
  checkIdentical(truth, twoWayAlphabetFrequency(x, y))
  
  ## CASE: some, XString, TRUE
  truth <- truth / sum(truth)
  checkIdentical(truth, twoWayAlphabetFrequency(x, y, as.prob = TRUE))
  
  ## CASE: XString, DNAString
  x <- BString(some.a)
  y <- DNAString(some.b)
  
  truth <- testMatrix(x, y, some.mat)
  checkIdentical(truth, twoWayAlphabetFrequency(x, y))

  ## CASE: XString, DNAString, baseOnly = TRUE
  truth <- testMatrix(x, y, some.mat, baseOnly = TRUE)
  checkIdentical(truth, twoWayAlphabetFrequency(x, y, baseOnly = TRUE))
  
  ## CASE: some, DNAString, NA, FALSE
  x <- DNAString(some.a)
  y <- DNAString(some.b)
  
  checkException(twoWayAlphabetFrequency(x, y, as.prob = NA),
                 silent = TRUE)  
  
  ## CASE: some, DNAString, "foo", FALSE
  checkException(twoWayAlphabetFrequency(x, y, as.prob = "foo"),
                 silent = TRUE)

  
  ## CASE: some, DNAString, FALSE, NA
  checkException(twoWayAlphabetFrequency(x, y, baseOnly = NA),
                 silent = TRUE)

  
  ## CASE: some, DNAString, FALSE, "foo"
  checkException(twoWayAlphabetFrequency(x, y, baseOnly = "foo"),
                 silent = TRUE)

  ## CASE: some, DNAString, TRUE, FALSE
  truth <- testMatrix(x, y, some.mat, as.prob = TRUE)
  checkIdentical(truth, twoWayAlphabetFrequency(x, y, as.prob = TRUE))
  
  ## CASE: some, DNAString, FALSE, TRUE
  truth <- testMatrix(x, y, some.mat, baseOnly = TRUE)
  checkIdentical(truth, twoWayAlphabetFrequency(x, y, baseOnly = TRUE))
  
  ## CASE: some, DNAString, TRUE, TRUE
  ## CASE: every, DNAString, FALSE, FALSE
  ## CASE: DNAStringSet, DNAString
  ## CASE: different length, DNAStringSet
  ## CASE: different length, XStringSet
  ## CASE: different widths, DNAStringSet
  ## CASE: different widths, XStringSet
  ## CASE: length==0, DNAString, FALSE, FALSE, FALSE
  ## CASE: length==0, DNAString, FALSE, FALSE, TRUE
  ## CASE: length==1, some, DNAStringSet, FALSE, FALSE, FALSE
  ## CASE: length==1, some, DNAStringSet, FALSE, FALSE, TRUE
  ## CASE: length >1, some, XStringSet, FALSE, , FALSE
  ## CASE: length >1, some, XStringSet, TRUE, , FALSE
  ## CASE: length >1, some, XStringSet, FALSE, , TRUE
  ## CASE: length >1, some, XStringSet, TRUE, , TRUE
  ## CASE: length >1, some, DNAStringSet, FALSE, FALSE, NA
  ## CASE: length >1, some, DNAStringSet, FALSE, FALSE, "foo"
  ## CASE: length >1, some, DNAStringSet, FALSE, FALSE, FALSE
  ## CASE: length >1, some, DNAStringSet, TRUE, FALSE, FALSE
  ## CASE: length >1, some, DNAStringSet, FALSE, FALSE, TRUE
  ## CASE: length >1, some, DNAStringSet, TRUE, FALSE, TRUE
  ## CASE: length >1, some, DNAStringSet, FALSE, TRUE, FALSE
  ## CASE: length >1, some, DNAStringSet, TRUE, TRUE, FALSE
  ## CASE: length >1, some, DNAStringSet, FALSE, TRUE, TRUE
  ## CASE: length >1, some, DNAStringSet, TRUE, TRUE, TRUE

  ## CASE: XString, DNAStringSet

  ## CASE: XStringSet, DNAString
  
}
