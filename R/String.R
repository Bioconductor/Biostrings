setClass("BioString",
         representation(alphabet="BioAlphabet",
                        end="integer",
                        start="integer",
                        values="externalptr",
                        initialized="logical"))

if (!isGeneric("BioStringNewValues"))
    setGeneric("BioStringNewValues",
               function(alphabet, length.string)
               standardGeneric("BioStringNewValues"))

setMethod("BioStringNewValues",
          signature(alphabet="NucleotideAlphabet",
                    length.string="numeric"),
          function (alphabet, length.string)
      {
          if (nchar(alphabet@letters) != 5)
              stop("incorrect nucleotide alphabet")
          if (length.string > .Machine$integer.max)
              stop("string too long")
          .Call("BioStringValues", 5, length.string)
      })

setMethod("BioStringNewValues",
          signature(alphabet="AminoAcidAlphabet",
                    length.string="numeric"),
          function (alphabet, length.string)
      {
          if (length.string > .Machine$integer.max)
              stop("string too long")
          .Call("BioStringValues", nchar(alphabet@letters), length.string)
      })

setMethod("BioStringNewValues",
          signature(alphabet="BioPatternAlphabet",
                    length.string="numeric"),
          function (alphabet, length.string)
      {
          alphabet <- alphabet@baseAlphabet
          callGeneric()
      })

setMethod("initialize",
          signature(.Object = "BioString"),
          function (.Object, alphabet,
                    end = 0,
                    start = 0,
                    values=BioStringNewValues(alphabet, end),
                    initialized=!missing(values), ...)
      {
          if (start > end || start < 0)
              stop("invalid range")
          .Object@alphabet <- alphabet
          .Object@start <- as.integer(start)
          .Object@end <- as.integer(end)
          .Object@values <- values
          .Object@initialized <- initialized
          .Object
      })

NucleotideString <- function(src,
                             type=c("DNA", "RNA"),
                             srctype=c("character", "connection"),
                             alphabet=if (type == "DNA") DNAAlphabet()
                                      else RNAAlphabet(),
                             gap=alphabet@gap)
{
    srctype <- match.arg(srctype)
    if (srctype != "character")
        stop("source type not implemented")
    type <- match.arg(type)
    if (!is.character(src) || length(src) != 1)
        stop("src must be a character string")
    if (!is.character(gap) || length(gap) != 1 || nchar(gap) != 1)
        stop("gap must be a single character")
    ans <- new("BioString", alphabet)
    if (gap != alphabet@gap) {
        alphgap <- alphabet@gap
        gapindex <- match(alphgap, names(ans@alphabet@mapping))
        names(ans@alphabet@mapping)[gapindex] <- gap
        substr(ans@alphabet@letters, gapindex, gapindex) <- gap
        ans@alphabet@gap <- gap
        ans <- .Call("setBioString", ans, src)
        ans@initialized <- TRUE
        names(ans@alphabet@mapping)[gapindex] <- alphgap
        substr(ans@alphabet@letters, gapindex, gapindex) <- alphgap
        ans@alphabet@gap <- alphgap
    } else {
        ans <- .Call("setBioString", ans, src)
        ans@initialized <- TRUE
    }
    ans
}

DNAString <- function(src="", gap='-')
{
    NucleotideString(src=src, gap=gap)
}

setMethod("substr",
          signature(x="BioString"),
          function (x, start, stop)
      {
          if (x@start == 0 || start > stop) {
              x@start <- 0
              x@end <- 0
          } else {
              first <- x@start+start-1
              last <- x@start+stop-1
              if (first < x@start)
                  first = x@start;
              if (last >= x@end)
                  last = x@end
              x@start <- as.integer(first)
              x@end <- as.integer(last)
          }
          x
      })

setMethod("as.character",
          signature(x = "BioString"),
          function (x)
      {
          .Call("BioStringToRString", x)
      })

setMethod("nchar",
          signature(x = "BioString"),
          function (x)
      {
          ans <- x@end-x@start+1
          if (ans < 0)
              ans <- 0
          ans
      })

setMethod("show",
          signature(object = "BioString"),
          function (object)
      {
          n <- nchar(object)
          if (n == 0) {
              cat("  Empty biological sequence with alphabet",
                  object@alphabet@letters)
              cat("\n")
          } else {
              cat("  Biological sequence of length", n, "with alphabet",
                  object@alphabet@letters, "and\n")
              if (n > 40) {
                  object <- substr(object, 1, 40)
                  cat(" begining with: ")
              } else {
                  cat(" with value: ")
              }
              cat(as.character(object))
              cat('\n')
          }
      })

if (!isGeneric("matchDNAPattern"))
    setGeneric("matchDNAPattern",
               function(pattern, x, algorithm)
               standardGeneric("matchDNAPattern"))

setMethod("matchDNAPattern",
          signature(pattern="character"),
          function (pattern, x, algorithm)
      {
          pattern <- NucleotideString(pattern,
                                      alphabet=DNAPatternAlphabet())
          callGeneric()
      })

setMethod("matchDNAPattern",
          signature(x="character"),
          function (pattern, x, algorithm)
      {
          x <- NucleotideString(x,
                                alphabet=DNAAlphabet())
          callGeneric()
      })

setMethod("matchDNAPattern",
          signature(pattern="BioString", x="BioString"),
          function (pattern, x, algorithm)
      {
          patternWithPatternAlphabet <-
              is(pattern@alphabet, "BioPatternAlphabet")
          if (!(patternWithPatternAlphabet &&
                identical(pattern@alphabet@baseAlphabet, x@alphabet)) &&
              !(!patternWithPatternAlphabet &&
               identical(pattern@alphabet, x@alphabet)))
              stop("The pattern and the string are based on different alphabets")
          if (pattern@start == pattern@end) {
              patstart <- .Call("LengthOne_exactMatch", pattern, x)
              ans <- matrix(as.integer(0), nrow=length(patstart), ncol=2,
                            dimnames=list(NULL, c("start", "end")))
              ans[, 1] <- patstart
              ans[, 2] <- patstart
              ans
          } else {
              algorithm <-
                  if (missing(algorithm))
                      "boyre-moore"
                  else match.arg(algorithm,
                                 c("boyre-moore"))
              patstart <- switch (algorithm,
                                  "boyre-moore"=.Call("BoyerMoore_exactMatch",
                                  pattern, x),
                                  stop("Unknown algorithm"))
              cbind(start=patstart,
                    end=patstart+(pattern@end-pattern@start))
          }
      })

reverseComplement <-
    function(x)
{
    if (!is(x, "BioString"))
        stop("argument must be a BioString")
    new("BioString", x@alphabet, end=x@end, start=x@start,
        values=.Call("reverseComplementBioString", x))
}
