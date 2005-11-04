#Copyright (C) 2003-2004 by Saikat DebRoy
setClass("BioString",
         representation(alphabet="BioAlphabet",
                        offsets="matrix",
                        values="externalptr",
                        initialized="logical"))

setMethod("as.matrix",
          signature(x = "BioString"),
          function (x)
          x@offsets)

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
          .Call("BioStringValues", 5, length.string,
                PACKAGE="Biostrings")
      })

setMethod("BioStringNewValues",
          signature(alphabet="AminoAcidAlphabet",
                    length.string="numeric"),
          function (alphabet, length.string)
      {
          if (length.string > .Machine$integer.max)
              stop("string too long")
          .Call("BioStringValues", nchar(alphabet@letters),
                length.string, PACKAGE="Biostrings")
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
                    offsets = cbind(1, 0),
                    values=BioStringNewValues(alphabet, end),
                    initialized=!missing(values), ...)
      {
          if (!is.matrix(offsets) || ncol(offsets) != 2)
              stop("invalid offsets")
          storage.mode(offsets) <- "integer"
          start <- offsets[, 1]
          end <- offsets[, 2]
          offsets[, 1] <- start <- ifelse(start <= 0, 1, start)
          offsets[start > end, ] <- 1:0
          storage.mode(offsets) <- "integer"
          .Object@alphabet <- alphabet
          .Object@offsets <- offsets
          .Object@values <- values
          .Object@initialized <- initialized
          .Object
      })

NucleotideString <- function(src,
                             type=c("DNA", "RNA"),
                             srctype=c("character", "connection"),
                             alphabet=if (type == "DNA") DNAPatternAlphabet()
                                      else RNAPatternAlphabet(),
                             gap=alphabet@gap)
{
    srctype <- match.arg(srctype)
    if (srctype != "character")
        stop("source type not implemented")
    type <- match.arg(type)
    if (!is.character(src))
        stop("src must be a character string")
    if (!is.character(gap) || length(gap) != 1 || nchar(gap) != 1)
        stop("gap must be a single character")
    if (gap != alphabet@gap) {
        savealph <- alphabet
        if (gap %in% names(alphabet@mapping)) {
            if (!is(alphabet, "BioPatternAlphabet") ||
                gap %in% names(alphabet@baseAlphabet@mapping))
                stop("gap character conflicts with other alphabet characters")
            gindex <- match(gap, names(alphabet@mapping))
            names(alphabet@mapping)[gindex] <- '0'
            substr(alphabet@letters, gindex, gindex) <- '0'
        }
        alphgap <- alphabet@gap
        gapletter(alphabet) <- gap
        ans <- .Call("setBioString", new("BioString", alphabet), src,
                     PACKAGE="Biostrings")
        ans@initialized <- TRUE
        ans@alphabet <- savealph
    } else {
        ans <- .Call("setBioString", new("BioString", alphabet), src,
                     PACKAGE="Biostrings")
        ans@initialized <- TRUE
    }
    ans
}

DNAString <- function(src="", gap='-')
{
    NucleotideString(src=src, gap=gap)
}

setMethod("length",
          signature(x = "BioString"),
          function (x)
      {
          nrow(x@offsets)
      })

setMethod("[",
          signature(x = "BioString"),
          function (x, i, j, ..., drop)
      {
          if (nargs() == 1 || nargs() == 2 && missing(i))
              x
          else if (nargs() == 2 && !missing(i)) {
              x@offsets <- x@offsets[i, , drop = FALSE]
              x
          } else stop("invalid subsetting")
      })

setMethod("[[",
          signature(x = "BioString"),
          function (x, i, j, ..., drop)
      {
          if (nargs() != 2 && (nargs() != 3 || !missing(j)))
              stop("incorrect number of subscripts")
          if (missing(i))
              stop("invalid subscript type")
          x@offsets <- cbind(x@offsets[[i, 1]], x@offsets[[i, 2]])
          storage.mode(x@offsets) <- "integer"
          x
      })

setMethod("substr",
          signature(x = "BioString"),
          function (x, start, stop)
      {
          .Call("BioString_substring", x, start, stop, FALSE,
                PACKAGE="Biostrings")
      })

setMethod("substring",
          signature(text = "BioString"),
          function (text, first, last)
      {
          .Call("BioString_substring", text, first, last, TRUE,
                PACKAGE="Biostrings")
      })

setMethod("as.character",
          signature(x = "BioString"),
          function (x)
      {
          .Call("BioStringToRString", x, PACKAGE="Biostrings")
      })

##FIXME: nchar now does more careful things with multibyte character
##strings. I think we can ignore this here - but just in case...
setMethod("nchar",
          signature(x = "BioString", type="character"),
          function (x, type)
      {
          ans <- x@offsets[,2]-x@offsets[,1]+1
          as.integer(ifelse(ans < 0, as.integer(0), ans))
      })

setMethod("nchar",
          signature(x = "BioString", type="missing"),
          function (x, type)
      {
          ans <- x@offsets[,2]-x@offsets[,1]+1
          as.integer(ifelse(ans < 0, as.integer(0), ans))
      })

setMethod("show",
          signature(object = "BioString"),
          function (object)
      {
          nvec <- length(object)
          if (nvec == 0) {
              cat("    Object of class BioString with\n")
              if (is(object@alphabet, "BioPatternAlphabet"))
                  cat("Pattern alphabet: ")
              else cat("Alphabet: ")
              cat(object@alphabet@letters)
              cat("\nLength: 0\n")
          } else {
              cat("    Object of class BioString with\n")
              if (is(object@alphabet, "BioPatternAlphabet"))
                  cat("Pattern alphabet: ")
              else cat("Alphabet: ")
              cat(object@alphabet@letters)
              cat("\nValues:\n")
              printsome <- nvec > 30
              i <- 1
              while (i <= nvec) {
                  if (printsome && i == 15) {
                      cat("\n\n")
                      i <- nvec-13
                  }
                  if (i >= 10 || nvec < 10)
                      cat(paste('[', i, ']', sep=''))
                  else cat(paste(' [', i, ']', sep=''))
                  cat(' ')
                  tmp <- object[i]
                  n <- nchar(tmp)
                  if (n > 40) {
                      cat(as.character(substr(tmp, 1, 18)))
                      cat("....")
                      cat(as.character(substr(tmp, n-17, n)))
                  } else cat(as.character(tmp))
                  cat('\n')
                  i <- i+1
              }
          }
      })

if (!isGeneric("matchDNAPattern"))
    setGeneric("matchDNAPattern",
               function(pattern, x, algorithm, mismatch)
               standardGeneric("matchDNAPattern"))

setMethod("matchDNAPattern",
          signature(pattern="character"),
          function (pattern, x, algorithm, mismatch)
      {
          pattern <- NucleotideString(pattern,
                                      alphabet=DNAPatternAlphabet())
          callGeneric()
      })

setMethod("matchDNAPattern",
          signature(x="character"),
          function (pattern, x, algorithm, mismatch)
      {
          x <- DNAString(x)
          callGeneric()
      })

setMethod("matchDNAPattern",
          signature(pattern="BioString", x="BioString"),
          function (pattern, x, algorithm, mismatch)
      {
          if (missing(mismatch))
              mismatch <- 0
          else if (mismatch < 0)
              stop("mismatch must be a non-negative integer")
          pattern <- NucleotideString(as.character(pattern),
                                      alphabet=pattern@alphabet)
          patalph <- pattern@alphabet
          if (is(patalph, "BioPatternAlphabet"))
              patalph <- patalph@baseAlphabet
          xalph <- x@alphabet
          if (is(xalph, "BioPatternAlphabet"))
              xalph <- xalph@baseAlphabet
          if (!identical(patalph, xalph))
              stop("The pattern and the string are based on different alphabets")
          algorithm <-
              if (missing(algorithm)) {
                  if (mismatch == 0 && nchar(pattern) > .Machine$sizeof.long*8)
                      "boyer-moore"
                  else "shift-or"
              } else match.arg(algorithm,
                             c("boyer-moore",
                               "forward-search",
                               "shift-or"))
          if (mismatch > 0 && algorithm != "shift-or")
              stop("only mismatch algorithm is shift-or")
          switch(algorithm,
                 "boyer-moore"=.Call("BoyerMoore_exactMatch", pattern,
                                     x, PACKAGE="Biostrings"),
                 "forward-search"=.Call("ForwardSearch_exactMatch",
                                        pattern, x, PACKAGE="Biostrings"),
                 "shift-or"=.Call('ShiftOr_inexactMatch', pattern, x,
                 mismatch, mismatch, mismatch, mismatch, PACKAGE="Biostrings"),
                 stop("Unknown algorithm"))
      })

reverseComplement <-
    function(x)
{
    if (!is(x, "BioString"))
        stop("argument must be a BioString")
    .Call("reverseComplementBioString", x, PACKAGE="Biostrings")
}

if (!isGeneric("allSameLetter"))
    setGeneric("allSameLetter",
               function(x, letter)
               standardGeneric("allSameLetter"))

setMethod("allSameLetter",
          signature(x="BioString", letter="character"),
          function (x, letter)
      {
          callGeneric(x, NucleotideString(letter, alphabet=x@alphabet))
      })

setMethod("allSameLetter",
          signature(x="character"),
          function (x, letter)
      {
          callGeneric(DNAString(x), letter)
      })

setMethod("allSameLetter",
          signature(x = "BioString", letter="BioString"),
          function (x, letter)
      {
          .Call("allSameLetter", x, letter, TRUE, PACKAGE="Biostrings")
      })

if (!isGeneric("anySameLetter"))
    setGeneric("anySameLetter",
               function(x, letter)
               standardGeneric("anySameLetter"))

setMethod("anySameLetter",
          signature(x="BioString", letter="character"),
          function (x, letter)
      {
          callGeneric(x, NucleotideString(letter, alphabet=x@alphabet))
      })

setMethod("anySameLetter",
          signature(x="character"),
          function (x, letter)
      {
          callGeneric(DNAString(x), letter)
      })

setMethod("anySameLetter",
          signature(x = "BioString", letter="BioString"),
          function (x, letter)
      {
          .Call("allSameLetter", x, letter, FALSE, PACKAGE="Biostrings")
      })

alphabetFrequency <-
    function (x, baseOnly=TRUE)
{
    xtable <- .Call('baseFrequency', x, PACKAGE='Biostrings')
    if (baseOnly && is(x@alphabet, "BioPatternAlphabet")) {
        baseLetters <- sort(names(x@alphabet@baseAlphabet@mapping))
        xtable <- c(xtable[baseLetters],
                    sum(xtable)-sum(xtable[baseLetters]))
        names(xtable) <- c(baseLetters, 'Other')
    }
    xtable
}

sortDNAString <-
    function(x, prefixLength = max(nchar(x)))
{
    .Call('SortDNAString', x, prefixLength, PACKAGE="Biostrings")
}

DNASuffixArray <-
    function(x, prefixLength = max(nchar(x)))
{
    .Call('DNASuffixArray', x, prefixLength, PACKAGE="Biostrings")
}

setClass("LongestCommonPrefix",
         representation(abovediag="integer"))

LongestCommonPrefix <-
    function(x)
{
    new("LongestCommonPrefix",
        abovediag=.Call('longestCommonPrefixSuffixArray', x,
        PACKAGE='Biostrings'))
}

setMethod("[",
          signature(x = "LongestCommonPrefix", i = "matrix",
                    j = "missing", drop = "missing"),
          function (x, i, j, ..., drop)
      {
          if ((storage.mode(i) != "integer" && storage.mode(i) !=
               "double") || ncol(i) != 2)
              stop("invalid matrix index")
          ans <- integer(nrow(i))
          for (j in seq(along=ans)) {
              starti <- i[j, 1]
              endi <- i[j, 2]
              if (starti != endi) {
                  if (starti > endi) {
                      starti <- endi
                      endi <- i[j, 1]
                  }
                  ans[j] <- as.integer(min(x@abovediag[starti:(endi-1)]))
              }
          }
          ans
      })

