### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Compare strings
###

setGeneric("compareStrings", signature = c("pattern", "subject"),
           function(pattern, subject)  standardGeneric("compareStrings"))
setMethod("compareStrings", signature = c(pattern = "character", subject = "character"),
          function(pattern, subject) {
              if (nchar(pattern) != nchar(subject))
                  stop("'pattern' and 'subject' must have the same number of characters")
              patternCodes <- charToRaw(pattern)
              subjectCodes <- charToRaw(subject)
              nonGaps <- (patternCodes != charToRaw("-")) & (subjectCodes != charToRaw("-"))
              matches <- sum(patternCodes[nonGaps] == subjectCodes[nonGaps])
			  c(matches = matches, mismatches = sum(nonGaps) - matches)
          })
setMethod("compareStrings", signature = c(pattern = "AlignedXString", subject = "AlignedXString"),
          function(pattern, subject) {
              compareStrings(as.character(pattern), as.character(subject))
          })
setMethod("compareStrings", signature = c(pattern = "PairwiseAlignment", subject = "missing"),
		  function(pattern, subject) {
			  compareStrings(pattern@pattern, pattern@subject)
		  })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Alignment consensus matrix
###

setGeneric("consmat", signature="x", function(x, freq=TRUE)  standardGeneric("consmat"))

setMethod("consmat", "character",
    function(x, freq=TRUE)
    {
        nrow <- length(x)
        if (nrow == 0)
            stop("'x' must contain at least 1 string")
        if (any(is.na(x)))
            stop("NAs are not allowed in 'x'")
        nchars <- nchar(x)
        ncol <- nchars[1]
        if (!all(nchars == ncol))
            stop("'x' elements are not equal-length strings")
        allletters <- unlist(strsplit(x, NULL))
        pos <- rep.int(seq_len(ncol), nrow)
        ans <- table(letter=allletters, pos=pos)
        if (freq)
            ans <- ans / nrow
        ans
    }
)

setMethod("consmat", "matrix",
    function(x, freq=TRUE)
    {
        cat("coming soon...")
    }
)

### 'x' must be a list of FASTA records as one returned by readFASTA()
setMethod("consmat", "list",
    function(x, freq=TRUE)
    {
        consmat(FASTArecordsToCharacter(x, use.names=FALSE), freq=freq)
    }
)

setMethod("consmat", "XStringSet",
    function(x, freq=TRUE)
    {
        consmat(as.character(x, use.names=FALSE), freq=freq)
    }
)

setMethod("consmat", "XStringViews",
    function(x, freq=TRUE)
    {
        consmat(as.character(x, use.names=FALSE), freq=freq)
    }
)

setMethod("consmat", "AlignedXString",
    function(x, freq=TRUE)
    {
        consmat(as.character(x), freq=freq)
    }
)

