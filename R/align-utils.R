### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Compare strings
###

setGeneric("compareStrings", signature = c("pattern", "subject"),
           function(pattern, subject)  standardGeneric("compareStrings"))
setMethod("compareStrings", signature = c(pattern = "character", subject = "character"),
          function(pattern, subject) {
              if (length(pattern) != length(subject))
                  stop("'pattern' and 'subject' must have the same length")
              if (any(nchar(pattern) != nchar(subject)))
                  stop("'pattern' and 'subject' must have the same number of characters")
              output <- rep("", length(pattern))
              for (i in 1:length(pattern)) {
                  patternCodes <- charToRaw(pattern[i])
                  subjectCodes <- charToRaw(subject[i])
                  insertionLocations <- subjectCodes == charToRaw("-")
                  deletionLocations <- patternCodes == charToRaw("-")
                  nonIndels <- (!insertionLocations & !deletionLocations)
                  current <- rawToChar(subjectCodes, multiple = TRUE)
                  current[insertionLocations] <- "+"
                  current[deletionLocations] <- "-"
                  current[nonIndels & patternCodes != subjectCodes] <- "?"
                  output[i] <- paste(current, collapse = "")
              }
              output
          })
setMethod("compareStrings", signature = c(pattern = "AlignedXStringSet", subject = "AlignedXStringSet"),
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

setMethod("consmat", "PairwiseAlignment",
    function(x, freq=TRUE)
    {
        consmat(as.character(x), freq=freq)
    }
)

