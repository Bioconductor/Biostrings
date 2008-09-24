### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Compare strings
###

setGeneric("compareStrings", signature = c("pattern", "subject"),
           function(pattern, subject)  standardGeneric("compareStrings"))
setMethod("compareStrings", signature = c(pattern = "character", subject = "character"),
          function(pattern, subject) {
              if (length(pattern) != length(subject))
                  stop("'pattern' and 'subject' must have the same length")
              ncharPattern <- nchar(pattern)
              if (any(ncharPattern != nchar(subject)))
                  stop("'pattern' and 'subject' must have the same number of characters")
              .Call("align_compareStrings", pattern, subject, max(ncharPattern), "+", "-", "?", PACKAGE="Biostrings")
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

setGeneric("consmat", signature="x", function(x, ...)  standardGeneric("consmat"))

setMethod("consmat", "character",
    function(x, freq=FALSE)
    {
        consmat(BStringSet(x), freq=freq)
    }
)

setMethod("consmat", "matrix",
    function(x, freq=FALSE)
    {
        consmat(BStringSet(apply(x, 1, paste, collapse="")), freq=freq)
    }
)

### 'x' must be a list of FASTA records as one returned by readFASTA()
setMethod("consmat", "list",
    function(x, freq=FALSE)
    {
        consmat(BStringSet(FASTArecordsToCharacter(x, use.names=FALSE)), freq=freq)
    }
)

setMethod("consmat", "XStringSet",
    function(x, baseOnly=FALSE, freq=FALSE)
    {
        codes <- codes(super(x), baseOnly=baseOnly)
        if (is.null(names(codes))) {
            names(codes) <- intToUtf8(codes, multiple = TRUE)
            removeUnused <- TRUE
        } else {
            removeUnused <- FALSE
        }
        freq <- .normargFreq(freq)
        ans <- .Call("XStringSet_char_frequency_by_pos",
                x, codes, baseOnly, PACKAGE="Biostrings")
        if (removeUnused) {
            ans <- ans[rowSums(ans) > 0, , drop=FALSE]
        }
        if (freq) {
            ans <- t(t(ans) / colSums(ans))
        }
        ans
    }
)

setMethod("consmat", "XStringViews",
    function(x, baseOnly=FALSE, freq=FALSE)
    {
        y <- XStringViewsToSet(x, use.names=FALSE, verbose=FALSE)
        consmat(y, baseOnly=baseOnly, freq=freq)
    }
)

setMethod("consmat", "AlignedXStringSet",
    function(x, freq=FALSE)
    {
        consmat(aligned(x), freq=freq)
    }
)

setMethod("consmat", "PairwiseAlignment",
    function(x, freq=FALSE)
    {
        consmat(aligned(x), freq=freq)
    }
)
