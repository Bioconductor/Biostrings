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

### 'x' must be a list of FASTA records as one returned by readFASTA()
setMethod("consmat", "list",
    function(x, freq=TRUE)
    {
        consmat(FASTArecordsToCharacter(x, use.names=FALSE), freq=freq)
    }
)

setMethod("consmat", "BStringViews",
    function(x, freq=TRUE)
    {
        consmat(as.character(x, use.names=FALSE), freq=freq)
    }
)

