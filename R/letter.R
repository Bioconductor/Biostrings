### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "letter" new generic.
###

### Return a character vector.
setGeneric("letter", signature="x",
    function(x, i) standardGeneric("letter")
)

### Return a character vector of the same length as 'x'
### where all non-NA elements have length(i) characters.
setMethod("letter", "character",
    function(x, i)
    {
        if (!is.numeric(i) || anyNA(i))
            stop("subscript 'i' must be an NA-free numeric vector")
        if (length(x) == 0)
            return(character(0))
        if (length(i) == 0)
            return(character(length(x)))
        noNA_x <- x[!is.na(x)]
        imax <- min(nchar(noNA_x))
        if (!all(i >= 1) || !all(i <= imax))
            stop("subscript out of bounds")
        ## Which one is faster? Looping on i with substr or looping
        ## on noNA_x with substring? It depends...
        if (length(x) >= length(i)) {
            tmp <- lapply(i, function(i1) substr(noNA_x, i1, i1))
            x[!is.na(x)] <- do.call(paste, c(tmp, sep=""))
        } else {
            x[!is.na(x)] <- sapply(
                                noNA_x,
                                function(x1) paste(substring(x1, i, i), collapse=""),
                                USE.NAMES=FALSE
                            )
        }
        x
    }
)

### Return a character vector of length 1.
setMethod("letter", "XString",
    function(x, i)
    {
        if (!is.numeric(i) || anyNA(i))
            stop("subscript 'i' must be an NA-free numeric vector")
        if (!is.integer(i))
            i <- as.integer(i)
        imax <- length(x)
        if (!all(i >= 1) || !all(i <= imax))
            stop("subscript out of bounds")
        extract_character_from_XString_by_positions(x, i, collapse=TRUE)
    }
)

### Return a character vector of the same length as 'x'.
setMethod("letter", "XStringViews",
    function(x, i)
    {
        if (!is.numeric(i) || anyNA(i))
            stop("subscript 'i' must be an NA-free numeric vector")
        if (!is.integer(i))
            i <- as.integer(i)
        if (length(x) == 0)
            return(character(0))
        imax <- min(width(x))
        if (!all(i >= 1L) || !all(i <= imax))
            stop("subscript out of bounds")
        sapply(x, letter, i)
    }
)

setMethod("letter", "MaskedXString",
    function(x, i)
        letter(unmasked(x), i)
)

