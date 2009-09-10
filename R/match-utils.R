### =========================================================================
### Miscellaneous utility functions operating on the matches returned by a
### high-level matching function like matchPattern(), matchPDict(), etc...
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mismatch()
###

### Helper function used by .mismatch()
### Returns a vector of the positions of mismatches of 'pattern'
### in a view on 'subject' starting at 'start' and whose width is length(pattern).
.bsMismatch <- function(pattern, subject, start, fixed)
{
    mm <- integer(0)
    j0 <- start - as.integer(1)
    for (i in seq_len(length(pattern))) {
        j <- j0 + i
        if (j < 1 || j > length(subject)) {
            mm <- c(mm, i)
        } else {
            l <- subseq(pattern, start=i, end=i)
            cp <- isMatchingAt(l, subject, at=j, fixed=fixed)
            if (cp == 0)
                mm <- c(mm, i)
        }
    }
    mm
}

.mismatch <- function(pattern, x, fixed)
{
    if (length(x) == 0)
        return(list())
    if (any(width(x) != length(pattern)))
        warning("some views in 'x' have a width that differs from 'length(pattern)'")
    lapply(1:length(x),
           function(i) .bsMismatch(pattern, subject(x), start(x)[i], fixed))
}

setGeneric("mismatch", signature=c("pattern", "x"),
    function(pattern, x, fixed=TRUE) standardGeneric("mismatch")
)

### Typical use:
###   mp <- matchPattern("TGA", DNAString("GTGACGTGCAT"), max.mismatch=2)
###   mismatch("TGA", mp)
### Dispatch on 'x' (see signature of generic).
setMethod("mismatch", c(pattern="ANY", x="XStringViews"),
    function(pattern, x, fixed)
    {
        pattern <- normargPattern(pattern, x)
        .mismatch(pattern, x, fixed)
    }
)

setGeneric("nmatch", signature=c("pattern", "x"),
    function(pattern, x, fixed=TRUE) standardGeneric("nmatch")
)

setMethod("nmatch", c(pattern="ANY", x="XStringViews"),
    function(pattern, x, fixed)
    {
        funCall <- match.call()
        funCall[[1]] <- as.name("nmismatch")
        nchar(pattern) - eval(funCall, sys.parent())
    }
)

setGeneric("nmismatch", signature=c("pattern", "x"),
    function(pattern, x, fixed=TRUE) standardGeneric("nmismatch")
)

setMethod("nmismatch", c(pattern="ANY", x="XStringViews"),
    function(pattern, x, fixed)
    {
        funCall <- match.call()
        funCall[[1]] <- as.name("mismatch")
        mismatches <- eval(funCall, sys.parent())
        elementLengths(mismatches)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "coverage" methods.
###

setMethod("coverage", "MaskedXString",
    function(x, start=NA, end=NA, shift=0L, width=NULL, weight=1L)
        coverage(masks(x), start=start, end=end, shift=shift, width=width, weight=weight)
)

setMethod("coverage", "MIndex",
    function(x, start=NA, end=NA, shift=0L, width=NULL, weight=1L)
        coverage(unlist(x), start=start, end=end, shift=shift, width=width, weight=weight)
)

