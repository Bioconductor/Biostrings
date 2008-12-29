### =========================================================================
### Functions for fast substring extraction
### -------------------------------------------------------------------------


subXString <- function(x, start=NA, end=NA, length=NA)
{
    .Deprecated("subseq")
    if (!is(x, "XString")) {
        if (!isSingleString(x))
            stop("'x' must be an XString object or a single string")
        x <- XString(NULL, x)
    }
    subseq(x, start=start, end=end, width=length)
}

subBString <- function(x, start=NA, end=NA, length=NA)
    .Defunct("subseq")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "substr" and "substring" generics and methods.
###

setMethod("substr", "XString",
    function(x, start=NA, stop=NA) subseq(x, start=start, end=stop)
)

setMethod("substring", "XString",
    function(text, first=NA, last=NA) subseq(text, start=first, end=last)
)

