### =========================================================================
### Functions for fast substring extraction
### -------------------------------------------------------------------------


### An attempt to have an substr() equivalent for XString objects with some
### important differences though:
###   - subBString() is NOT vectorized with respect to 'x'
###   - the arguments have different names
###   - there is one extra argument: the 'length' arg
###   - NAs are valid starts, ends or lengths
###   - using subBString() for replacement is not supported
### We deliberately choose the "NA trick" over defaulting 'start' and 'end'
### to '1' and 'length(x)' because we want to be consistent with the "views"
### function.
subXString <- function(x, start=NA, end=NA, length=NA)
{
    if (!is(x, "XString")) {
        if (!isSingleString(x))
            stop("'x' must be an XString object or a single string")
        x <- XString(NULL, x)
    }
    if (!isNumericOrNAs(start) || length(start) != 1)
        stop("'start' must be a single number")
    if (!isNumericOrNAs(end) || length(end) != 1)
        stop("'end' must be a single number")
    if (!isNumericOrNAs(length) || length(length) != 1)
        stop("'length' must be a single number")
    if (!is.na(length)) {
        if (is.na(start) && is.na(end))
            stop("you must specify 'start' or 'end'")
        if (is.na(end))
            end <- start + length - 1
        else if (is.na(start))
            start <- end - length + 1
        else if (length != end - start + 1)
            stop("incompatible 'start', 'end' and 'length' values")
    } else {
        if (is.na(start))
            start <- 1
        if (is.na(end))
            end <- x@length
    }
    ## This is NA-proof (well, 'start' and 'end' can't be NAs anymore...)
    if (!isTRUE(1 <= start && start <= end && end <= length(x)))
        stop("'start' and 'end' must verify '1 <= start <= end <= length(x)'")
    XString.substr(x, as.integer(start), as.integer(end))
}

subBString <- function(x, start=NA, end=NA, length=NA)
{
    .Deprecated("subXString")
    subXString(x, start=start, end=end, length=length)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "substr" and "substring" generics and methods.
###

setGeneric("substr", signature="x",
    function(x, start=NA, stop=NA) standardGeneric("substr")
)
setMethod("substr", "XString",
    function(x, start=NA, stop=NA) subXString(x, start, stop)
)
setMethod("substr", "MaskedXString",
    function(x, start=NA, stop=NA) subXString(unmasked(x), start, stop)
)

setGeneric("substring", signature="text",
    function(text, first=NA, last=NA) standardGeneric("substring")
)
setMethod("substring", "XString",
    function(text, first=NA, last=NA) subXString(text, first, last)
)
setMethod("substring", "MaskedXString",
    function(text, first=NA, last=NA) subXString(unmasked(text), first, last)
)
