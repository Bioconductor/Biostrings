### =========================================================================
### The mask() generic and methods
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###
### Some notes:
###
###   o BStringViews object 'x' is "normal" iff the underlying IRanges object
###     is normal and 'x' has no "out of limits" views.
###
###   o It can been shown that
###       length(reduce(x)) <= (nchar(subject(x)) + 1) / 2
###
### A formal definition of mask(x):
###
###   For any given BStringViews object 'x', 'y <- mask(x)' is the unique
###   normal BStringViews object that:
###     (a) subject(y) == subject(x);
###     (b) 'y' views cover the regions in the subject that are not covered
###         by 'x' views.
###
###   Note that 'mask(mask(x))' is reduce(x).
###

### We wan't to be able to do this:
###   > mask(BString("AbcbcbDE"), 2, 6)   # masking by position
### or this:
###   > mask(BString("AbcbcbDE"), "bcb")  # masking by content
setMethod("mask", "XString",
    function(x, start, end, pattern)
    {
        if (missing(pattern)) {
            if (missing(start))
                start <- NA
            if (isNumericOrNAs(start)) {
                if (missing(end))
                    end <- NA
                ans <- as(x, paste("Masked", class(x), sep=""))
                ans@mask <- toNormalIRanges(views(x, start, end))
                return(ans)
            }
            if (!missing(end))
                stop("invalid 'start' argument")
            pattern <- start
        } else {
            if (!missing(start) || !missing(end))
                stop("can't give 'start' (or 'end') when 'pattern' is given")
        }
        m <- matchPattern(pattern, x)
        mask(x, start(m), end(m))
    }
)

setMethod("mask", "BStringViews",
    function(x, start, end, ...)
    {
        if (missing(start))
            start <- 1L
        if (missing(end))
            end <- nchar(subject(x))
        callNextMethod(x, start, end, ...)
    }
)

setMethod("mask", "MaskedXString",
    function(x, start, end, pattern)
    {
        if (!missing(start) || !missing(end) || !missing(pattern))
            stop("calling mask() on a MaskedXString object with more than 1 arg is not supported yet")
        x@mask <- mask(x@mask, 1L, length(x))
    }
)

setMethod("mask", "character",
    function(x, start, end, pattern)
    {
        if (length(x) != 1 || is.na(x))
            stop("can't mask a character vector that is not a single string")
        if (missing(pattern))
            mask(BString(x), start, end)
        else
            mask(BString(x), start, end, pattern)
    }
)

