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

### We wan't to be able to do this:
###   > mask(BString("AbcbcbDE"), 2, 6)
### or this:
###   > mask(BString("AbcbcbDE"), "bcb")
setMethod("mask", "XString",
    function(x, start, end, pattern)
    {
        if (missing(pattern)) {
            if (missing(start))
                start <- NA
            if (isNumericOrNAs(start)) {
                if (missing(end))
                    end <- NA
                return(mask(views(x, start, end)))
            }
            if (!missing(end))
                stop("invalid 'start' argument")
            pattern <- start
        } else {
            if (!missing(start) || !missing(end))
                stop("can't give 'start' (or 'end') when 'pattern' is given")
        }
        mask(matchPattern(pattern, x))
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

