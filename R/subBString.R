### =========================================================================
### Functions for fast substring extraction
### -------------------------------------------------------------------------


### The safe (and exported) version of "BString.substr". Not vectorized.
### We deliberately choose the "NA trick" over defaulting 'start' and 'end'
### to '1' and 'length(x)' because we want to be consistent with the "views"
### function.
subBString <- function(x, start=NA, end=NA, length=NA)
{
    if (!is(x, "BString")) {
        if (!is.character(x) || length(x) != 1 || is.na(x))
            stop("'x' must be a BString (or derived) object or a single string")
        x <- BString(x)
    }
    if (!isNumericOrNAs(start) || length(start) != 1)
        stop("'start' is not a single numeric")
    if (!isNumericOrNAs(end) || length(end) != 1)
        stop("'end' is not a single numeric")
    if (!isNumericOrNAs(length) || length(length) != 1)
        stop("'length' is not a single numeric")
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
    BString.substr(x, as.integer(start), as.integer(end))
}

### For internal use ONLY!
### Solve the views defined by 'start' and 'end' (i.e. get rid of any NAs
### found in 'start' or 'end') by using the information passed in the
### remaining arguments.
### 'start', 'end' and 'width' must be integer vectors of same length N
### (N >= 1) eventually with NAs.
### 'default_start' and 'default_end' must be integer vectors of length >= 1
### with _no_ NAs.
### Return a "valid views data frame" with N views (see definition of the
### "BStringViews" class for the meaning of "valid views data frame").
solveViews <- function(start, end, width, default_start, default_end)
{
    if (!all(is.na(width) | (is.na(start) != is.na(end))))
        stop("invalid start/end/width values")
    ## Solve 'start'
    if (any(is.na(start))) {
        index0 <- is.na(start) & !is.na(width)
        if (any(index0)) {
            start[index0] <- end[index0] - width[index0] + 1L
        }
        if (any(is.na(start))) {
            ## This should be more efficient than recycling 'default_start'
            ## to length N
            i <- (which(is.na(start)) - 1L) %% length(default_start) + 1L
            start[is.na(start)] <- default_start[i]
        }
    }
    ## Solve 'end'
    if (any(is.na(end))) {
        index0 <- is.na(end) & !is.na(width)
        if (any(index0)) {
            end[index0] <- start[index0] + width[index0] - 1L
        }
        if (any(is.na(end))) {
            ## This should be more efficient than recycling 'default_end'
            ## to length N
            i <- (which(is.na(end)) - 1L) %% length(default_end) + 1L
            end[is.na(end)] <- default_end[i]
        }
    }
    if (!all(start <= end))
        stop("'start' and 'end' must verify 'start <= end'")
    data.frame(start=start, end=end)
}

subviews <- function(x, start=NA, end=NA, width=NA, check.limits=TRUE)
{
    if (!is(x, "BStringViews"))
        stop("'x' must be a BStringViews object")
    if (!isNumericOrNAs(start))
        stop("'start' must be a numeric vector")
    if (!isNumericOrNAs(end))
        stop("'end' must be a numeric vector")
    if (!isNumericOrNAs(width))
        stop("'width' must be a numeric vector")
    if (length(x) == 0)
        return(x)
    if (length(start) == 0 || length(start) > length(x))
        stop("length of 'start' must be != 0 and <= length of 'x'")
    if (length(end) == 0 || length(end) > length(x))
        stop("length of 'end' must be != 0 and <= length of 'x'")
    if (length(width) == 0 || length(width) > length(x))
        stop("length of 'width' must be != 0 and <= length of 'x'")
    if (!is.integer(start))
        start <- as.integer(start)
    if (!is.integer(end))
        end <- as.integer(end)
    if (!is.integer(width))
        width <- as.integer(width)
    if (length(start) < length(x))
        start <- rep(start, length.out=length(x))
    if (length(end) < length(x))
        end <- rep(end, length.out=length(x))
    if (length(width) < length(x))
        width <- rep(width, length.out=length(x))
    views <- solveViews(start, end, width, 1L, width(x))
    start <- start(x) + views$start - 1L
    end <- start(x) + views$end - 1L
    if (check.limits && (any(start < 1L) || any(nchar(subject(x)) < end)))
        stop("result contains \"out of limits\" views")
    x@views <- data.frame(start=start, end=end)
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "substr" and "substring" generics and methods.
###

setGeneric(
    "substr", signature="x",
    function(x, start=NA, stop=NA) standardGeneric("substr")
)
setMethod("substr", "BString",
    function(x, start, stop) subBString(x, start, stop)
)
setMethod("substr", "BStringViews",
    function(x, start, stop) subviews(x, start, stop)
)

setGeneric(
    "substring", signature="text",
    function(text, first=NA, last=NA) standardGeneric("substring")
)
setMethod("substring", "BString",
    function(text, first, last) subBString(text, first, last)
)
setMethod("substring", "BStringViews",
    function(text, first, last) subviews(text, first, last)
)

