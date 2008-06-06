### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions (not exported) used by matching functions from other
### files (like matchPattern(), matchPDict(), etc...) to check and normalize
### their arguments.
###

normargMaxMismatch <- function(max.mismatch)
{
    if (!isSingleNumber(max.mismatch))
        stop("'max.mismatch' must be a single integer")
    max.mismatch <- as.integer(max.mismatch)
    if (max.mismatch < 0)
        stop("'max.mismatch' must be a non-negative integer")
    max.mismatch
}

### Return a logical vector of length 2.
normargFixed <- function(fixed, subjectClass)
{
    if (!is.logical(fixed) && !is.character(fixed))
        stop("'fixed' not a logical or character vector")
    if (is.logical(fixed)) {
        if (any(is.na(fixed)))
            stop("'fixed' has NAs")
        fixed_names <- names(fixed)
        if (is.null(fixed_names)) {
            if (!(length(fixed) %in% 1:2))
                stop("when an unnamed logical vector, ",
                     "'fixed' fixed must be of length 1 or 2")
            if (length(fixed) == 1)
                fixed <- c(fixed, fixed)
        } else {
            if (length(fixed) != 2)
                stop("when a named logical vector, 'fixed' must be of length 2")
            if (!setequal(fixed_names, c("pattern", "subject")))
                stop("'fixed' names must be \"pattern\" and \"subject\"")
            fixed <- c(fixed["pattern"], fixed["subject"])
        }
    } else if (is.character(fixed)) {
        if (any(duplicated(fixed)) || !all(fixed %in% c("pattern", "subject")))
            stop("when a character vector, 'fixed' must be ",
                 "a subset of 'c(\"pattern\", \"subject\")' ",
                 "with no duplicated")
        fixed <- c("pattern" %in% fixed, "subject" %in% fixed)
    }
    if (!all(fixed) && !extends(subjectClass, "DNAString") && !extends(subjectClass, "RNAString"))
        stop("'fixed' value only supported for a DNAString or RNAString subject ",
             "(you can only use 'fixed=TRUE' with your subject)")
    fixed
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "nmismatchStartingAt" and "nmismatchEndingAt" generic and methods.
###
### 'starting.at' (or 'ending.at') must be integer vectors containing the
### starting (or ending) positions of the pattern relatively to the subject.
### The two functions return an integer vector of the same length as
### 'starting.at' (or 'ending.at').
###

setGeneric("nmismatchStartingAt", signature="subject",
    function(pattern, subject, starting.at=1, fixed=TRUE)
        standardGeneric("nmismatchStartingAt")
)

setGeneric("nmismatchEndingAt", signature="subject",
    function(pattern, subject, ending.at=1, fixed=TRUE)
        standardGeneric("nmismatchEndingAt")
)

### If 'starting=TRUE' then 'at' contains starting positions, otherwise it
### contains ending positions.
.nmismatchAt <- function(pattern, subject, starting, at, fixed)
{
    if (!is(subject, "XString"))
        subject <- XString(NULL, subject)
    if (class(pattern) != class(subject))
        pattern <- XString(class(subject), pattern)
    if (!is.numeric(at)) {
        what <- if (starting) "starting.at" else "ending.at"
        stop("'", what, "'  must be a vector of integers")
    }
    if (!is.integer(at))
        at <- as.integer(at)
    fixed <- normargFixed(fixed, class(subject))
    .Call("nmismatch_at", pattern, subject,
          starting, at, fixed,
          PACKAGE="Biostrings")
}

### Dispatch on 'subject' (see signature of generic).
setMethod("nmismatchStartingAt", "character",
    function(pattern, subject, starting.at=1, fixed=TRUE)
        .nmismatchAt(pattern, subject, TRUE, starting.at, fixed)
)
setMethod("nmismatchEndingAt", "character",
    function(pattern, subject, ending.at=1, fixed=TRUE)
        .nmismatchAt(pattern, subject, FALSE, ending.at, fixed)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("nmismatchStartingAt", "XString",
    function(pattern, subject, starting.at=1, fixed=TRUE)
        .nmismatchAt(pattern, subject, TRUE, starting.at, fixed)
)
setMethod("nmismatchEndingAt", "XString",
    function(pattern, subject, ending.at=1, fixed=TRUE)
        .nmismatchAt(pattern, subject, FALSE, ending.at, fixed)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "isMatching" generic and methods.
###
### Return a logical vector of the same length as 'start'.
###

setGeneric("isMatching", signature="subject",
    function(pattern, subject, start=1, max.mismatch=0, fixed=TRUE)
        standardGeneric("isMatching")
)

.isMatching <- function(pattern, subject, start, max.mismatch, fixed)
{
    if (!is(subject, "XString"))
        subject <- XString(NULL, subject)
    if (class(pattern) != class(subject))
        pattern <- XString(class(subject), pattern)
    if (!is.numeric(start))
        stop("'start' must be a vector of integers")
    if (!is.integer(start))
        start <- as.integer(start)
    max.mismatch <- normargMaxMismatch(max.mismatch)
    fixed <- normargFixed(fixed, class(subject))
    .Call("is_matching", pattern, subject, start,
          max.mismatch, fixed,
          PACKAGE="Biostrings")
}

### Dispatch on 'subject' (see signature of generic).
setMethod("isMatching", "character",
    function(pattern, subject, start=1, max.mismatch=0, fixed=TRUE)
        .isMatching(pattern, subject, start, max.mismatch, fixed)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("isMatching", "XString",
    function(pattern, subject, start=1, max.mismatch=0, fixed=TRUE)
        .isMatching(pattern, subject, start, max.mismatch, fixed)
)


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
            l <- XString.substr(pattern, i, i)
            cp <- isMatching(l, subject, j, max.mismatch=0, fixed=fixed)
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
        warning("views in 'x' don't have a width equal to pattern length")
    lapply(1:length(x),
           function(i) .bsMismatch(pattern, subject(x), start(x)[i], fixed))
}

setGeneric("mismatch", signature="x",
    function(pattern, x, fixed=TRUE) standardGeneric("mismatch")
)

### Typical use:
###   mp <- matchPattern("TGA", DNAString("GTGACGTGCAT"), max.mismatch=2)
###   mismatch("TGA", mp)
### Dispatch on 'x' (see signature of generic).
setMethod("mismatch", "XStringViews",
    function(pattern, x, fixed)
    {
        if (class(pattern) != class(x@subject))
            pattern <- XString(class(x@subject), pattern)
        .mismatch(pattern, x, fixed)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "coverage" generic and methods.
###

setGeneric("coverage", signature="x",
    function(x, start=NA, end=NA) standardGeneric("coverage")
)

setMethod("coverage", "IRanges",
    function(x, start=NA, end=NA)
    {
        if (!isSingleNumber(start))
            stop("'start' must be a single integer")
        if (!is.integer(start))
            start <- as.integer(start)
        if (!isSingleNumber(end))
            stop("'end' must be a single integer")
        if (!is.integer(end))
            end <- as.integer(end)
        width <- end - start + 1L
        if (width < 0)
            stop("'end' must be >= 'start' - 1")
        x1 <- shift(restrict(x, start=start, end=end), 1L - start)
        ans <- integer(width)
        for (i in seq_len(length(x1))) {
            start1 <- start(x1)[i]
            end1 <- end(x1)[i]
            if (end1 < start1)
                next
            ii <- start1:end1
            ans[ii] <- ans[ii] + 1L
        }
        ans
    }
)

setMethod("coverage", "MaskCollection",
    function(x, start=NA, end=NA)
    {
        if (!isSingleNumberOrNA(start))
            stop("'start' must be a single integer or NA")
        if (!is.integer(start))
            start <- as.integer(start)
        if (is.na(start))
            start <- 1L
        if (!isSingleNumberOrNA(end))
            stop("'end' must be a single integer or NA")
        if (!is.integer(end))
            end <- as.integer(end)
        if (is.na(end))
            end <- width(x)
        width <- end - start + 1L
        if (width < 0)
            stop("'end' must be >= 'start' - 1")
        ans <- integer(width)
        for (i in seq_len(length(x)))
            ans <- ans + coverage(x[[i]], start=start, end=end)
        ans
    }
)

setMethod("coverage", "XStringViews",
    function(x, start=NA, end=NA)
    {
        if (!isSingleNumberOrNA(start))
            stop("'start' must be a single integer or NA")
        if (!is.integer(start))
            start <- as.integer(start)
        if (is.na(start))
            start <- 1L
        if (!isSingleNumberOrNA(end))
            stop("'end' must be a single integer or NA")
        if (!is.integer(end))
            end <- as.integer(end)
        if (is.na(end))
            end <- length(subject(x))
        callNextMethod(x, start=start, end=end)
    }
)

setMethod("coverage", "MaskedXString",
    function(x, start=NA, end=NA)
        coverage(masks(x), start=start, end=end)
)

setMethod("coverage", "MIndex",
    function(x, start=NA, end=NA)
    {
        if (!isSingleNumber(start))
            stop("'start' must be a single integer")
        if (!is.integer(start))
            start <- as.integer(start)
        if (!isSingleNumber(end))
            stop("'end' must be a single integer")
        if (!is.integer(end))
            end <- as.integer(end)
        width <- end - start + 1L
        if (width < 0)
            stop("'end' must be >= 'start' - 1")
        if (is(x, "ByPos_MIndex"))
            return(.Call("ByPos_MIndex_coverage",
                         x@ends, x@width, start, end,
                         PACKAGE="Biostrings"))
        if (is(x, "ByName_MIndex"))
            return(.Call("ByName_MIndex_coverage",
                         x@length, x@ends_envir, x@width, start, end,
                         PACKAGE="Biostrings"))
        stop("Biostrings internal error: unknown MIndex subtype ", class(x))
    }
)

