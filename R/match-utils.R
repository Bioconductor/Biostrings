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

normargWithIndels <- function(with.indels)
{
    if (!isTRUEorFALSE(with.indels))
        stop("'with.indels' must be TRUE or FALSE")
    with.indels
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
### The "neditStartingAt", "neditEndingAt", "isMatchingStartingAt" and
### "isMatchingEndingAt" generic and methods.
###
### 'starting.at' (or 'ending.at') must be integer vectors containing the
### starting (or ending) positions of the pattern relatively to the subject.
### All these functions return a vector of the same length as 'starting.at'
### (or 'ending.at'): the first 2 functions return an integer vector and the
### last 2 functions return a logical vector.
###

setGeneric("neditStartingAt", signature="subject",
    function(pattern, subject, starting.at=1, with.indels=FALSE, fixed=TRUE)
        standardGeneric("neditStartingAt")
)

setGeneric("neditEndingAt", signature="subject",
    function(pattern, subject, ending.at=1, with.indels=FALSE, fixed=TRUE)
        standardGeneric("neditEndingAt")
)

neditAt <- function(pattern, subject, at=1, with.indels=FALSE, fixed=TRUE)
{
    if (!is.numeric(at))
        stop("'at' must be a vector of integers")
    neditStartingAt(pattern, subject, starting.at=at, with.indels=with.indels, fixed=fixed)
}

setGeneric("isMatchingStartingAt", signature="subject",
    function(pattern, subject, starting.at=1,
             max.mismatch=0, with.indels=FALSE, fixed=TRUE)
        standardGeneric("isMatchingStartingAt")
)

setGeneric("isMatchingEndingAt", signature="subject",
    function(pattern, subject, ending.at=1,
             max.mismatch=0, with.indels=FALSE, fixed=TRUE)
        standardGeneric("isMatchingEndingAt")
)

isMatchingAt <- function(pattern, subject, at=1,
                         max.mismatch=0, with.indels=FALSE, fixed=TRUE)
{
    if (!is.numeric(at))
        stop("'at' must be a vector of integers")
    isMatchingStartingAt(pattern, subject, starting.at=at,
                         max.mismatch=max.mismatch, with.indels=with.indels, fixed=fixed)
}

### If 'at.type == 0' then 'at' contains starting positions, otherwise it
### contains ending positions.
.matchPatternAt <- function(pattern, subject, at, at.type,
                            max.mismatch, with.indels, fixed, ans.type)
{
    if (!is(subject, "XString"))
        subject <- XString(NULL, subject)
    if (class(pattern) != class(subject))
        pattern <- XString(class(subject), pattern)
    if (!is.numeric(at)) {
        what <- if (at.type == 0) "starting.at" else "ending.at"
        stop("'", what, "'  must be a vector of integers")
    }
    if (!is.integer(at))
        at <- as.integer(at)
    if (ans.type == 0)
        max.mismatch <- normargMaxMismatch(max.mismatch)
    else
        max.mismatch <- length(pattern)
    with.indels <- normargWithIndels(with.indels)
    fixed <- normargFixed(fixed, class(subject))
    .Call("match_pattern_at",
          pattern, subject, at, at.type,
          max.mismatch, with.indels, fixed, ans.type,
          PACKAGE="Biostrings")
}

.vmatchPatternAt <- function(pattern, subject, at, at.type,
                             max.mismatch, with.indels, fixed, ans.type)
{
    if (!is(subject, "XStringSet"))
        subject <- XStringSet(NULL, subject)
    if (class(pattern) != baseXStringSubtype(subject))
        pattern <- XString(baseXStringSubtype(subject), pattern)
    if (!is.numeric(at)) {
        what <- if (at.type == 0) "starting.at" else "ending.at"
        stop("'", what, "'  must be a vector of integers")
    }
    if (!is.integer(at))
        at <- as.integer(at)
    if (ans.type == 0)
        max.mismatch <- normargMaxMismatch(max.mismatch)
    else
        max.mismatch <- length(pattern)
    with.indels <- normargWithIndels(with.indels)
    fixed <- normargFixed(fixed, baseXStringSubtype(subject))
    .Call("vmatch_pattern_at",
          pattern, subject, at, at.type,
          max.mismatch, with.indels, fixed, ans.type,
          PACKAGE="Biostrings")
}

### Dispatch on 'subject' (see signature of generic).

setMethod("neditStartingAt", "character",
    function(pattern, subject, starting.at=1, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, starting.at, 0L,
                        NA, with.indels, fixed, 1L)
)

setMethod("neditStartingAt", "XString",
    function(pattern, subject, starting.at=1, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, starting.at, 0L,
                        NA, with.indels, fixed, 1L)
)

setMethod("neditStartingAt", "XStringSet",
    function(pattern, subject, starting.at=1, with.indels=FALSE, fixed=TRUE)
        .vmatchPatternAt(pattern, subject, starting.at, 0L,
                         NA, with.indels, fixed, 1L)
)

setMethod("neditEndingAt", "character",
    function(pattern, subject, ending.at=1, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, ending.at, 1L,
                        NA, with.indels, fixed, 1L)
)

setMethod("neditEndingAt", "XString",
    function(pattern, subject, ending.at=1, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, ending.at, 1L,
                        NA, with.indels, fixed, 1L)
)

setMethod("neditEndingAt", "XStringSet",
    function(pattern, subject, ending.at=1, with.indels=FALSE, fixed=TRUE)
        .vmatchPatternAt(pattern, subject, ending.at, 1L,
                         NA, with.indels, fixed, 1L)
)

setMethod("isMatchingStartingAt", "character",
    function(pattern, subject, starting.at=1,
             max.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, starting.at, 0L,
                        max.mismatch, with.indels, fixed, 0L)
)

setMethod("isMatchingStartingAt", "XString",
    function(pattern, subject, starting.at=1,
             max.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, starting.at, 0L,
                        max.mismatch, with.indels, fixed, 0L)
)

setMethod("isMatchingStartingAt", "XStringSet",
    function(pattern, subject, starting.at=1,
             max.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .vmatchPatternAt(pattern, subject, starting.at, 0L,
                         max.mismatch, with.indels, fixed, 0L)
)

setMethod("isMatchingEndingAt", "character",
    function(pattern, subject, ending.at=1,
             max.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, ending.at, 1L,
                        max.mismatch, with.indels, fixed, 0L)
)

setMethod("isMatchingEndingAt", "XString",
    function(pattern, subject, ending.at=1,
             max.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, ending.at, 1L,
                        max.mismatch, with.indels, fixed, 0L)
)

setMethod("isMatchingEndingAt", "XStringSet",
    function(pattern, subject, ending.at=1,
             max.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .vmatchPatternAt(pattern, subject, ending.at, 1L,
                         max.mismatch, with.indels, fixed, 0L)
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
        if (class(pattern) != class(x@subject))
            pattern <- XString(class(x@subject), pattern)
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
        sapplyLength(mismatches)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "coverage" methods.
###

setMethod("coverage", "XStringViews",
    function(x, start=NA, end=NA, weight=1L)
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
        if (!is.numeric(weight) || !(length(weight) %in% c(1, length(x))))
            stop("'weight' must be an integer vector with length 1 or 'length(x)'")
        if (!is.integer(weight))
            weight <- as.integer(weight)
        callNextMethod(x, start=start, end=end, weight=weight)
    }
)

setMethod("coverage", "MaskedXString",
    function(x, start=NA, end=NA, weight=1L)
        coverage(masks(x), start=start, end=end, weight=weight)
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
		ans <- rep.int(0L, end - start + 1L)
        if (is(x, "ByPos_MIndex"))
            .Call("ByPos_MIndex_coverage",
                  endIndex(x), x@width, start, ans,
                  PACKAGE="Biostrings")
        else if (is(x, "ByName_MIndex"))
            .Call("ByName_MIndex_coverage",
                  x@ends_envir, x@width, start, ans,
                  PACKAGE="Biostrings")
        else
            stop("Biostrings internal error: unknown MIndex subtype ", class(x))
        Rle(ans)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Old stuff (Defunct or Deprecated).
###

nmismatchStartingAt <- function(pattern, subject, starting.at=1, fixed=TRUE)
{
    .Deprecated("neditStartingAt")
    neditStartingAt(pattern, subject, starting.at=starting.at, with.indels=FALSE, fixed=fixed)
}

nmismatchEndingAt <- function(pattern, subject, ending.at=1, fixed=TRUE)
{
    .Deprecated("neditEndingAt")
    neditEndingAt(pattern, subject, ending.at=ending.at, with.indels=FALSE, fixed=fixed)
}

isMatching <- function(pattern, subject, start=1, max.mismatch=0, fixed=TRUE)
{
    .Deprecated("isMatchingAt")
    isMatchingAt(pattern, subject, at=start, max.mismatch=max.mismatch, with.indels=FALSE, fixed=fixed)
}

