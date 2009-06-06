### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions (not exported) used by matching functions from other
### files (like matchPattern(), matchPDict(), etc...) to check and normalize
### their arguments.
###

normargPattern <- function(pattern, subject, argname="pattern")
{
    subject_baseclass <- xsbaseclass(subject)
    if (is(pattern, "XString")) {
        if (xsbaseclass(pattern) == subject_baseclass)
            return(pattern)
    } else if (!isSingleString(pattern))
        stop("'", argname, "' must be a single string or an XString object")
    pattern <- try(XString(xsbasetype(subject), pattern))
    if (is(pattern, "try-error"))
        stop("could not turn '", argname, "' into a ", subject_baseclass, " instance")
    pattern
}

normargMaxMismatch <- function(max.mismatch, argname="max.mismatch")
{
    if (!isSingleNumber(max.mismatch))
        stop("'", argname, "' must be a single integer")
    max.mismatch <- as.integer(max.mismatch)
    if (max.mismatch < 0)
        stop("'", argname, "' must be a non-negative integer")
    max.mismatch
}

normargWithIndels <- function(with.indels, argname="with.indels")
{
    if (!isTRUEorFALSE(with.indels))
        stop("'", argname, "' must be TRUE or FALSE")
    with.indels
}

### Return a logical vector of length 2.
normargFixed <- function(fixed, subject, argname="fixed")
{
    if (!is.logical(fixed) && !is.character(fixed))
        stop("'", argname, "' not a logical or character vector")
    if (is.logical(fixed)) {
        if (any(is.na(fixed)))
            stop("'", argname, "' has NAs")
        fixed_names <- names(fixed)
        if (is.null(fixed_names)) {
            if (!(length(fixed) %in% 1:2))
                stop("when an unnamed logical vector, '", argname,
                     "' fixed must be of length 1 or 2")
            if (length(fixed) == 1)
                fixed <- c(fixed, fixed)
        } else {
            if (length(fixed) != 2)
                stop("when a named logical vector, '", argname, "' must be of length 2")
            if (!setequal(fixed_names, c("pattern", "subject")))
                stop("'", argname, "' names must be \"pattern\" and \"subject\"")
            fixed <- c(fixed["pattern"], fixed["subject"])
        }
    } else if (is.character(fixed)) {
        if (any(duplicated(fixed)) || !all(fixed %in% c("pattern", "subject")))
            stop("when a character vector, '", argname, "' must be ",
                 "a subset of 'c(\"pattern\", \"subject\")' ",
                 "with no duplicated")
        fixed <- c("pattern" %in% fixed, "subject" %in% fixed)
    }
    if (!all(fixed) && !(xsbasetype(subject) %in% c("DNA", "RNA")))
        stop("'", argname, "' value only supported for a DNA or RNA subject ",
             "(you can only use 'fixed=TRUE' with your subject)")
    fixed
}

normargCollapse <- function(collapse)
{
    if (identical(collapse, FALSE))
        return(0L)
    if (!isSingleNumber(collapse) || !(collapse %in% 0:2))
        stop("'collapse' must be FALSE, 1 or 2")
    as.integer(collapse)
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
        subject <- as(subject, "XString")
    pattern <- normargPattern(pattern, subject)
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
    fixed <- normargFixed(fixed, subject)
    .Call("XString_match_pattern_at",
          pattern, subject, at, at.type,
          max.mismatch, with.indels, fixed, ans.type,
          PACKAGE="Biostrings")
}

.vmatchPatternAt <- function(pattern, subject, at, at.type,
                             max.mismatch, with.indels, fixed, ans.type)
{
    if (!is(subject, "XStringSet"))
        subject <- as(subject, "XStringSet")
    pattern <- normargPattern(pattern, subject)
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
    fixed <- normargFixed(fixed, subject)
    .Call("XStringSet_vmatch_pattern_at",
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
### hasLetterAt()
###

hasLetterAt <- function(x, letter, at, fixed=TRUE)
{
    if (!is(x, "XStringSet")) {
        if (!is.character(x) && !is(x, "XString"))
            stop("'x' must be a character vector, or an XString or XStringSet object")
        x <- as(x, "XStringSet")
    }
    if (!is(letter, "XString")) {
        if (!isSingleString(letter))
            stop("'letter' must be a character string or an XString object")
        letter <- XString(xsbasetype(x), letter)
    } else {
        if (xsbasetype(letter) != xsbasetype(x))
            stop("'x' and 'letter' must have the same XString base type")
    }
    if (!is.numeric(at))
        stop("'at' must be a vector of integers")
    if (length(at) != length(letter))
        stop("'letter' and 'at' must have the same length")
    if (!is.integer(at))
        at <- as.integer(at)
    if (any(is.na(at)))
        stop("'at' cannot contain NAs")
    fixed <- normargFixed(fixed, x)

    .hasLetterAt1 <- function(x, l1, at1)
    {
        ans <- .Call("XStringSet_vmatch_pattern_at",
                     l1, x, at1, 0L,
                     0L, FALSE, fixed, 0L,
                     PACKAGE="Biostrings")
        ans[at1 < 1 | at1 > width(x)] <- NA
        ans
    }
    sapply(seq_len(length(letter)),
           function(i)
               .hasLetterAt1(x, subseq(letter, start=i, width=1L), at[i]))
}


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

