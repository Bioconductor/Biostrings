### =========================================================================
### Low-level matching functions
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions (not exported) used by matching functions from other
### files (like matchPattern(), matchPDict(), etc...) to check and normalize
### their arguments.
###

normargSubject <- function(subject, argname="subject")
{
    if (is(subject, "XString") || is(subject, "XStringSet"))
        return(subject)
    if (!is.character(subject))
        stop("'", argname, "' must be a character vector, ",
             "or an XString or XStringSet object")
    if (length(subject) == 1L)
        subject <- as(subject, "XString")
    else
        subject <- as(subject, "XStringSet")
    subject
}

normargPattern <- function(pattern, subject, argname="pattern")
{
    subject_baseclass <- xsbaseclass(subject)
    if (is(pattern, "XString")) {
        if (xsbaseclass(pattern) == subject_baseclass)
            return(pattern)
    } else if (!isSingleString(pattern))
        stop("'", argname, "' must be a single string or an XString object")
    pattern <- try(XString(seqtype(subject), pattern))
    if (is(pattern, "try-error"))
        stop("could not turn '", argname, "' into a ",
             subject_baseclass, " instance")
    pattern
}

normargMaxMismatch <- function(max.mismatch, argname="max.mismatch")
{
    if (!isSingleNumber(max.mismatch))
        stop("'", argname, "' must be a single integer")
    max.mismatch <- as.integer(max.mismatch)
    if (max.mismatch < 0L)
        stop("'", argname, "' must be a non-negative integer")
    max.mismatch
}

normargMinMismatch <- function(min.mismatch, max.mismatch, argname="min.mismatch")
{
    if (!isSingleNumber(min.mismatch))
        stop("'", argname, "' must be a single integer")
    min.mismatch <- as.integer(min.mismatch)
    if (min.mismatch < 0L)
        stop("'", argname, "' must be a non-negative integer")
    if (min.mismatch > max.mismatch)
        stop("'", argname, "' must be <= 'max.mismatch'")
    min.mismatch
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
    if (!all(fixed) && !(seqtype(subject) %in% c("DNA", "RNA")))
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
### .matchPatternAt()
###
### This is the horse-power behind the neditStartingAt(), neditEndingAt(),
### isMatchingStartingAt(), isMatchingEndingAt(), which.isMatchingStartingAt()
### and which.isMatchingEndingAt() low-level matching functions defined later
### in this file.
### If 'at.type == 0' then 'at' contains starting positions, otherwise it
### contains ending positions.
###

.matchPatternAt <- function(pattern, subject, at, at.type,
                            max.mismatch, min.mismatch, with.indels, fixed,
                            ans.type, auto.reduce.pattern=FALSE)
{
    subject <- normargSubject(subject)
    pattern <- normargPattern(pattern, subject)
    if (!is.numeric(at)) {
        what <- if (at.type == 0L) "starting.at" else "ending.at"
        stop("'", what, "' must be a vector of integers")
    }
    if (!is.integer(at))
        at <- as.integer(at)

    if (auto.reduce.pattern) {
        at.length <- length(at)
        P <- nchar(pattern)
        if (at.length == 1)
            at <- rep.int(at, P)
        else if (at.length != P || length(unique(at)) > 1)
            stop("With 'auto.reduce.pattern', 'at' must be a single integer")
    }

    if (ans.type == 0L) {
        max.mismatch <- length(pattern)
    } else {
        if (!is.numeric(max.mismatch))
            stop("'max.mismatch' must be a vector of integers")
        if (!is.integer(max.mismatch))
            max.mismatch <- as.integer(max.mismatch)
        if (!is.numeric(min.mismatch))
            stop("'min.mismatch' must be a vector of integers")
        if (!is.integer(min.mismatch))
            min.mismatch <- as.integer(min.mismatch)
    }
    with.indels <- normargWithIndels(with.indels)
    fixed <- normargFixed(fixed, subject)
    if (is(subject, "XString"))
        .Call2("XString_match_pattern_at",
              pattern, subject, at, at.type,
              max.mismatch, min.mismatch, with.indels, fixed, ans.type,
              auto.reduce.pattern, PACKAGE="Biostrings")
    else
        .Call2("XStringSet_vmatch_pattern_at",
              pattern, subject, at, at.type,
              max.mismatch, min.mismatch, with.indels, fixed, ans.type,
              auto.reduce.pattern, PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### neditStartingAt(), neditEndingAt() and neditAt().
###
### 'starting.at' (or 'ending.at') must be integer vectors containing the
### starting (or ending) positions of the pattern relatively to the subject.
### These functions return an integer vector of the same length as
### 'starting.at' (or 'ending.at').
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

### Dispatch on 'subject' (see signature of generic).

setMethod("neditStartingAt", "character",
    function(pattern, subject, starting.at=1, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, starting.at, 0L,
                        NA, NA, with.indels, fixed, 0L)
)

setMethod("neditStartingAt", "XString",
    function(pattern, subject, starting.at=1, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, starting.at, 0L,
                        NA, NA, with.indels, fixed, 0L)
)

setMethod("neditStartingAt", "XStringSet",
    function(pattern, subject, starting.at=1, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, starting.at, 0L,
                        NA, NA, with.indels, fixed, 0L)
)

setMethod("neditEndingAt", "character",
    function(pattern, subject, ending.at=1, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, ending.at, 1L,
                        NA, NA, with.indels, fixed, 0L)
)

setMethod("neditEndingAt", "XString",
    function(pattern, subject, ending.at=1, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, ending.at, 1L,
                        NA, NA, with.indels, fixed, 0L)
)

setMethod("neditEndingAt", "XStringSet",
    function(pattern, subject, ending.at=1, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, ending.at, 1L,
                        NA, NA, with.indels, fixed, 0L)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isMatchingStartingAt(), isMatchingEndingAt() and isMatchingAt().
###
### 'starting.at' (or 'ending.at') must be integer vectors containing the
### starting (or ending) positions of the pattern relatively to the subject.
### These functions return a logical vector of the same length as
### 'starting.at' (or 'ending.at').
###

setGeneric("isMatchingStartingAt", signature="subject",
    function(pattern, subject, starting.at=1,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        standardGeneric("isMatchingStartingAt")
)

setGeneric("isMatchingEndingAt", signature="subject",
    function(pattern, subject, ending.at=1,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        standardGeneric("isMatchingEndingAt")
)

isMatchingAt <- function(pattern, subject, at=1,
                         max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
{
    if (!is.numeric(at))
        stop("'at' must be a vector of integers")
    isMatchingStartingAt(pattern, subject, starting.at=at,
                         max.mismatch=max.mismatch, min.mismatch=min.mismatch,
                         with.indels=with.indels, fixed=fixed)
}

### Dispatch on 'subject' (see signature of generic).

setMethod("isMatchingStartingAt", "character",
    function(pattern, subject, starting.at=1,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, starting.at, 0L,
                        max.mismatch, min.mismatch, with.indels, fixed, 1L)
)

setMethod("isMatchingStartingAt", "XString",
    function(pattern, subject, starting.at=1,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, starting.at, 0L,
                        max.mismatch, min.mismatch, with.indels, fixed, 1L)
)

setMethod("isMatchingStartingAt", "XStringSet",
    function(pattern, subject, starting.at=1,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, starting.at, 0L,
                        max.mismatch, min.mismatch, with.indels, fixed, 1L)
)

setMethod("isMatchingEndingAt", "character",
    function(pattern, subject, ending.at=1,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, ending.at, 1L,
                        max.mismatch, min.mismatch, with.indels, fixed, 1L)
)

setMethod("isMatchingEndingAt", "XString",
    function(pattern, subject, ending.at=1,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, ending.at, 1L,
                        max.mismatch, min.mismatch, with.indels, fixed, 1L)
)

setMethod("isMatchingEndingAt", "XStringSet",
    function(pattern, subject, ending.at=1,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .matchPatternAt(pattern, subject, ending.at, 1L,
                        max.mismatch, min.mismatch, with.indels, fixed, 1L)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### which.isMatchingStartingAt(), which.isMatchingEndingAt() and
### which.isMatchingAt().
###
### 'starting.at' (or 'ending.at') must be integer vectors containing the
### starting (or ending) positions of the pattern relatively to the subject.
### These functions return the lowest *index* in 'starting.at' (or 'ending.at')
### for which a match occurred (or NA if no match occurred).
###

setGeneric("which.isMatchingStartingAt", signature="subject",
    function(pattern, subject, starting.at=1,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             follow.index=FALSE, auto.reduce.pattern=FALSE)
        standardGeneric("which.isMatchingStartingAt")
)

setGeneric("which.isMatchingEndingAt", signature="subject",
    function(pattern, subject, ending.at=1,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             follow.index=FALSE, auto.reduce.pattern=FALSE)
        standardGeneric("which.isMatchingEndingAt")
)

which.isMatchingAt <- function(pattern, subject, at=1,
                               max.mismatch=0, min.mismatch=0,
                               with.indels=FALSE, fixed=TRUE,
                               follow.index=FALSE, auto.reduce.pattern=FALSE)
{
    if (!is.numeric(at))
        stop("'at' must be a vector of integers")
    which.isMatchingStartingAt(pattern, subject, starting.at=at,
                               max.mismatch=max.mismatch, min.mismatch=min.mismatch,
                               with.indels=with.indels, fixed=fixed,
                               follow.index=follow.index,
                               auto.reduce.pattern=auto.reduce.pattern)
}

.to.ans.type <- function(follow.index)
{
    if (!isTRUEorFALSE(follow.index))
        stop("'follow.index' must be TRUE or FALSE")
    if (follow.index)
        return(3L)
    return(2L)
}

### Dispatch on 'subject' (see signature of generic).

setMethod("which.isMatchingStartingAt", "character",
    function(pattern, subject, starting.at=1,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             follow.index=FALSE, auto.reduce.pattern=FALSE)
        .matchPatternAt(pattern, subject, starting.at, 0L,
                        max.mismatch, min.mismatch, with.indels, fixed,
                        .to.ans.type(follow.index), auto.reduce.pattern)
)

setMethod("which.isMatchingStartingAt", "XString",
    function(pattern, subject, starting.at=1,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             follow.index=FALSE, auto.reduce.pattern=FALSE)
        .matchPatternAt(pattern, subject, starting.at, 0L,
                        max.mismatch, min.mismatch, with.indels, fixed,
                        .to.ans.type(follow.index), auto.reduce.pattern)
)

setMethod("which.isMatchingStartingAt", "XStringSet",
    function(pattern, subject, starting.at=1,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             follow.index=FALSE, auto.reduce.pattern=FALSE)
        .matchPatternAt(pattern, subject, starting.at, 0L,
                        max.mismatch, min.mismatch, with.indels, fixed,
                        .to.ans.type(follow.index), auto.reduce.pattern)
)

setMethod("which.isMatchingEndingAt", "character",
    function(pattern, subject, ending.at=1,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             follow.index=FALSE, auto.reduce.pattern=FALSE)
        .matchPatternAt(pattern, subject, ending.at, 1L,
                        max.mismatch, min.mismatch, with.indels, fixed,
                        .to.ans.type(follow.index), auto.reduce.pattern)
)

setMethod("which.isMatchingEndingAt", "XString",
    function(pattern, subject, ending.at=1,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             follow.index=FALSE, auto.reduce.pattern=FALSE)
        .matchPatternAt(pattern, subject, ending.at, 1L,
                        max.mismatch, min.mismatch, with.indels, fixed,
                        .to.ans.type(follow.index), auto.reduce.pattern)
)

setMethod("which.isMatchingEndingAt", "XStringSet",
    function(pattern, subject, ending.at=1,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             follow.index=FALSE, auto.reduce.pattern=FALSE)
        .matchPatternAt(pattern, subject, ending.at, 1L,
                        max.mismatch, min.mismatch, with.indels, fixed,
                        .to.ans.type(follow.index), auto.reduce.pattern)
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
        letter <- XString(seqtype(x), letter)
    } else {
        if (seqtype(letter) != seqtype(x))
            stop("'x' and 'letter' must have the same sequence type")
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
        ans <- .Call2("XStringSet_vmatch_pattern_at",
                     l1, x, at1, 0L, 0L, 0L, FALSE, fixed, 1L, FALSE,
                     PACKAGE="Biostrings")
        ans[at1 < 1 | at1 > width(x)] <- NA
        ans
    }
    sapply(seq_len(length(letter)),
           function(i)
               .hasLetterAt1(x, subseq(letter, start=i, width=1L), at[i]))
}

