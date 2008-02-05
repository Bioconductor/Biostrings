### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions (not exported) used by matching functions from other
### files (like matchPattern(), matchPDict(), etc...) to check and normalize
### their arguments.
###

normalize.max.mismatch <- function(max.mismatch)
{
    if (!isSingleNumber(max.mismatch))
        stop("'max.mismatch' must be a single integer")
    max.mismatch <- as.integer(max.mismatch)
    if (max.mismatch < 0)
        stop("'max.mismatch' must be a non-negative integer")
    max.mismatch
}

### Return a logical vector of length 2.
normalize.fixed <- function(fixed, subjectClass)
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
### The "isMatching" generic and methods.
###
### Return a logical vector of the same length as 'start'.
###

.isMatching <- function(pattern, subject, start, max.mismatch, fixed)
{
    if (!is(subject, "BString"))
        subject <- BString(subject)
    if (class(pattern) != class(subject))
        pattern <- mkBString(class(subject), pattern)
    if (!is.numeric(start))
        stop("'start' must be a vector of integers")
    if (!is.integer(start))
        start <- as.integer(start)
    max.mismatch <- normalize.max.mismatch(max.mismatch)
    fixed <- normalize.fixed(fixed, class(subject))
    .Call("is_matching", pattern, subject, start,
          max.mismatch, fixed,
          PACKAGE="Biostrings")
}

setGeneric("isMatching", signature="subject",
    function(pattern, subject, start=1, max.mismatch=0, fixed=TRUE)
        standardGeneric("isMatching")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("isMatching", "character",
    function(pattern, subject, start=1, max.mismatch=0, fixed=TRUE)
    {
        .isMatching(pattern, subject, start, max.mismatch, fixed)
    }
)

### Dispatch on 'subject' (see signature of generic).
setMethod("isMatching", "BString",
    function(pattern, subject, start=1, max.mismatch=0, fixed=TRUE)
    {
        .isMatching(pattern, subject, start, max.mismatch, fixed)
    }
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
            l <- BString.substr(pattern, i, i)
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
           function(i) .bsMismatch(pattern, x@subject, x@views$start[i], fixed))
}

setGeneric("mismatch", signature="x",
    function(pattern, x, fixed=TRUE) standardGeneric("mismatch")
)

### Typical use:
###   mp <- matchPattern("TGA", DNAString("GTGACGTGCAT"), max.mismatch=2)
###   mismatch("TGA", mp)
### Dispatch on 'x' (see signature of generic).
setMethod("mismatch", "BStringViews",
    function(pattern, x, fixed)
    {
        if (class(pattern) != class(x@subject))
            pattern <- mkBString(class(x@subject), pattern)
        .mismatch(pattern, x, fixed)
    }
)

