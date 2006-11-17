### =========================================================================
### The matchPattern() generic & related functions
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Boyer-Moore

### Must return an integer vector.
BoyerMoore <- function(pattern, subject, fixed, count.only)
{
    stop("\"boyer-moore\" algorithm will be back soon...")

}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### forward-search

### Must return an integer vector.
ForwardSearch <- function(pattern, subject, fixed, count.only)
{
    stop("\"forward-search\" algorithm will be back soon...")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### shift-or

debug_shiftor <- function()
{
    invisible(.Call("shiftor_debug", PACKAGE="Biostrings"))
}

### Must return an integer vector.
ShiftOr <- function(pattern, subject, mismatch, fixed, count.only)
{
    ## We treat the edge-cases at the R level
    p <- length(pattern)
    if (p <= mismatch) {
        if (count.only)
            return(length(subject) + p - as.integer(1))
        return((1-p):(length(subject)-1))
    }
    if (p > mismatch + length(subject)) {
        if (count.only)
            return(as.integer(0))
        return(integer(0))
    }
    .Call("shiftor",
          pattern@data@xp, pattern@offset, pattern@length,
          subject@data@xp, subject@offset, subject@length,
          mismatch, fixed, count.only,
          PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Dispatch function & user interface

.matchPattern <- function(pattern, subject, algorithm, mismatch,
                          fixed, count.only=FALSE)
{
    if (class(pattern) != class(subject))
        pattern <- new(class(subject), pattern)
    algo <- match.arg(algorithm, c("boyer-moore", "forward-search", "shift-or"))
    if (length(mismatch) != 1)
        stop("'mismatch' must be a single integer")
    if (mismatch < 0)
        stop("'mismatch' must be a non-negative integer")
    if (mismatch > 0 && algo != "shift-or")
        stop("only \"shift-or\" algorithm supports 'mismatch > 0'")
    if (!is.logical(fixed) || length(fixed) != 1)
        stop("'fixed' must be a single logical")
    if (is.na(fixed)) {
        if (class(subject) == "BString")
            fixed <- TRUE
        else
            fixed <- FALSE
    }
    if (length(count.only) != 1)
        stop("'count.only' must be a single logical")
    count.only <- as.logical(count.only)
    ans <- switch(algo,
        "boyer-moore"=BoyerMoore(pattern, subject, fixed, count.only),
        "forward-search"=ForwardSearch(pattern, subject, fixed, count.only),
        "shift-or"=ShiftOr(pattern, subject, as.integer(mismatch),
                           fixed, count.only)
    )
    if (count.only)
        return(ans)
    new("BStringViews", subject, ans + as.integer(1), ans + pattern@length)
}

### Typical use:
###   matchPattern("TG", DNAString("GTGACGTGCAT"))
###   matchPattern("TG", DNAString("GTGACGTGCAT"), algo="shift", mis=1)
### Edge cases:
###   matchPattern("---", DNAString("ACGTGCA"), mismatch=3)
###   matchPattern("---", DNAString("A"))
setGeneric(
    "matchPattern",
    function(pattern, subject, algorithm="shift-or", mismatch=0, fixed=NA)
        standardGeneric("matchPattern")
)
setMethod("matchPattern", signature(subject="character"),
    function(pattern, subject, algorithm, mismatch, fixed)
    {
        subject <- BString(subject)
        .matchPattern(pattern, subject, algorithm, mismatch, fixed)
    }
)
setMethod("matchPattern", signature(subject="BString"),
    function(pattern, subject, algorithm, mismatch, fixed)
    {
        .matchPattern(pattern, subject, algorithm, mismatch, fixed)
    }
)
setMethod("matchPattern", signature(subject="BStringViews"),
    function(pattern, subject, algorithm, mismatch, fixed)
    {
        if (length(subject) != 1)
            stop("'subject' must have a single view")
        .matchPattern(pattern, subject[[1]], algorithm, mismatch, fixed)
    }
)

matchDNAPattern <- function(...)
{
    stop("matchDNAPattern() is DEPRECATED, please use matchPattern() instead")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### countPattern() is a slightly faster equivalent to length(matchPattern())

### Typical use:
###   countPattern("TG", DNAString("GTGACGTGCAT"))
###   countPattern("TG", DNAString("GTGACGTGCAT"), algo="shift", mis=1)
### Edge cases:
###   countPattern("---", DNAString("ACGTGCA"), mismatch=3)
###   countPattern("---", DNAString("A"))
setGeneric(
    "countPattern",
    function(pattern, subject, algorithm="shift-or", mismatch=0, fixed=NA)
        standardGeneric("countPattern")
)
setMethod("countPattern", signature(subject="character"),
    function(pattern, subject, algorithm, mismatch, fixed)
    {
        subject <- BString(subject)
        .matchPattern(pattern, subject, algorithm, mismatch, fixed, TRUE)
    }
)
setMethod("countPattern", signature(subject="BString"),
    function(pattern, subject, algorithm, mismatch, fixed)
    {
        .matchPattern(pattern, subject, algorithm, mismatch, fixed, TRUE)
    }
)
setMethod("countPattern", signature(subject="BStringViews"),
    function(pattern, subject, algorithm, mismatch, fixed)
    {
        if (length(subject) != 1)
            stop("'subject' must have a single view")
        .matchPattern(pattern, subject[[1]], algorithm, mismatch, fixed, TRUE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mismatch()

### Helper function used by .mismatch()
### Returns a vector of the positions of mismatches of 'pattern'
### in a view on 'subject' starting at 'start' and whose width is length(pattern).
bsMismatch <- function(pattern, subject, start, fixed)
{
    mm <- integer(0)
    j0 <- start - as.integer(1)
    for (i in 1:length(pattern)) {
        j <- j0 + i
        if (j < 1 || j > length(subject)) {
            mm <- c(mm, i)
        } else {
            l <- BString.substring(pattern, i, i)
            r <- BString.substring(subject, j, j)
            cp <- countPattern(l, r, algo="shift-or", mismatch=0, fixed)
            if (cp == 0)
                mm <- c(mm, i)
        }
    }
    mm
}

.mismatch <- function(pattern, x, fixed)
{
    if (any(width(x) != length(pattern)))
        warning("views in 'x' don't have a width equal to pattern length")
    lapply(1:length(x),
           function(i) bsMismatch(pattern, x@subject, x@start[i], fixed))
}

setGeneric(
    "mismatch",
    function(pattern, x, fixed=NA) standardGeneric("mismatch")
)

### Typical use:
###   mp <- matchPattern("TGA", DNAString("GTGACGTGCAT"), mismatch=2)
###   mismatch("TGA", mp)
setMethod("mismatch", signature(x="BStringViews"),
    function(pattern, x, fixed)
    {
        if (class(pattern) != class(x@subject)) {
            pattern <- new(class(x@subject), pattern)
        }
        .mismatch(pattern, x, fixed)
    }
)
