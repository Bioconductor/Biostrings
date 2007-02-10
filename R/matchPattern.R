### =========================================================================
### The matchPattern() generic & related functions
### -------------------------------------------------------------------------


.Clongint.nbits <- function()
{
    .Call("bits_per_long", PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "naive" method

debug_naive <- function()
{
    invisible(.Call("match_naive_debug", PACKAGE="Biostrings"))
}

### Must return an integer vector.
.matchNaive <- function(pattern, subject, count.only)
{
    .Call("match_naive",
          pattern@data@xp, pattern@offset, pattern@length,
          subject@data@xp, subject@offset, subject@length,
          count.only,
          PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Boyer-Moore

debug_boyermoore <- function()
{
    invisible(.Call("match_boyermoore_debug", PACKAGE="Biostrings"))
}

### Must return an integer vector.
.matchBoyerMoore <- function(pattern, subject, count.only)
{
    ## We treat the edge-cases at the R level
    p <- length(pattern)
    if (p > length(subject)) {
        if (count.only)
            return(as.integer(0))
        return(integer(0))
    }
    .Call("match_boyermoore",
          pattern@data@xp, pattern@offset, pattern@length,
          subject@data@xp, subject@offset, subject@length,
          count.only,
          PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### forward-search

### Must return an integer vector.
.matchForwardSearch <- function(pattern, subject, fixed, count.only)
{
    stop("\"forward-search\" algorithm will be back soon...")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### shift-or

debug_shiftor <- function()
{
    invisible(.Call("match_shiftor_debug", PACKAGE="Biostrings"))
}

### Must return an integer vector.
.matchShiftOr <- function(pattern, subject, mismatch, fixed, count.only)
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
    .Call("match_shiftor",
          pattern@data@xp, pattern@offset, pattern@length,
          subject@data@xp, subject@offset, subject@length,
          mismatch, fixed, count.only,
          PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Dispatch function & user interface

.matchPattern <- function(pattern, subject, algorithm, mismatch, fixed,
                          count.only=FALSE)
{
    if (class(pattern) != class(subject))
        pattern <- new(class(subject), pattern)
    if (length(pattern) > 10000)
        stop("patterns with more than 10000 letters are not supported, sorry")
    algo <- match.arg(algorithm, c("auto", "naive", "boyer-moore", "forward-search", "shift-or"))
    if (!is.numeric(mismatch) || length(mismatch) != 1 || is.na(mismatch))
        stop("'mismatch' must be a single integer")
    mismatch <- as.integer(mismatch)
    if (mismatch < 0)
        stop("'mismatch' must be a non-negative integer")
    if (!is.logical(fixed) || length(fixed) != 1 || is.na(fixed))
        stop("'fixed' must be TRUE or FALSE")
    if (!fixed && !(class(subject) %in% c("DNAString", "RNAString")))
        stop("'fixed=FALSE' is only supported for DNAString or RNAString objects")
    if (!is.logical(count.only) || length(count.only) != 1 || is.na(count.only))
        stop("'count.only' must be TRUE or FALSE")
    if (algo == "auto") {
        ## We try to choose the best algo
        if (mismatch != 0 || !fixed) {
            algo <- "shift-or"
        } else {
            algo <- "boyer-moore"
        }
    } else {
        ## We check that the algo choosen by the user is a valid choice
        if (mismatch != 0 && algo != "shift-or")
            stop("only the \"shift-or\" algorithm supports 'mismatch != 0'")
        if (!fixed && algo != "shift-or")
            stop("only the \"shift-or\" algorithm supports 'fixed=FALSE'")
    }
    ## At this point we have the guarantee that if 'mismatch' is != 0 or 'fixed'
    ## is FALSE then 'algo' is "shift-or"
    if (algo == "shift-or" && length(pattern) > .Clongint.nbits())
        stop("your system can only support patterns up to ",
             .Clongint.nbits(), " letters\n",
             "        when 'algo' is \"shift-or\" ",
             "or 'mismatch' is != 0 or 'fixed' is FALSE")
    ans <- switch(algo,
        "naive"=.matchNaive(pattern, subject, count.only),
        "boyer-moore"=.matchBoyerMoore(pattern, subject, count.only),
        "forward-search"=.matchForwardSearch(pattern, subject, fixed, count.only),
        "shift-or"=.matchShiftOr(pattern, subject, mismatch,
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
    function(pattern, subject, algorithm="auto", mismatch=0, fixed=TRUE)
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
    function(pattern, subject, algorithm="auto", mismatch=0, fixed=TRUE)
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
    function(pattern, x, fixed=TRUE) standardGeneric("mismatch")
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
