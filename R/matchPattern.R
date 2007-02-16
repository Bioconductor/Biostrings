### =========================================================================
### The matchPattern() generic & related functions
### -------------------------------------------------------------------------


.Clongint.nbits <- function()
{
    .Call("bits_per_long", PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "gregexpr" and "gregexpr2" methods (undocumented)

### This method can miss matches (see below why)
.match.gregexpr <- function(pattern, subject, count.only)
{
    matches <- gregexpr(pattern, subject, fixed=TRUE)[[1]]
    if (length(matches) == 1 && matches == -1)
        matches <- integer(0)
    else
        attr(matches, "match.length") <- NULL
    if (count.only)
        return(length(matches))
    matches
}

### Standard R gregexpr() misses matches when they are overlapping:
###   > gregexpr("aa", c("XaaaYaa", "a"), fixed=TRUE)
###   [[1]]
###   [1] 2 6
###   attr(,"match.length")
###   [1] 2 2
###
###   [[2]]
###   [1] -1
###   attr(,"match.length")
###   [1] -1
###
### gregexpr2() is a modified version of gregexpr() that returns _all_
### matches but it only works in 'fixed=TRUE' mode (i.e. for exact matching,
### no regular expression):
###   > gregexpr2("aa", c("XaaaYaa", "a"))
###   [[1]]
###   [1] 2 3 6
###
###   [[2]]
###   [1] -1
###
### Note that, unlike gregexpr(), gregexpr2() doesn't attach a "match.length"
### attribute to each element of the returned list because, since it only works
### in 'fixed=TRUE' mode, then all the matches have the length of the pattern.
### Another difference with gregexpr() is that with gregexpr2(), the 'pattern'
### argument must be a single (non-NA, non-empty) string.

gregexpr2 <- function(pattern, text)
{
    if (!is.character(pattern) || length(pattern) != 1
      || is.na(pattern) || nchar(pattern) == 0)
        stop("invalid pattern")
    matches <- gregexpr(pattern, text, fixed=TRUE)
    nP <- nchar(pattern)
    for (i in 1:length(text)) {
        mi <- matches[[i]]
        if (length(mi) == 1 && mi == -1) {
            attr(matches[[i]], "match.length") <- NULL
        } else {
            subtexts <- substring(text[i], mi + 1, mi + 2*nP - 2)
            missing_matches <- gregexpr2(pattern, subtexts)
            for (j in 1:length(mi)) {
                mj <- missing_matches[[j]]
                if (length(mj) != 1 || mj != -1)
                    matches[[i]] <- c(matches[[i]], mi[j] + mj)
            }
            matches[[i]] <- sort(matches[[i]])
        }
    }
    matches
}

.match.gregexpr2 <- function(pattern, subject, count.only)
{
    matches <- gregexpr2(pattern, subject)[[1]]
    if (length(matches) == 1 && matches == -1)
        matches <- integer(0)
    if (count.only)
        return(length(matches))
    matches
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "naive" method

debug_naive <- function()
{
    invisible(.Call("match_naive_debug", PACKAGE="Biostrings"))
}

### Must return an integer vector.
.match.naive <- function(pattern, subject, count.only)
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
.match.boyermoore <- function(pattern, subject, count.only)
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
.match.forwardsearch <- function(pattern, subject, fixed, count.only)
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
.match.shiftor <- function(pattern, subject, mismatch, fixed, count.only)
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
    algo <- match.arg(algorithm, c("auto", "gregexpr", "gregexpr2",
                                   "naive", "boyer-moore", "forward-search", "shift-or"))
    if (algo == "gregexpr" || algo == "gregexpr2") {
        if (!is.character(subject))
            stop("'algorithms \"gregexpr\" and \"gregexpr2\" are only ",
                 "supported for character strings")
        if (length(subject) != 1 || is.na(subject) || nchar(subject) == 0)
            stop("'subject' must be a single (non-NA) string")
        if (!is.character(pattern) || length(pattern) != 1
         || is.na(pattern) || nchar(pattern) == 0)
            stop("'pattern' must be a single (non-NA, non-empty) string")
    } else {
        if (is.character(subject))
            subject <- BString(subject)
        if (class(pattern) != class(subject))
            pattern <- new(class(subject), pattern)
        if (nchar(pattern) > 10000)
            stop("patterns with more than 10000 letters are not supported, sorry")
    }
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
    if (algo == "shift-or" && nchar(pattern) > .Clongint.nbits())
        stop("your system can only support patterns up to ",
             .Clongint.nbits(), " letters\n",
             "        when 'algo' is \"shift-or\" ",
             "or 'mismatch' is != 0 or 'fixed' is FALSE")
    ans <- switch(algo,
        "gregexpr"=.match.gregexpr(pattern, subject, count.only),
        "gregexpr2"=.match.gregexpr2(pattern, subject, count.only),
        "naive"=.match.naive(pattern, subject, count.only),
        "boyer-moore"=.match.boyermoore(pattern, subject, count.only),
        "forward-search"=.match.forwardsearch(pattern, subject, fixed, count.only),
        "shift-or"=.match.shiftor(pattern, subject, mismatch,
                                fixed, count.only)
    )
    if (count.only)
        return(ans)
    if (algo == "gregexpr" || algo == "gregexpr2")
        return(ans)
    new("BStringViews", subject, ans, ans + pattern@length - as.integer(1))
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
    if (length(x) == 0)
        return(list())
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
