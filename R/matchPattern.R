# ===========================================================================
# The matchPattern() generic
# ---------------------------------------------------------------------------


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Boyer-Moore

# Must return an integer vector.
BoyerMoore <- function(pattern, x, fixed, count.only)
{
    stop("\"boyer-moore\" algorithm will be back soon...")

}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# forward-search

# Must return an integer vector.
ForwardSearch <- function(pattern, x, fixed, count.only)
{
    stop("\"forward-search\" algorithm will be back soon...")
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# shift-or

debug_shiftor <- function()
{
    invisible(.Call("shiftor_debug", PACKAGE="Biostrings"))
}

# Must return an integer vector.
ShiftOr <- function(pattern, x, mismatch, fixed, count.only)
{
    # We treat the edge-cases at the R level
    p <- length(pattern)
    if (p <= mismatch) {
        if (count.only)
            return(length(x) + p - as.integer(1))
        return((1-p):(length(x)-1))
    }
    if (p > mismatch + length(x)) {
        if (count.only)
            return(as.integer(0))
        return(integer(0))
    }
    .Call("shiftor",
          pattern@data@xp, pattern@offset, pattern@length,
          x@data@xp, x@offset, x@length,
          mismatch, fixed, count.only,
          PACKAGE="Biostrings")
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Dispatch function & user interface

.matchPattern <- function(pattern, x, algorithm, mismatch,
                          fixed, count.only=FALSE)
{
    algo <- match.arg(algorithm, c("boyer-moore", "forward-search", "shift-or"))
    if (length(mismatch) != 1)
        stop("'mismatch' must be a single integer")
    if (mismatch < 0)
        stop("'mismatch' must be a non-negative integer")
    if (mismatch > 0 && algo != "shift-or")
        stop("only \"shift-or\" algorithm supports 'mismatch > 0'")
    if (length(fixed) != 1)
        stop("'fixed' must be a single logical")
    fixed <- as.logical(fixed)
    if (length(count.only) != 1)
        stop("'count.only' must be a single logical")
    count.only <- as.logical(count.only)
    ans <- switch(algo,
        "boyer-moore"=BoyerMoore(pattern, x, fixed, count.only),
        "forward-search"=ForwardSearch(pattern, x, fixed, count.only),
        "shift-or"=ShiftOr(pattern, x, as.integer(mismatch), fixed, count.only)
    )
    if (count.only)
        return(ans)
    BStringViews(x, ans + as.integer(1), ans + length(pattern))
}

setGeneric(
    "matchPattern",
    function(pattern, x, algorithm="shift-or", mismatch=0, fixed=FALSE)
        standardGeneric("matchPattern")
)

# Typical use:
#   matchPattern("TG", DNAString("GTGACGTGCAT"))
#   matchPattern("TG", DNAString("GTGACGTGCAT"), algo="shift", mis=1)
# Edge cases:
#   matchPattern("---", DNAString("ACGTGCA"), mismatch=3)
#   matchPattern("---", DNAString("A"))
setMethod("matchPattern", signature(x="DNAString"),
    function(pattern, x, algorithm, mismatch, fixed)
    {
        if (class(pattern) != "DNAString")
            pattern <- DNAString(pattern)
        .matchPattern(pattern, x, algorithm, mismatch, fixed)
    }
)

setMethod("matchPattern", signature(x="RNAString"),
    function(pattern, x, algorithm, mismatch, fixed)
    {
        if (class(pattern) != "RNAString")
            pattern <- RNAString(pattern)
        .matchPattern(pattern, x, algorithm, mismatch, fixed)
    }
)

matchDNAPattern <- function(...)
{
    cat("This function has been deprecated.\n")
    cat("Please use matchPattern() instead.\n")
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# countPattern() is a fast equivalent to length(matchPattern())

setGeneric(
    "countPattern",
    function(pattern, x, algorithm="shift-or", mismatch=0, fixed=FALSE)
        standardGeneric("countPattern")
)

# Typical use:
#   countPattern("TG", DNAString("GTGACGTGCAT"))
#   countPattern("TG", DNAString("GTGACGTGCAT"), algo="shift", mis=1)
# Edge cases:
#   countPattern("---", DNAString("ACGTGCA"), mismatch=3)
#   countPattern("---", DNAString("A"))
setMethod("countPattern", signature(x="DNAString"),
    function(pattern, x, algorithm, mismatch, fixed)
    {
        if (class(pattern) != "DNAString")
            pattern <- DNAString(pattern)
        .matchPattern(pattern, x, algorithm, mismatch, fixed, TRUE)
    }
)

setMethod("countPattern", signature(x="RNAString"),
    function(pattern, x, algorithm, mismatch, fixed)
    {
        if (class(pattern) != "RNAString")
            pattern <- RNAString(pattern)
        .matchPattern(pattern, x, algorithm, mismatch, fixed, TRUE)
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mismatch()

# Helper function used by .mismatch()
# Returns a vector of the positions of mismatches of 'pattern'
# in a view on 'x' starting at 'first' and whose width is length(pattern).
bsMismatch <- function(pattern, x, first, fixed)
{
    mm <- integer(0)
    j0 <- first - as.integer(1)
    for (i in 1:length(pattern)) {
        j <- j0 + i
        if (j < 1 || j > length(x)) {
            mm <- c(mm, i)
        } else {
            l <- bsSubstr(pattern, i, i)
            r <- bsSubstr(x, j, j)
            cp <- countPattern(l, r, algo="shift-or", mismatch=0, fixed)
            if (cp == 0)
                mm <- c(mm, i)
        }
    }
    mm
}

.mismatch <- function(pattern, x, fixed)
{
    if (fixed)
        stop("'fixed=TRUE' not supported yet")
    if (any(width(x) != length(pattern)))
        warning("views in 'x' don't have a width equal to pattern length")
    lapply(1:length(x),
           function(i) bsMismatch(pattern, x@subject, x@first[i], fixed))
}

setGeneric(
    "mismatch",
    function(pattern, x, fixed=FALSE) standardGeneric("mismatch")
)

# Typical use:
#   mp <- matchPattern("TGA", DNAString("GTGACGTGCAT"), mismatch=2)
#   mismatch("TGA", mp)
setMethod("mismatch", signature(x="BStringViews"),
    function(pattern, x, fixed)
    {
        if (class(pattern) != class(x@subject)) {
            pattern <- new(class(x@subject), pattern)
        }
        .mismatch(pattern, x, fixed)
    }
)
