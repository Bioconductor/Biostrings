### =========================================================================
### Miscellaneous helper/utility functions related to string matching.
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some helper functions used by matchPattern(), matchPDict(), etc...
###

.CHARACTER.ALGOS <- c("gregexpr", "gregexpr2")

isCharacterAlgo <- function(algo) {algo %in% .CHARACTER.ALGOS}

.ALL.ALGOS <- c(
    "auto",
    "naive-exact",
    "naive-inexact",
    "boyer-moore",
    "shift-or",
    "indels",
    .CHARACTER.ALGOS
)

normargAlgorithm <- function(algorithm)
{
    if (!isSingleString(algorithm))
        stop("'algorithm' must be a single string")
    match.arg(algorithm, .ALL.ALGOS)
}

### Return a character vector containing the valid algos (best suited first)
### for the given search criteria (the search criteria is described by the
### values of 'pattern', 'max.mismatch', 'with.indels' and 'fixed').
### Raise an error if the problem "doesn't make sense".
### All its arguments must have been normalized (thru the normarg*() functions)
### before they are passed to .valid.algos().
.valid.algos <- function(pattern, max.mismatch, min.mismatch,
                         with.indels, fixed)
{
    if (is(pattern, "XString")) {
        pattern_min_length <- pattern_max_length <- length(pattern)
    } else if (is(pattern, "XStringSet")) {
        pattern_min_length <- min(width(pattern))
        pattern_max_length <- max(width(pattern))
    } else {
        stop("'pattern' not an XString or XStringSet object")
    }
    if (pattern_min_length == 0L)
        stop("empty patterns are not supported")
    ## Some arbitrary limit just to limit the size of the dynamic buffers
    ## used by the Boyer Moore algo to preprocess the pattern.
    if (pattern_max_length > 20000L)
        stop("patterns with more than 20000 letters are not supported")
    if (max.mismatch != 0L && with.indels) {
        if (min.mismatch != 0L)
            stop("'min.mismatch' must be 0 when 'with.indels' is TRUE")
        return("indels")
    }
    algos <- character(0)
    if (max.mismatch == 0L && all(fixed)) {
        algos <- c(algos, "boyer-moore")
        if (pattern_max_length <= .Clongint.nbits())
            algos <- c(algos, "shift-or")
        algos <- c(algos, "naive-exact")
    } else {
        if (min.mismatch == 0L && fixed[1] == fixed[2]
         && pattern_max_length <= .Clongint.nbits())
            algos <- c(algos, "shift-or")
    }
    c(algos, "naive-inexact") # "naive-inexact" is universal but slow
}

selectAlgo <- function(algo, pattern, max.mismatch, min.mismatch,
                       with.indels, fixed)
{
    algos <- .valid.algos(pattern, max.mismatch, min.mismatch,
                          with.indels, fixed)
    if (algo == "auto")
        return(algos[1])
    if (!(algo %in% algos))
        stop("valid algos for your string matching problem (best suited first): ",
             paste(paste("\"", algos, "\"", sep=""), collapse=", "))
    algo
}


### =========================================================================
### Some utility functions for operating on the matches returned by a
### high-level matching function like matchPattern(), matchPDict(), etc...
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mismatch()
###

### Returns a vector of the positions of mismatches of 'pattern'
### in a view on 'subject' starting at 'start' and whose width is
### length(pattern).
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
        elementNROWS(mismatches)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "coverage" methods.
###

setMethod("coverage", "MaskedXString",
    function(x, shift=0L, width=NULL, weight=1L)
        coverage(masks(x), shift=shift, width=width, weight=weight)
)

setMethod("coverage", "MIndex",
    function(x, shift=0L, width=NULL, weight=1L)
        coverage(unlist(x), shift=shift, width=width, weight=weight)
)

