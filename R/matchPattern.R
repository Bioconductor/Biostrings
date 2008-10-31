### =========================================================================
### The matchPattern() generic & related functions
### -------------------------------------------------------------------------


.Clongint.nbits <- function()
{
    .Call("bits_per_long", PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### matchPattern algos for standard character vectors.
###
### Note that these algos are not documented.
###

.CHARACTER.ALGOS <- c("gregexpr", "gregexpr2")

.is.character.algo <- function(algo)
{
    algo %in% .CHARACTER.ALGOS
}

### This matchPattern algo can miss matches (see below why).
.matchPattern.gregexpr <- function(pattern, subject)
{
    matches <- gregexpr(pattern, subject, fixed=TRUE)[[1]]
    if (length(matches) == 1 && matches == -1)
        matches <- integer(0)
    else
        attr(matches, "match.length") <- NULL
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

.matchPattern.gregexpr2 <- function(pattern, subject)
{
    matches <- gregexpr2(pattern, subject)[[1]]
    if (length(matches) == 1 && matches == -1)
        matches <- integer(0)
    matches
}

.character.matchPattern <- function(pattern, subject, algo,
                                    max.mismatch, fixed, count.only)
{
    if (!isSingleString(pattern) || nchar(pattern) == 0)
        stop("'pattern' must be a single (non-empty) string ",
             "for this algorithm")
    if (!isSingleString(subject) || nchar(subject) == 0)
        stop("'subject' must be a single (non-empty) string ",
             "for this algorithm")
    max.mismatch <- normargMaxMismatch(max.mismatch)
    ## We need to cheat on normargFixed()
    fixed <- normargFixed(fixed, "DNAString")
    if (!(max.mismatch == 0 && all(fixed)))
        stop("this algorithm only supports exact matching ",
             "(i.e. 'max.mismatch=0' and 'fixed=TRUE')")
    if (!isTRUEorFALSE(count.only))
        stop("'count.only' must be TRUE or FALSE")
    matches <- switch(algo,
                      "gregexpr"=.matchPattern.gregexpr(pattern, subject),
                      "gregexpr2"=.matchPattern.gregexpr2(pattern, subject))
    if (count.only) length(matches) else matches
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .XString.matchPattern() and .XStringViews.matchPattern()
###

.ALL.ALGOS <- c(
    "auto",
    "naive-exact",
    "naive-inexact",
    "boyer-moore",
    "shift-or",
    .CHARACTER.ALGOS
)

.normargAlgorithm <- function(algorithm)
{
    if (!isSingleString(algorithm))
        stop("'algorithm' must be a single string")
    match.arg(algorithm, .ALL.ALGOS)
}

### Return a character vector containing the valid algos (best suited first)
### for the given problem (problem is described by the values of 'pattern',
### 'max.mismatch' and 'fixed').
### Raise an error if the problem "doesn't make sense".
### Make sure that:
###   1. 'pattern' is of the same class as 'subject'
###   2. the 'max.mismatch' ans 'fixed' args have been normalized
### before you call .valid.algos()
.valid.algos <- function(pattern, max.mismatch, fixed)
{
    if (length(pattern) == 0)
        stop("empty pattern")
    if (length(pattern) > 20000)
        stop("patterns with more than 20000 letters are not supported, sorry")
    algos <- character(0)
    if (max.mismatch == 0 && all(fixed)) {
        algos <- c(algos, "boyer-moore")
        if (length(pattern) <= .Clongint.nbits())
            algos <- c(algos, "shift-or")
        algos <- c(algos, "naive-exact")
    } else {
        if (fixed[1] == fixed[2] && length(pattern) <= .Clongint.nbits())
            algos <- c(algos, "shift-or")
    }
    c(algos, "naive-inexact") # "naive-inexact" is universal but slow
}

.select.algo <- function(algo, pattern, max.mismatch, fixed)
{
    algos <- .valid.algos(pattern, max.mismatch, fixed)
    if (algo == "auto")
        return(algos[1])
    if (!(algo %in% algos))
        stop("valid algos for your problem (best suited first): ",
             paste(paste("\"", algos, "\"", sep=""), collapse=", "))
    algo
}

.XString.matchPattern <- function(pattern, subject, algorithm,
                                  max.mismatch, fixed, count.only=FALSE)
{
    algo <- .normargAlgorithm(algorithm)
    if (.is.character.algo(algo))
        return(.character.matchPattern(pattern, subject, algo,
                                       max.mismatch, fixed, count.only))
    if (!is(subject, "XString"))
        subject <- XString(NULL, subject)
    if (class(pattern) != class(subject))
        pattern <- XString(class(subject), pattern)
    max.mismatch <- normargMaxMismatch(max.mismatch)
    fixed <- normargFixed(fixed, class(subject))
    if (!isTRUEorFALSE(count.only))
        stop("'count.only' must be TRUE or FALSE")
    algo <- .select.algo(algo, pattern, max.mismatch, fixed)
    matches <- .Call("XString_match_pattern",
                     pattern, subject, algo,
                     max.mismatch, fixed,
                     count.only,
                     PACKAGE="Biostrings")
    if (count.only)
        return(matches)
    ans_width <- rep.int(length(pattern), length(matches))
    unsafe.newXStringViews(subject, matches, ans_width)
}

.XStringViews.matchPattern <- function(pattern, subject, algorithm,
                                       max.mismatch, fixed, count.only=FALSE)
{
    algo <- .normargAlgorithm(algorithm)
    if (.is.character.algo(algo))
        stop("'subject' must be a single (non-empty) string ",
             "for this algorithm")
    if (class(pattern) != class(subject(subject)))
        pattern <- XString(class(subject(subject)), pattern)
    max.mismatch <- normargMaxMismatch(max.mismatch)
    fixed <- normargFixed(fixed, class(subject(subject)))
    if (!isTRUEorFALSE(count.only))
        stop("'count.only' must be TRUE or FALSE")
    algo <- .select.algo(algo, pattern, max.mismatch, fixed)
    matches <- .Call("XStringViews_match_pattern",
                     pattern,
                     subject(subject), start(subject), width(subject),
                     algo, max.mismatch, fixed,
                     count.only,
                     PACKAGE="Biostrings")
    if (count.only)
        return(matches)
    ans_width <- rep.int(length(pattern), length(matches))
    unsafe.newXStringViews(subject(subject), matches, ans_width)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "matchPattern" generic and methods.
###
### Typical use:
###   matchPattern("TG", DNAString("GTGACGTGCAT"))
###   matchPattern("TG", DNAString("GTGACGTGCAT"), algo="shift", max.mismatch=1)
### Edge cases:
###   matchPattern("---", DNAString("ACGTGCA"), max.mismatch=3)
###   matchPattern("---", DNAString("A"))
###

setGeneric("matchPattern", signature="subject",
    function(pattern, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)
        standardGeneric("matchPattern")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPattern", "character",
    function(pattern, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)
        .XString.matchPattern(pattern, subject, algorithm, max.mismatch, fixed)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPattern", "XString",
    function(pattern, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)
        .XString.matchPattern(pattern, subject, algorithm, max.mismatch, fixed)
)

### Dispatch on 'subject' (see signature of generic).
### WARNING: Unlike the other "matchPattern" methods, the XStringViews object
### returned by this method is not guaranteed to have its views ordered from
### left to right in general! One important particular case where this is
### guaranteed though is when 'isNormal(subject)' is TRUE (i.e. 'subject' is
### a normal XStringViews object) and 'max.mismatch=0' (no "out of limits"
### matches).
setMethod("matchPattern", "XStringViews",
    function(pattern, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)
        .XStringViews.matchPattern(pattern, subject, algorithm,
                                   max.mismatch, fixed)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPattern", "MaskedXString",
    function(pattern, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)
        matchPattern(pattern, toXStringViewsOrXString(subject),
                     algorithm=algorithm, max.mismatch=max.mismatch, fixed=fixed)
)

matchDNAPattern <- function(...) .Defunct("matchPattern")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "countPattern" generic and methods.
###
### countPattern() is equivalent to length(matchPattern()) but should be
### slightly faster, especially when there is a high number of matches.
###
### Typical use:
###   countPattern("TG", DNAString("GTGACGTGCAT"))
###   countPattern("TG", DNAString("GTGACGTGCAT"), max.mismatch=1)
### Edge cases:
###   countPattern("---", DNAString("ACGTGCA"), max.mismatch=3)
###   countPattern("---", DNAString("A"))
###

setGeneric("countPattern", signature="subject",
    function(pattern, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)
        standardGeneric("countPattern")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "character",
    function(pattern, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)
        .XString.matchPattern(pattern, subject, algorithm,
                              max.mismatch, fixed, count.only=TRUE)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "XString",
    function(pattern, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)
        .XString.matchPattern(pattern, subject, algorithm,
                              max.mismatch, fixed, count.only=TRUE)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "XStringViews",
    function(pattern, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)
        .XStringViews.matchPattern(pattern, subject, algorithm,
                                   max.mismatch, fixed, count.only=TRUE)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "MaskedXString",
    function(pattern, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)
        countPattern(pattern, toXStringViewsOrXString(subject),
                     algorithm=algorithm, max.mismatch=max.mismatch, fixed=fixed)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "vmatchPattern" and "vcountPattern" generic and methods.
###
### These are vectorized versions of matchPattern() and countPattern().
### vmatchPattern() returns an MIndex object and vcountPattern() an integer
### vector (like matchPDict() and countPDict() do).
###

.XStringSet.vmatchPattern <- function(pattern, subject, algorithm,
                                      max.mismatch, fixed,
                                      count.only=FALSE)
{
    if (!is(subject, "XStringSet"))
        subject <- XStringSet(NULL, subject)
    algo <- .normargAlgorithm(algorithm)
    if (.is.character.algo(algo)) 
        stop("'subject' must be a single (non-empty) string ", 
             "for this algorithm")
    if (class(pattern) != class(super(subject)))
        pattern <- XString(class(super(subject)), pattern)
    max.mismatch <- normargMaxMismatch(max.mismatch)
    fixed <- normargFixed(fixed, class(super(subject)))
    if (!isTRUEorFALSE(count.only)) 
        stop("'count.only' must be TRUE or FALSE")
    algo <- .select.algo(algo, pattern, max.mismatch, fixed)
    matches <- .Call("XStringSet_vmatch_pattern", pattern, subject,
                     algo, max.mismatch, fixed, count.only,
                     PACKAGE = "Biostrings")
    if (count.only)
        return(matches)
    new("ByPos_MIndex", ends=matches, width=length(pattern))
}

setGeneric("vmatchPattern", signature="subject",
    function(pattern, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)
        standardGeneric("vmatchPattern")
)

setMethod("vmatchPattern", "character",
    function(pattern, subject, algorithm="auto", max.mismatch=0L, fixed=TRUE)
        .XStringSet.vmatchPattern(pattern, subject, algorithm,
                                  max.mismatch, fixed)
)

setMethod("vmatchPattern", "XStringSet",
    function(pattern, subject, algorithm="auto", max.mismatch=0L, fixed=TRUE)
        .XStringSet.vmatchPattern(pattern, subject, algorithm,
                                  max.mismatch, fixed)
)

# TODO: Add method for XStringViews.

setGeneric("vcountPattern", signature="subject",
    function(pattern, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)
        standardGeneric("vcountPattern")
)

setMethod("vcountPattern", "character",
    function(pattern, subject, algorithm="auto", max.mismatch=0L, fixed=TRUE)
        .XStringSet.vmatchPattern(pattern, subject, algorithm,
                                  max.mismatch, fixed, count.only=TRUE)
)

setMethod("vcountPattern", "XStringSet",
    function(pattern, subject, algorithm="auto", max.mismatch=0L, fixed=TRUE)
        .XStringSet.vmatchPattern(pattern, subject, algorithm,
                                  max.mismatch, fixed, count.only=TRUE)
)

# TODO: Add method for XStringViews.

