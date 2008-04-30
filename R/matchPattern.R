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

.character.matchPattern <- function(pattern, subject, algo)
{
    switch(algo,
        "gregexpr"=.matchPattern.gregexpr(pattern, subject),
        "gregexpr2"=.matchPattern.gregexpr2(pattern, subject))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .matchPattern()
###

.normalize.algorithm <- function(algorithm)
{
    if (!isSingleString(algorithm))
        stop("'algorithm' must be a single string")
    match.arg(algorithm, c("auto", "gregexpr", "gregexpr2",
                           "naive-exact", "naive-inexact",
                           "boyer-moore", "shift-or"))
}

### Return a character vector containing the valid algos (best suited first)
### for the given problem (problem is described by the values of 'pattern',
### 'max.mismatch' and 'fixed').
### Raise an error if the problem "doesn't make sense".
### Make sure that:
###   1. 'pattern' is of the same class as 'subject'
###   2. 'max.mismatch' ans 'fixed' have been normalized
### before you call .valid.algos()
.valid.algos <- function(pattern, max.mismatch, fixed)
{
    if (nchar(pattern) == 0)
        stop("empty pattern")
    if (nchar(pattern) > 20000)
        stop("patterns with more than 20000 letters are not supported, sorry")
    algos <- character(0)
    if (max.mismatch == 0 && all(fixed)) {
        algos <- c(algos, "boyer-moore")
        if (nchar(pattern) <= .Clongint.nbits())
            algos <- c(algos, "shift-or")
        algos <- c(algos, "naive-exact")
    } else {
        if (fixed[1] == fixed[2] && nchar(pattern) <= .Clongint.nbits())
            algos <- c(algos, "shift-or")
    }
    c(algos, "naive-inexact") # "naive-inexact" is universal but slow
}

.matchPattern <- function(pattern, subject, algorithm, max.mismatch, fixed,
                          count.only=FALSE)
{
    algo <- .normalize.algorithm(algorithm)
    if (.is.character.algo(algo)) {
        if (!isSingleString(subject) || nchar(subject) == 0)
            stop("'subject' must be a single (and non-empty) string ",
                 "for this algorithm")
        if (!isSingleString(pattern) || nchar(pattern) == 0)
            stop("'pattern' must be a single (and non-empty) string ",
                 "for this algorithm")
    } else {
        if (!is(subject, "XString"))
            subject <- BString(subject)
        if (class(pattern) != class(subject))
            pattern <- XString(class(subject), pattern)
    }
    max.mismatch <- normalize.max.mismatch(max.mismatch)
    fixed <- normalize.fixed(fixed, class(subject))
    if (!isTRUEorFALSE(count.only))
        stop("'count.only' must be TRUE or FALSE")
    if (.is.character.algo(algo)) {
        if (!(max.mismatch == 0 && all(fixed)))
            stop("this algorithm only supports exact matching ",
                 "(i.e. 'max.mismatch=0' and 'fixed=TRUE')")
        matches <- .character.matchPattern(pattern, subject, algo)
        if (count.only)
            matches <- length(matches)
        return(matches)
    } 
    algos <- .valid.algos(pattern, max.mismatch, fixed)
    if (algo == "auto")
        algo <- algos[1]
    else if (!(algo %in% algos))
        stop("valid algos for your problem (best suited first): ",
             paste(paste("\"", algos, "\"", sep=""), collapse=", "))
    matches <- .Call("match_pattern",
                     pattern, subject, algo,
                     max.mismatch, fixed,
                     count.only,
                     PACKAGE="Biostrings")
    if (count.only)
        return(matches)
    ans_width <- rep.int(nchar(pattern), length(matches))
    new("XStringViews", subject,
        start=matches, width=ans_width, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "matchPattern" generic and methods.
###
### Typical use:
###   matchPattern("TG", DNAString("GTGACGTGCAT"))
###   matchPattern("TG", DNAString("GTGACGTGCAT"), algo="shift", mis=1)
### Edge cases:
###   matchPattern("---", DNAString("ACGTGCA"), max.mismatch=3)
###   matchPattern("---", DNAString("A"))

setGeneric("matchPattern", signature="subject",
    function(pattern, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)
        standardGeneric("matchPattern")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPattern", "character",
    function(pattern, subject, algorithm, max.mismatch, fixed)
    {
        .matchPattern(pattern, subject, algorithm, max.mismatch, fixed)
    }
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPattern", "XString",
    function(pattern, subject, algorithm, max.mismatch, fixed)
    {
        .matchPattern(pattern, subject, algorithm, max.mismatch, fixed)
    }
)

### Dispatch on 'subject' (see signature of generic).
### WARNING: Unlike the other "matchPattern" methods, the XStringViews object
### returned by this method is not guaranteed to have its views ordered from
### left to right in general! One important particular case where this is
### guaranteed though is when 'subject' is a normalized XStringViews object
### and 'max.mismatch=0' (no "out of limits" matches).
setMethod("matchPattern", "XStringViews",
    function(pattern, subject, algorithm, max.mismatch, fixed)
    {
        ans_start <- ans_width <- integer(0)
        for (i in seq_len(length(subject))) {
            pm <- .matchPattern(pattern, subject[[i]], algorithm, max.mismatch, fixed)
            offset <- start(subject)[i] - 1L
            ans_start <- c(ans_start, offset + start(pm))
            ans_width <- c(ans_width, width(pm))
        }
        new("XStringViews", subject(subject),
            start=ans_start, width=ans_width, check=FALSE)
    }
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPattern", "MaskedXString",
    function(pattern, subject, algorithm, max.mismatch, fixed)
        matchPattern(pattern, as(subject, "XStringViews"), algorithm, max.mismatch, fixed)
)

matchDNAPattern <- function(...) .Defunct("matchPattern")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### countPattern() is a slightly faster equivalent to length(matchPattern())

### Typical use:
###   countPattern("TG", DNAString("GTGACGTGCAT"))
###   countPattern("TG", DNAString("GTGACGTGCAT"), algo="shift", mis=1)
### Edge cases:
###   countPattern("---", DNAString("ACGTGCA"), max.mismatch=3)
###   countPattern("---", DNAString("A"))
setGeneric("countPattern", signature="subject",
    function(pattern, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)
        standardGeneric("countPattern")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "character",
    function(pattern, subject, algorithm, max.mismatch, fixed)
    {
        .matchPattern(pattern, subject, algorithm, max.mismatch, fixed, count.only=TRUE)
    }
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "XString",
    function(pattern, subject, algorithm, max.mismatch, fixed)
    {
        .matchPattern(pattern, subject, algorithm, max.mismatch, fixed, count.only=TRUE)
    }
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "XStringViews",
    function(pattern, subject, algorithm, max.mismatch, fixed)
    {
        sum(
            sapply(seq_len(length(subject)),
                   function(i) .matchPattern(pattern, subject[[i]],
                                             algorithm, max.mismatch, fixed,
                                             count.only=TRUE)
            )
        )
    }
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "MaskedXString",
    function(pattern, subject, algorithm, max.mismatch, fixed)
        countPattern(pattern, as(subject, "XStringViews"), algorithm, max.mismatch, fixed)
)

