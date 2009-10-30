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
                                    max.mismatch, fixed,
                                    count.only)
{
    if (!isSingleString(pattern) || nchar(pattern) == 0)
        stop("'pattern' must be a single (non-empty) string ",
             "for this algorithm")
    if (!isSingleString(subject) || nchar(subject) == 0)
        stop("'subject' must be a single (non-empty) string ",
             "for this algorithm")
    max.mismatch <- normargMaxMismatch(max.mismatch)
    ## we cheat on normargFixed() to keep it quiet
    fixed <- normargFixed(fixed, DNAString())
    if (!(max.mismatch == 0L && all(fixed)))
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
    "indels",
    .CHARACTER.ALGOS
)

.normargAlgorithm <- function(algorithm)
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
.valid.algos <- function(pattern, max.mismatch, min.mismatch, with.indels, fixed)
{
    if (length(pattern) == 0L)
        stop("empty pattern")
    if (length(pattern) > 20000L)
        stop("patterns with more than 20000 letters are not supported, sorry")
    if (max.mismatch != 0L && with.indels) {
        if (min.mismatch != 0L)
            stop("'min.mismatch' must be 0 when 'with.indels' is TRUE")
        return("indels")
    }
    algos <- character(0)
    if (max.mismatch == 0L && all(fixed)) {
        algos <- c(algos, "boyer-moore")
        if (length(pattern) <= .Clongint.nbits())
            algos <- c(algos, "shift-or")
        algos <- c(algos, "naive-exact")
    } else {
        if (min.mismatch == 0L && fixed[1] == fixed[2]
         && length(pattern) <= .Clongint.nbits())
            algos <- c(algos, "shift-or")
    }
    c(algos, "naive-inexact") # "naive-inexact" is universal but slow
}

.select.algo <- function(algo, pattern, max.mismatch, min.mismatch, with.indels, fixed)
{
    algos <- .valid.algos(pattern, max.mismatch, min.mismatch, with.indels, fixed)
    if (algo == "auto")
        return(algos[1])
    if (!(algo %in% algos))
        stop("valid algos for your problem (best suited first): ",
             paste(paste("\"", algos, "\"", sep=""), collapse=", "))
    algo
}

.XString.matchPattern <- function(pattern, subject, algorithm,
                                  max.mismatch, min.mismatch, with.indels, fixed,
                                  count.only=FALSE)
{
    algo <- .normargAlgorithm(algorithm)
    if (.is.character.algo(algo))
        return(.character.matchPattern(pattern, subject, algo,
                                       max.mismatch, fixed, count.only))
    if (!is(subject, "XString"))
        subject <- XString(NULL, subject)
    pattern <- normargPattern(pattern, subject)
    max.mismatch <- normargMaxMismatch(max.mismatch)
    min.mismatch <- normargMinMismatch(min.mismatch, max.mismatch)
    with.indels <- normargWithIndels(with.indels)
    fixed <- normargFixed(fixed, subject)
    if (!isTRUEorFALSE(count.only))
        stop("'count.only' must be TRUE or FALSE")
    algo <- .select.algo(algo, pattern, max.mismatch, min.mismatch, with.indels, fixed)
    C_ans <- .Call("XString_match_pattern",
                   pattern, subject, algo,
                   max.mismatch, min.mismatch, with.indels, fixed,
                   count.only,
                   PACKAGE="Biostrings")
    if (count.only)
        return(C_ans)
    unsafe.newXStringViews(subject, start(C_ans), width(C_ans))
}

.XStringViews.matchPattern <- function(pattern, subject, algorithm,
                                       max.mismatch, min.mismatch, with.indels, fixed,
                                       count.only=FALSE)
{
    algo <- .normargAlgorithm(algorithm)
    if (.is.character.algo(algo))
        stop("'subject' must be a single (non-empty) string ",
             "for this algorithm")
    pattern <- normargPattern(pattern, subject)
    max.mismatch <- normargMaxMismatch(max.mismatch)
    min.mismatch <- normargMinMismatch(min.mismatch, max.mismatch)
    with.indels <- normargWithIndels(with.indels)
    fixed <- normargFixed(fixed, subject)
    if (!isTRUEorFALSE(count.only))
        stop("'count.only' must be TRUE or FALSE")
    algo <- .select.algo(algo, pattern, max.mismatch, min.mismatch, with.indels, fixed)
    C_ans <- .Call("XStringViews_match_pattern",
                   pattern, subject(subject), start(subject), width(subject),
                   algo,
                   max.mismatch, min.mismatch, with.indels, fixed,
                   count.only,
                   PACKAGE="Biostrings")
    if (count.only)
        return(C_ans)
    unsafe.newXStringViews(subject(subject), start(C_ans), width(C_ans))
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
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        standardGeneric("matchPattern")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPattern", "character",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .XString.matchPattern(pattern, subject, algorithm,
                              max.mismatch, min.mismatch, with.indels, fixed)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPattern", "XString",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .XString.matchPattern(pattern, subject, algorithm,
                              max.mismatch, min.mismatch, with.indels, fixed)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPattern", "XStringSet",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        stop("please use vmatchPattern() when 'subject' is an XStringSet object (multiple sequence)")
)

### Dispatch on 'subject' (see signature of generic).
### WARNING: Unlike the other "matchPattern" methods, the XStringViews object
### returned by this method is not guaranteed to have its views ordered from
### left to right in general! One important particular case where this is
### guaranteed though is when 'isNormal(subject)' is TRUE (i.e. 'subject' is
### a normal XStringViews object) and 'max.mismatch=0' (no "out of limits"
### matches).
setMethod("matchPattern", "XStringViews",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .XStringViews.matchPattern(pattern, subject, algorithm,
                                   max.mismatch, min.mismatch, with.indels, fixed)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPattern", "MaskedXString",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        matchPattern(pattern, toXStringViewsOrXString(subject), algorithm=algorithm,
                     max.mismatch=max.mismatch, min.mismatch=min.mismatch,
                     with.indels=with.indels, fixed=fixed)
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
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        standardGeneric("countPattern")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "character",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .XString.matchPattern(pattern, subject, algorithm,
                              max.mismatch, min.mismatch, with.indels, fixed,
                              count.only=TRUE)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "XString",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .XString.matchPattern(pattern, subject, algorithm,
                              max.mismatch, min.mismatch, with.indels, fixed,
                              count.only=TRUE)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "XStringSet",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        stop("please use vcountPattern() when 'subject' is an XStringSet object (multiple sequence)")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "XStringViews",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .XStringViews.matchPattern(pattern, subject, algorithm,
                                   max.mismatch, min.mismatch, with.indels, fixed,
                                   count.only=TRUE)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "MaskedXString",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        countPattern(pattern, toXStringViewsOrXString(subject), algorithm=algorithm,
                     max.mismatch=max.mismatch, min.mismatch=min.mismatch,
                     with.indels=with.indels, fixed=fixed)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "vmatchPattern" and "vcountPattern" generic and methods.
###
### These are vectorized versions of matchPattern() and countPattern().
### vmatchPattern() returns an MIndex object and vcountPattern() an integer
### vector (like matchPDict() and countPDict(), respectively).
###

.XStringSet.vmatchPattern <- function(pattern, subject, algorithm,
                                      max.mismatch, min.mismatch, with.indels, fixed,
                                      count.only=FALSE)
{
    if (!is(subject, "XStringSet"))
        subject <- XStringSet(NULL, subject)
    algo <- .normargAlgorithm(algorithm)
    if (.is.character.algo(algo)) 
        stop("'subject' must be a single (non-empty) string ", 
             "for this algorithm")
    pattern <- normargPattern(pattern, subject)
    max.mismatch <- normargMaxMismatch(max.mismatch)
    min.mismatch <- normargMinMismatch(min.mismatch, max.mismatch)
    with.indels <- normargWithIndels(with.indels)
    fixed <- normargFixed(fixed, subject)
    if (!isTRUEorFALSE(count.only)) 
        stop("'count.only' must be TRUE or FALSE")
    algo <- .select.algo(algo, pattern, max.mismatch, min.mismatch, with.indels, fixed)
    # because MIndex objects do not support variable-width matches yet
    if (algo == "indels" && !count.only)
        stop("vmatchPattern() does not support indels yet")
    C_ans <- .Call("XStringSet_vmatch_pattern", pattern, subject,
                   algo, max.mismatch, min.mismatch, with.indels, fixed, count.only,
                   PACKAGE = "Biostrings")
    if (count.only)
        return(C_ans)
    # C_ans is a list of IRanges objects
    ans_width0 <- rep.int(length(pattern), length(subject))
    new("ByPos_MIndex", width0=ans_width0, NAMES=names(subject), ends=lapply(C_ans, end))
}

setGeneric("vmatchPattern", signature="subject",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE, ...)
        standardGeneric("vmatchPattern")
)

setMethod("vmatchPattern", "character",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .XStringSet.vmatchPattern(pattern, subject, algorithm,
                                  max.mismatch, min.mismatch, with.indels, fixed)
)

setMethod("vmatchPattern", "XString",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        stop("please use matchPattern() when 'subject' is an XString object (single sequence)")
)

setMethod("vmatchPattern", "XStringSet",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .XStringSet.vmatchPattern(pattern, subject, algorithm,
                                  max.mismatch, min.mismatch, with.indels, fixed)
)

# TODO: Add a "vmatchPattern" method for XStringViews objects.
# Note that the start/end of the matches need to be returned as relative
# to subject(subject).
setMethod("vmatchPattern", "XStringViews",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        stop("XStringViews objects are not supported yet, sorry")
)

setMethod("vmatchPattern", "MaskedXString",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        stop("please use matchPattern() when 'subject' is a MaskedXString object (single sequence)")
)

setGeneric("vcountPattern", signature="subject",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE, ...)
        standardGeneric("vcountPattern")
)

setMethod("vcountPattern", "character",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .XStringSet.vmatchPattern(pattern, subject, algorithm,
                                  max.mismatch, min.mismatch, with.indels, fixed,
                                  count.only=TRUE)
)

setMethod("vcountPattern", "XString",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        stop("please use countPattern() when 'subject' is an XString object (single sequence)")
)

setMethod("vcountPattern", "XStringSet",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0L, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        .XStringSet.vmatchPattern(pattern, subject, algorithm,
                                  max.mismatch, min.mismatch, with.indels, fixed,
                                  count.only=TRUE)
)

setMethod("vcountPattern", "XStringViews",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        vcountPattern(pattern, XStringViewsToSet(subject, FALSE, verbose=FALSE),
                      algorithm=algorithm,
                      max.mismatch=max.mismatch, min.mismatch=min.mismatch,
                      with.indels=with.indels, fixed=fixed)
)

setMethod("vcountPattern", "MaskedXString",
    function(pattern, subject, algorithm="auto",
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
        stop("please use countPattern() when 'subject' is a MaskedXString object (single sequence)")
)

