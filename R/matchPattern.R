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
### The "naive" algos.
###

### Must return an integer vector.
.match.naive.exact <- function(pattern, subject, count.only)
{
    .Call("match_pattern",
          pattern, subject, "naive-exact",
          0L, c(TRUE, TRUE),
          count.only,
          PACKAGE="Biostrings")
}

### Must return an integer vector.
.match.naive.inexact <- function(pattern, subject, max.mismatch, fixed, count.only)
{
    .Call("match_pattern",
          pattern, subject, "naive-inexact",
          max.mismatch, fixed,
          count.only,
          PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Boyer-Moore
###

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
          pattern@xdata@xp, pattern@offset, pattern@length,
          subject@xdata@xp, subject@offset, subject@length,
          count.only,
          PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### shift-or

### Must return an integer vector.
.match.shiftor <- function(pattern, subject, max.mismatch, fixed, count.only)
{
    ## We treat the edge-cases at the R level
    p <- length(pattern)
    if (p <= max.mismatch) {
        if (count.only)
            return(length(subject) + p - as.integer(1))
        return((1-p):length(subject))
    }
    if (p > max.mismatch + length(subject)) {
        if (count.only)
            return(as.integer(0))
        return(integer(0))
    }
    .Call("match_shiftor",
          pattern@xdata@xp, pattern@offset, pattern@length,
          subject@xdata@xp, subject@offset, subject@length,
          max.mismatch, fixed, count.only,
          PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .matchPattern()
###

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
    if (!isSingleString(algorithm))
        stop("'algorithm' must be a single string")
    algo <- match.arg(algorithm, c("auto", "gregexpr", "gregexpr2",
                                   "naive-exact", "naive-inexact",
                                   "boyer-moore", "shift-or"))
    if (algo %in% c("gregexpr", "gregexpr2")) {
        if (!isSingleString(subject) || nchar(subject) == 0)
            stop("for algorithms \"gregexpr\" and \"gregexpr2\" ",
                 "'subject' must be a single (and non-empty) string")
        if (!isSingleString(pattern) || nchar(pattern) == 0)
            stop("'pattern' must be a single (and non-empty) string")
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
    if (algo %in% c("gregexpr", "gregexpr2")) {
        if (max.mismatch != 0 || !all(fixed))
            stop("algorithms \"gregexpr\" and \"gregexpr2\" only support ",
                 "'max.mismatch=0' and 'fixed=TRUE'")
    } else {
        algos <- .valid.algos(pattern, max.mismatch, fixed)
        if (algo == "auto") {
            algo <- algos[1]
        } else {
            if (!(algo %in% algos))
                stop("valid algos for your problem (best suited first): ",
                     paste(paste("\"", algos, "\"", sep=""), collapse=", "))
        }
    }
    matches <- switch(algo,
        "gregexpr"=.match.gregexpr(pattern, subject, count.only),
        "gregexpr2"=.match.gregexpr2(pattern, subject, count.only),
        "naive-exact"=.match.naive.exact(pattern, subject, count.only),
        "naive-inexact"=.match.naive.inexact(pattern, subject, max.mismatch, fixed,
                                         count.only),
        "boyer-moore"=.match.boyermoore(pattern, subject, count.only),
        "shift-or"=.match.shiftor(pattern, subject, max.mismatch, fixed,
                                  count.only)
    )
    if (count.only)
        return(matches)
    if (algo == "gregexpr" || algo == "gregexpr2")
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

