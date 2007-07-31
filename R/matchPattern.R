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
### The "naive" methods

debug_naive <- function()
{
    invisible(.Call("match_naive_debug", PACKAGE="Biostrings"))
}

### Must return an integer vector.
.match.naive.exact <- function(pattern, subject, count.only)
{
    .Call("match_naive_exact",
          pattern@data@xp, pattern@offset, pattern@length,
          subject@data@xp, subject@offset, subject@length,
          count.only,
          PACKAGE="Biostrings")
}

### Must return an integer vector.
.match.naive.fuzzy <- function(pattern, subject, mismatch, fixed, count.only)
{
    ## We treat the edge-cases at the R level
    p <- length(pattern)
    if (p <= mismatch) {
        if (count.only)
            return(length(subject) + p - as.integer(1))
        return((1-p):length(subject))
    }
    if (p > mismatch + length(subject)) {
        if (count.only)
            return(as.integer(0))
        return(integer(0))
    }
    .Call("match_naive_fuzzy",
          pattern@data@xp, pattern@offset, pattern@length,
          subject@data@xp, subject@offset, subject@length,
          mismatch, fixed, count.only,
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
        return((1-p):length(subject))
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
### .matchPattern()
###

### Return a logical vector of length 2.
.normalize.fixed <- function(fixed)
{
    if (!is.logical(fixed) && !is.character(fixed))
        stop("'fixed' not a logical or character vector")
    if (is.logical(fixed)) {
        if (any(is.na(fixed)))
            stop("'fixed' has NAs")
        fixed_names <- names(fixed)
        if (is.null(fixed_names)) {
            if (!(length(fixed) %in% 1:2))
                stop("when an unamed logical vector, ",
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
    fixed
}

### Return a character vector containing the valid algos (best suited first)
### for the given problem (problem is described by the values of 'pattern',
### 'mismatch' and 'fixed').
### Raise an error if the problem "doesn't make sense".
### Make sure that:
###   1. 'pattern' is of the same class as 'subject'
###   2. 'mismatch' is a non-negative integer
###   3. 'fixed' has been normalized
### before you call .valid.algos()
.valid.algos <- function(pattern, mismatch, fixed)
{
    if (!all(fixed) && !(class(pattern) %in% c("DNAString", "RNAString")))
        stop("'fixed' value only supported for a DNAString or RNAString subject ",
             "(you can only use 'fixed=TRUE' with your subject)")
    algos <- character(0)
    if (mismatch == 0 && all(fixed)) {
        algos <- c(algos, "boyer-moore")
        if (nchar(pattern) <= .Clongint.nbits())
            algos <- c(algos, "shift-or")
        algos <- c(algos, "naive-exact")
    } else {
        if (fixed[1] == fixed[2] && nchar(pattern) <= .Clongint.nbits())
            algos <- c(algos, "shift-or")
    }
    c(algos, "naive-fuzzy") # "naive-fuzzy" is universal but slow
}

.matchPattern <- function(pattern, subject, algorithm, mismatch, fixed,
                          count.only=FALSE)
{
    if (!is.character(algorithm) || length(algorithm) != 1 || is.na(algorithm))
        stop("'algorithm' must be a single string")
    algo <- match.arg(algorithm, c("auto", "gregexpr", "gregexpr2",
                                   "naive-exact", "naive-fuzzy",
                                   "boyer-moore", "shift-or"))
    if (algo %in% c("gregexpr", "gregexpr2")) {
        if (!is.character(subject))
            stop("algorithms \"gregexpr\" and \"gregexpr2\" are only ",
                 "supported for character strings")
        if (length(subject) != 1 || is.na(subject) || nchar(subject) == 0)
            stop("'subject' must be a single (non-NA, non-empty) string")
        if (!is.character(pattern) || length(pattern) != 1
         || is.na(pattern) || nchar(pattern) == 0)
            stop("'pattern' must be a single (non-NA, non-empty) string")
    } else {
        if (is.character(subject))
            subject <- BString(subject)
        if (class(pattern) != class(subject))
            pattern <- new(class(subject), pattern)
        if (nchar(pattern) > 20000)
            stop("patterns with more than 20000 letters are not supported, sorry")
    }
    if (!is.numeric(mismatch) || length(mismatch) != 1 || is.na(mismatch))
        stop("'mismatch' must be a single integer")
    mismatch <- as.integer(mismatch)
    if (mismatch < 0)
        stop("'mismatch' must be a non-negative integer")
    fixed <- .normalize.fixed(fixed)
    if (!is.logical(count.only) || length(count.only) != 1 || is.na(count.only))
        stop("'count.only' must be TRUE or FALSE")
    if (algo %in% c("gregexpr", "gregexpr2")) {
        if (mismatch != 0 || !all(fixed))
            stop("algorithms \"gregexpr\" and \"gregexpr2\" only support ",
                 "'mismatch=0' and 'fixed=TRUE'")
    } else {
        algos <- .valid.algos(pattern, mismatch, fixed)
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
        "naive-fuzzy"=.match.naive.fuzzy(pattern, subject, mismatch, fixed,
                                         count.only),
        "boyer-moore"=.match.boyermoore(pattern, subject, count.only),
        "shift-or"=.match.shiftor(pattern, subject, mismatch, fixed,
                                  count.only)
    )
    if (count.only)
        return(matches)
    if (algo == "gregexpr" || algo == "gregexpr2")
        return(matches)
    new("BStringViews", subject=subject,
                        views=data.frame(start=matches, end=matches+nchar(pattern)-1L))
}

### Typical use:
###   matchPattern("TG", DNAString("GTGACGTGCAT"))
###   matchPattern("TG", DNAString("GTGACGTGCAT"), algo="shift", mis=1)
### Edge cases:
###   matchPattern("---", DNAString("ACGTGCA"), mismatch=3)
###   matchPattern("---", DNAString("A"))
setGeneric(
    "matchPattern", signature="subject",
    function(pattern, subject, algorithm="auto", mismatch=0, fixed=TRUE)
        standardGeneric("matchPattern")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPattern", "character",
    function(pattern, subject, algorithm, mismatch, fixed)
    {
        .matchPattern(pattern, subject, algorithm, mismatch, fixed)
    }
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPattern", "BString",
    function(pattern, subject, algorithm, mismatch, fixed)
    {
        .matchPattern(pattern, subject, algorithm, mismatch, fixed)
    }
)

### Dispatch on 'subject' (see signature of generic).
### WARNING: Unlike the other "matchPattern" methods, the BStringViews object
### returned by this method is not guaranteed to have its views ordered from
### left to right in general! One important particular case where this is
### guaranteed though is when 'subject' is a normalized BStringViews object
### and 'mismatch=0' (no "out of limits" matches).
setMethod("matchPattern", "BStringViews",
    function(pattern, subject, algorithm, mismatch, fixed)
    {
        ans_start <- ans_end <- integer(0)
        for (i in seq_len(length(subject))) {
            pm <- .matchPattern(pattern, subject[[i]], algorithm, mismatch, fixed)
            offset <- start(subject)[i] - 1L
            ans_start <- c(ans_start, offset + start(pm))
            ans_end <- c(ans_end, offset + end(pm))
        }
        new("BStringViews", subject=subject(subject),
            views=data.frame(start=ans_start, end=ans_end))
    }
)

matchDNAPattern <- function(...) .Defunct("matchPattern")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### countPattern() is a slightly faster equivalent to length(matchPattern())

### Typical use:
###   countPattern("TG", DNAString("GTGACGTGCAT"))
###   countPattern("TG", DNAString("GTGACGTGCAT"), algo="shift", mis=1)
### Edge cases:
###   countPattern("---", DNAString("ACGTGCA"), mismatch=3)
###   countPattern("---", DNAString("A"))
setGeneric(
    "countPattern", signature="subject",
    function(pattern, subject, algorithm="auto", mismatch=0, fixed=TRUE)
        standardGeneric("countPattern")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "character",
    function(pattern, subject, algorithm, mismatch, fixed)
    {
        .matchPattern(pattern, subject, algorithm, mismatch, fixed, count.only=TRUE)
    }
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "BString",
    function(pattern, subject, algorithm, mismatch, fixed)
    {
        .matchPattern(pattern, subject, algorithm, mismatch, fixed, count.only=TRUE)
    }
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "BStringViews",
    function(pattern, subject, algorithm, mismatch, fixed)
    {
        sum(
            sapply(seq_len(length(subject)),
                   function(i) .matchPattern(pattern, subject[[i]],
                                             algorithm, mismatch, fixed,
                                             count.only=TRUE)
            )
        )
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
            l <- BString.substr(pattern, i, i)
            r <- BString.substr(subject, j, j)
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
           function(i) bsMismatch(pattern, x@subject, x@views$start[i], fixed))
}

setGeneric(
    "mismatch", signature="x",
    function(pattern, x, fixed=TRUE) standardGeneric("mismatch")
)

### Typical use:
###   mp <- matchPattern("TGA", DNAString("GTGACGTGCAT"), mismatch=2)
###   mismatch("TGA", mp)
### Dispatch on 'x' (see signature of generic).
setMethod("mismatch", "BStringViews",
    function(pattern, x, fixed)
    {
        if (class(pattern) != class(x@subject)) {
            pattern <- new(class(x@subject), pattern)
        }
        .mismatch(pattern, x, fixed)
    }
)

