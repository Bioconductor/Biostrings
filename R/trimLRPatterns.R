### =========================================================================
### The "trimLRPatterns" generic.
### -------------------------------------------------------------------------


setGeneric("trimLRPatterns", signature = "subject",
    function(Lpattern = "", Rpattern = "", subject,
             max.Lmismatch = 0, max.Rmismatch = 0,
             with.Lindels = FALSE, with.Rindels = FALSE,
             Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE)
        standardGeneric("trimLRPatterns")
)

.normarg_maxLmismatch <- function(max.Lmismatch, Lpattern_len, LorR="L")
{
    argname <- paste0("max.", LorR, "mismatch", collapse="")
    if (!is.numeric(max.Lmismatch))
        stop("'", argname, "' must be numeric")
    if (any(is.na(max.Lmismatch)))
        stop("'", argname, "' cannot contain NAs")
    if (length(max.Lmismatch) == 1L
     && max.Lmismatch >= 0 && max.Lmismatch < 1) {
        max.Lmismatch <- as.integer(max.Lmismatch * seq_len(Lpattern_len))
    } else if (length(max.Lmismatch) > Lpattern_len) {
        stop("'", argname, "' cannot be longer than '", LorR, "pattern'")
    } else {
        max.Lmismatch <- as.integer(max.Lmismatch)
        if (length(max.Lmismatch) < Lpattern_len) {
            max.Lmismatch <-
                c(rep.int(-1L, Lpattern_len - length(max.Lmismatch)),
                  max.Lmismatch)
        }
    }
    rev(max.Lmismatch)
}

### 'subject' must be an XStringSet object of length != 0.
### Returns an integer vector of the same length as 'subject' where the i-th
### value is guaranteed to be >= 1 and <= width(subject)[i] + 1.
.computeTrimStart <- function(Lpattern, subject,
                              max.Lmismatch, with.Lindels, Lfixed)
{
    Lpattern <- normargPattern(Lpattern, subject, argname="Lpattern")
    Lpattern_len <- length(Lpattern)
    if (Lpattern_len == 0L)
        return(rep.int(1L, length(subject)))
    max.Lmismatch <- .normarg_maxLmismatch(max.Lmismatch, Lpattern_len,
                                           LorR="L")
    ## Test the pattern "from the inside out" (moving it to the left).
    ii <- which.isMatchingStartingAt(Lpattern,
                                     subject,
                                     starting.at = 1L,
                                     max.mismatch = max.Lmismatch,
                                     with.indels = with.Lindels,
                                     fixed = Lfixed,
                                     auto.reduce.pattern = TRUE)
    ii[is.na(ii)] <- Lpattern_len + 1L
    start <- Lpattern_len + 2L - ii
    if (length(start) == 0L)
        return(start)
    ## For elements in 'subject' shorter than 'Lpattern', 'start' can be
    ## > width(subject) + 1L.
    pmin(start, width(subject) + 1L)
}

### 'subject' must be an XStringSet object of length != 0.
### Returns an integer vector of the same length as 'subject' where the i-th
### value is guaranteed to be >= 0 and <= width(subject)[i].
.computeTrimEnd <- function(Rpattern, subject,
                            max.Rmismatch, with.Rindels, Rfixed)
{
    Rpattern <- normargPattern(Rpattern, subject, argname="Rpattern")
    Rpattern_len <- length(Rpattern)
    if (Rpattern_len == 0L)
        return(width(subject))
    ## Because we want to use which.isMatchingEndingAt() with
    ## 'auto.reduce.pattern=TRUE', the 'ending.at' arg will need to be a
    ## single value. But that won't be possible if 'subject' is not
    ## rectangular hence the ugly trick.
    if (!isConstant(width(subject))) {
        tmp <- .computeTrimStart(reverse(Rpattern), reverse(subject),
                                 max.Rmismatch, with.Rindels, Rfixed)
        return(width(subject) - tmp + 1L)
    }
    max.Rmismatch <- .normarg_maxLmismatch(max.Rmismatch, Rpattern_len,
                                           LorR="R")
    ## Test the pattern "from the inside out" (moving it to the right).
    subject_width <- width(subject)[1L]
    ii <- which.isMatchingEndingAt(pattern=Rpattern,
                                   subject=subject,
                                   ending.at=subject_width,
                                   max.mismatch=max.Rmismatch,
                                   with.indels=with.Rindels,
                                   fixed=Rfixed,
                                   auto.reduce.pattern=TRUE)
    ii[is.na(ii)] <- Rpattern_len + 1L
    end <- subject_width - Rpattern_len - 1L + ii
    if (length(end) == 0L)
        return(end)
    ## For elements in 'subject' shorter than 'Rpattern', 'end' can be < 0L.
    pmax(end, 0L)
}

.XStringSet.trimLRPatterns <- function(Lpattern, Rpattern, subject,
                                       max.Lmismatch, max.Rmismatch,
                                       with.Lindels, with.Rindels,
                                       Lfixed, Rfixed, ranges)
{
    if (!isTRUEorFALSE(ranges))
        stop("'ranges' must be TRUE or FALSE")
    if (length(subject) == 0L) {
        if (ranges)
            return(IRanges())
        return(subject)
    }
    start <- .computeTrimStart(Lpattern, subject,
                               max.Lmismatch, with.Lindels, Lfixed)
    end <- .computeTrimEnd(Rpattern, subject,
                           max.Rmismatch, with.Rindels, Rfixed)
    ## For those invalid ranges where 'start > end + 1L', we arbitrarily
    ## decide to set the 'start' to 'end + 1' (another reasonable choice
    ## would have been to set the 'end' to 'start - 1').
    idx <- which(start > end + 1L)
    start[idx] <- end[idx] + 1L
    if (ranges)
        return(IRanges(start=start, end=end))
    return(narrow(subject, start=start, end=end))
}

### Dispatch on 'subject' (see signature of generic).
setMethod("trimLRPatterns", "XString",
    function(Lpattern = "", Rpattern = "", subject,
             max.Lmismatch = 0, max.Rmismatch = 0,
             with.Lindels = FALSE, with.Rindels = FALSE,
             Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE)
    {
        subject <- as(subject, "XStringSet")
        ans <- .XStringSet.trimLRPatterns(Lpattern, Rpattern, subject,
                                          max.Lmismatch, max.Rmismatch,
                                          with.Lindels, with.Rindels,
                                          Lfixed, Rfixed, ranges)
        if (is(ans, "XStringSet"))
            ans <- ans[[1L]]
        ans
    }
)

### Dispatch on 'subject' (see signature of generic).
setMethod("trimLRPatterns", "XStringSet",
    function(Lpattern = "", Rpattern = "", subject,
             max.Lmismatch = 0, max.Rmismatch = 0,
             with.Lindels = FALSE, with.Rindels = FALSE,
             Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE)
    {
        .XStringSet.trimLRPatterns(Lpattern, Rpattern, subject,
                                   max.Lmismatch, max.Rmismatch,
                                   with.Lindels, with.Rindels,
                                   Lfixed, Rfixed, ranges)
    }
)

### Dispatch on 'subject' (see signature of generic).
setMethod("trimLRPatterns", "character",
    function(Lpattern = "", Rpattern = "", subject,
             max.Lmismatch = 0, max.Rmismatch = 0,
             with.Lindels = FALSE, with.Rindels = FALSE,
             Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE)
    {
        subject <- as(subject, "XStringSet")
        ans <- .XStringSet.trimLRPatterns(Lpattern, Rpattern, subject,
                                          max.Lmismatch, max.Rmismatch,
                                          with.Lindels, with.Rindels,
                                          Lfixed, Rfixed, ranges)
        if (is(ans, "XStringSet"))
            ans <- as.character(ans)
        ans
    }
)

