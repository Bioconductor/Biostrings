### =========================================================================
### The "pmatchPattern" generic
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Longest Common Prefix: the "lcprefix" new generic
###

### 's1' and 's2' must be BString (or derived) objects of the same class.
### Return the length (integer) of the Longest Common Prefix.
BString.lcprefix <- function(s1, s2)
{
    .Call("lcprefix", s1@data@xp, s1@offset, s1@length,
                      s2@data@xp, s2@offset, s2@length,
                      PACKAGE="Biostrings")
}

setGeneric("lcprefix", signature=c("s1", "s2"),
    function(s1, s2) standardGeneric("lcprefix")
)
setMethod("lcprefix", signature(s1="character", s2="character"),
    function(s1, s2)
        BString.lcprefix(BString(s1), BString(s2))
)
setMethod("lcprefix", signature(s1="character", s2="BString"),
    function(s1, s2)
        BString.lcprefix(new(class(s2), s1), s2)
)
setMethod("lcprefix", signature(s1="BString", s2="character"),
    function(s1, s2)
        BString.lcprefix(s1, new(class(s1), s2))
)
setMethod("lcprefix", signature(s1="BString", s2="BString"),
    function(s1, s2)
    {
        if (class(s1) != class(s2))
            stop("'s1' and 's2' are not of the same class")
        BString.lcprefix(s1, s2)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Longest Common Suffix: the "lcsuffix" new generic
###

### 's1' and 's2' must be BString (or derived) objects of the same class.
### Return the length (integer) of the Longest Common Suffix.
BString.lcsuffix <- function(s1, s2)
{
    .Call("lcsuffix", s1@data@xp, s1@offset, s1@length,
                      s2@data@xp, s2@offset, s2@length,
                      PACKAGE="Biostrings")
}

setGeneric("lcsuffix", signature=c("s1", "s2"),
    function(s1, s2) standardGeneric("lcsuffix")
)
setMethod("lcsuffix", signature(s1="character", s2="character"),
    function(s1, s2)
        BString.lcsuffix(BString(s1), BString(s2))
)
setMethod("lcsuffix", signature(s1="character", s2="BString"),
    function(s1, s2)
        BString.lcsuffix(new(class(s2), s1), s2)
)
setMethod("lcsuffix", signature(s1="BString", s2="character"),
    function(s1, s2)
        BString.lcsuffix(s1, new(class(s1), s2))
)
setMethod("lcsuffix", signature(s1="BString", s2="BString"),
    function(s1, s2)
    {
        if (class(s1) != class(s2))
            stop("'s1' and 's2' are not of the same class")
        BString.lcsuffix(s1, s2)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "pmatchPattern" generic
###

### 'subject' must be a BString (or derived) object.
### Return a BStringPartialMatches object.
### Speed (on my Thinkpad):
###   > library(BSgenome.Scerevisiae.UCSC.sacCer1)
###   > Scerevisiae$chr1
###   > file <- system.file("extdata", "someORF.fsa", package="Biostrings")
###   > orf <- BStringViews(file(file), "DNAString")
###   > system.time(pmatchPattern(orf[[2]], Scerevisiae$chr1, max=1))
###      user  system elapsed
###     0.900   0.012   0.913
###   > system.time(pmatchPattern(orf[[2]], Scerevisiae$chr1, max=6))
###      user  system elapsed
###   102.910   2.460 105.689
.pmatchPattern.rec <- function(pattern, subject, maxlength.out)
{
    if (class(pattern) != class(subject))
        pattern <- new(class(subject), pattern)
    if (nchar(pattern) <= 20000) {
        sv <- matchPattern(pattern, subject)
        if (length(sv) >= maxlength.out || nchar(pattern) < 2L) {
            pviews <- data.frame(start=rep(1L, length(sv)), end=rep(nchar(pattern), length(sv)))
            pv <- new("BStringViews", subject=pattern, views=pviews)
            ans <- new("BStringPartialMatches", sv, subpatterns=pv)
            return(ans)
        }
    }
    ## Split 'pattern' in 2 parts
    Lnc <- nchar(pattern) %/% 2L
    Lpattern <- BString.substr(pattern, 1L, Lnc)
    Rpattern <- BString.substr(pattern, Lnc+1L, nchar(pattern))
    ## Recursive call on the left part
    Lpm <- .pmatchPattern.rec(Lpattern, subject, maxlength.out)
    Lpm_start <- start(Lpm)
    Lpm_end <- end(Lpm)
    Lpv <- subpatterns(Lpm)
    Lpv_end <- end(Lpv)
    Lindex <- which((Lpv_end == nchar(Lpattern)) & (Lpm_end < nchar(subject)))
    Loverlapping <- integer(0)
    for (i in Lindex) {
        overlap <- BString.lcprefix(Rpattern, BString.substr(subject, Lpm_end[i]+1L, nchar(subject)))
        if (overlap == 0L)
            next
        Loverlapping <- c(Loverlapping, i)
        Lpv_end[i] <- Lpv_end[i] + overlap
        Lpm_end[i] <- Lpm_end[i] + overlap
    }
    Lpv_start <- start(Lpv)
    ## Recursive call on the right part
    Rpm <- .pmatchPattern.rec(Rpattern, subject, maxlength.out)
    Rpm_start <- start(Rpm)
    Rpm_end <- end(Rpm)
    Rpv <- subpatterns(Rpm)
    Rpv_start <- start(Rpv)
    Rindex <- which((Rpv_start == 1L) & (Rpm_start > 1L))
    Roverlapping <- integer(0)
    for (i in Rindex) {
        overlap <- BString.lcsuffix(Lpattern, BString.substr(subject, 1L, Rpm_start[i]-1L))
        if (overlap == 0L)
            next
        Roverlapping <- c(Roverlapping, i)
        Rpv_start[i] <- Rpv_start[i] - overlap
        Rpm_start[i] <- Rpm_start[i] - overlap
    }
    Rpv_start <- Rpv_start + Lnc
    Rpv_end <- end(Rpv) + Lnc
    ## Remove duplicates (using a data frame might not be the most efficient
    ## way to achieve this)
    if (length(Loverlapping) != 0 && length(Roverlapping) != 0)
    {
        Loverlaps <- data.frame(start=Lpm_start[Loverlapping],
                                end=Lpm_end[Loverlapping],
                                spstart=Lpv_start[Loverlapping],
                                spend=Lpv_end[Loverlapping])
        Roverlaps <- data.frame(start=Rpm_start[Roverlapping],
                                end=Rpm_end[Roverlapping],
                                spstart=Rpv_start[Roverlapping],
                                spend=Rpv_end[Roverlapping])
        overlaps <- rbind(Loverlaps, Roverlaps)
        dup <- duplicated(overlaps)
        if (any(dup)) {
            which_dup <- Roverlapping[which(dup) - nrow(Loverlaps)]
            Rpm_start <- Rpm_start[-which_dup]
            Rpm_end <- Rpm_end[-which_dup]
            Rpv_start <- Rpv_start[-which_dup]
            Rpv_end <- Rpv_end[-which_dup]
        }
    }
    ## Merge left results with right results
    sviews <- data.frame(start=c(Lpm_start, Rpm_start), end=c(Lpm_end, Rpm_end))
    sv <- new("BStringViews", subject=subject, views=sviews)
    pviews <- data.frame(start=c(Lpv_start, Rpv_start), end=c(Lpv_end, Rpv_end))
    pv <- new("BStringViews", subject=pattern, views=pviews)
    ans <- new("BStringPartialMatches", sv, subpatterns=pv)
    ii <- order(width(ans), -start(ans), decreasing=TRUE)
    minwidth <- width(ans)[ii[1]] %/% 2
    ii <- ii[(width(ans)[ii] >= minwidth) | (seq_along(ii) <= maxlength.out)]
    ans[ii]
}

### 'subject' must be a BString (or derived) object.
### Return a BStringPartialMatches object.
.pmatchPattern <- function(pattern, subject, maxlength.out)
{
    ans <- .pmatchPattern.rec(pattern, subject, maxlength.out)
    if (length(ans) > maxlength.out)
        ans <- ans[seq_len(maxlength.out)]
    ans
}

setGeneric("pmatchPattern", signature="subject",
    function(pattern, subject, maxlength.out=1L)
        standardGeneric("pmatchPattern")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("pmatchPattern", "character",
    function(pattern, subject, maxlength.out)
    {
        .pmatchPattern(pattern, BString(subject), maxlength.out)
    }
)

### Dispatch on 'subject' (see signature of generic).
setMethod("pmatchPattern", "BString",
    function(pattern, subject, maxlength.out)
    {
        .pmatchPattern(pattern, subject, maxlength.out)
    }
)

### Dispatch on 'subject' (see signature of generic).
setMethod("pmatchPattern", "BStringViews",
    function(pattern, subject, maxlength.out)
    {
        if (length(subject) != 1)
            stop("'subject' must have a single view")
        .pmatchPattern(pattern, subject[[1]], maxlength.out)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Longest Common Substring: the "lcsubstr" new generic
###

### Implementation taken from 
###   http://en.wikibooks.org/wiki/Algorithm_implementation/Strings/Longest_common_substring
### 's1' and 's2' must be BString (or derived) objects of the same class.
BString.lcsubstr <- function(s1, s2)
{
    stop("coming soon...")
}

setGeneric("lcsubstr", signature=c("s1", "s2"),
    function(s1, s2) standardGeneric("lcsubstr")
)
setMethod("lcsubstr", signature(s1="character", s2="character"),
    function(s1, s2)
        BString.lcsubstr(BString(s1), BString(s2))
)
setMethod("lcsubstr", signature(s1="character", s2="BString"),
    function(s1, s2)
        BString.lcsubstr(new(class(s2), s1), s2)
)
setMethod("lcsubstr", signature(s1="BString", s2="character"),
    function(s1, s2)
        BString.lcsubstr(s1, new(class(s1), s2))
)
setMethod("lcsubstr", signature(s1="BString", s2="BString"),
    function(s1, s2)
    {
        if (class(s1) != class(s2))
            stop("'s1' and 's2' are not of the same class")
        BString.lcsubstr(s1, s2)
    }
)

