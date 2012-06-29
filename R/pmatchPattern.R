### =========================================================================
### The "pmatchPattern" generic
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Longest Common Prefix: the "lcprefix" new generic
###

### 's1' and 's2' must be XString objects containing sequences of the same
### type. Return the length (integer) of the Longest Common Prefix.
XString.lcprefix <- function(s1, s2)
{
    .Call2("lcprefix", s1@shared@xp, s1@offset, s1@length,
                      s2@shared@xp, s2@offset, s2@length,
                      PACKAGE="Biostrings")
}

setGeneric("lcprefix", signature=c("s1", "s2"),
    function(s1, s2) standardGeneric("lcprefix")
)
setMethod("lcprefix", signature(s1="character", s2="character"),
    function(s1, s2)
        XString.lcprefix(BString(s1), BString(s2))
)
setMethod("lcprefix", signature(s1="character", s2="XString"),
    function(s1, s2)
        XString.lcprefix(XString(seqtype(s2), s1), s2)
)
setMethod("lcprefix", signature(s1="XString", s2="character"),
    function(s1, s2)
        XString.lcprefix(s1, XString(seqtype(s1), s2))
)
setMethod("lcprefix", signature(s1="XString", s2="XString"),
    function(s1, s2)
    {
        if (seqtype(s1) != seqtype(s2))
            stop("'s1' and 's2' must be XString objects containing ",
                 "sequences of the same type")
        XString.lcprefix(s1, s2)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Longest Common Suffix: the "lcsuffix" new generic
###

### 's1' and 's2' must be XString objects containing sequences of the same
### type. Return the length (integer) of the Longest Common Suffix.
XString.lcsuffix <- function(s1, s2)
{
    .Call2("lcsuffix", s1@shared@xp, s1@offset, s1@length,
                      s2@shared@xp, s2@offset, s2@length,
                      PACKAGE="Biostrings")
}

setGeneric("lcsuffix", signature=c("s1", "s2"),
    function(s1, s2) standardGeneric("lcsuffix")
)
setMethod("lcsuffix", signature(s1="character", s2="character"),
    function(s1, s2)
        XString.lcsuffix(BString(s1), BString(s2))
)
setMethod("lcsuffix", signature(s1="character", s2="XString"),
    function(s1, s2)
        XString.lcsuffix(XString(seqtype(s2), s1), s2)
)
setMethod("lcsuffix", signature(s1="XString", s2="character"),
    function(s1, s2)
        XString.lcsuffix(s1, XString(seqtype(s1), s2))
)
setMethod("lcsuffix", signature(s1="XString", s2="XString"),
    function(s1, s2)
    {
        if (seqtype(s1) != seqtype(s2))
            stop("'s1' and 's2' must be XString objects containing ",
                 "sequences of the same type")
        XString.lcsuffix(s1, s2)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "pmatchPattern" generic
###

### 'subject' must be an XString object.
### Return an XStringPartialMatches object.
### Speed (on my Thinkpad):
###   > library(BSgenome.Scerevisiae.UCSC.sacCer1)
###   > Scerevisiae$chr1
###   > file <- system.file("extdata", "someORF.fa", package="Biostrings")
###   > orf <- readDNAStringSet(file)
###   > system.time(pmatchPattern(orf[[2]], Scerevisiae$chr1, max=1))
###      user  system elapsed
###     0.900   0.012   0.913
###   > system.time(pmatchPattern(orf[[2]], Scerevisiae$chr1, max=6))
###      user  system elapsed
###   102.910   2.460 105.689
.pmatchPattern.rec <- function(pattern, subject, maxlength.out)
{
    pattern <- normargPattern(pattern, subject)
    if (nchar(pattern) <= 20000) {
        sv <- matchPattern(pattern, subject)
        if (length(sv) >= maxlength.out || nchar(pattern) < 2L) {
            pv <- unsafe.newXStringViews(pattern,
                                         rep.int(1L, length(sv)),
                                         rep.int(nchar(pattern), length(sv)))
            ans <- new("XStringPartialMatches", sv, subpatterns=pv)
            return(ans)
        }
    }
    ## Split 'pattern' in 2 parts
    Lnc <- nchar(pattern) %/% 2L
    Lpattern <- subseq(pattern, start=1L, end=Lnc)
    Rpattern <- subseq(pattern, start=Lnc+1L, end=nchar(pattern))
    ## Recursive call on the left part
    Lpm <- .pmatchPattern.rec(Lpattern, subject, maxlength.out)
    Lpm_start <- start(Lpm)
    Lpm_end <- end(Lpm)
    Lpv <- subpatterns(Lpm)
    Lpv_end <- end(Lpv)
    Lindex <- which((Lpv_end == nchar(Lpattern)) & (Lpm_end < nchar(subject)))
    Loverlapping <- integer(0)
    for (i in Lindex) {
        overlap <- XString.lcprefix(Rpattern, subseq(subject, start=Lpm_end[i]+1L, end=nchar(subject)))
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
        overlap <- XString.lcsuffix(Lpattern, subseq(subject, start=1L, end=Rpm_start[i]-1L))
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
    sv_start <- c(Lpm_start, Rpm_start)
    sv_width <- c(Lpm_end, Rpm_end) - sv_start + 1L
    sv <- unsafe.newXStringViews(subject, sv_start, sv_width)
    pv_start <- c(Lpv_start, Rpv_start)
    pv_width <- c(Lpv_end, Rpv_end) - pv_start + 1L
    pv <- unsafe.newXStringViews(pattern, pv_start, pv_width)
    ans <- new("XStringPartialMatches", sv, subpatterns=pv)
    ii <- order(width(ans), -start(ans), decreasing=TRUE)
    minwidth <- width(ans)[ii[1]] %/% 2
    ii <- ii[(width(ans)[ii] >= minwidth) | (seq_along(ii) <= maxlength.out)]
    ans[ii]
}

### 'subject' must be an XString object.
### Return an XStringPartialMatches object.
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
setMethod("pmatchPattern", "XString",
    function(pattern, subject, maxlength.out)
    {
        .pmatchPattern(pattern, subject, maxlength.out)
    }
)

### Dispatch on 'subject' (see signature of generic).
setMethod("pmatchPattern", "XStringViews",
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

### 's1' and 's2' must be XString objects containing sequences of the same
### type. Implementation taken from 
###   http://en.wikibooks.org/wiki/Algorithm_implementation/Strings/Longest_common_substring
XString.lcsubstr <- function(s1, s2)
{
    stop("coming soon...")
}

setGeneric("lcsubstr", signature=c("s1", "s2"),
    function(s1, s2) standardGeneric("lcsubstr")
)
setMethod("lcsubstr", signature(s1="character", s2="character"),
    function(s1, s2)
        XString.lcsubstr(BString(s1), BString(s2))
)
setMethod("lcsubstr", signature(s1="character", s2="XString"),
    function(s1, s2)
        XString.lcsubstr(XString(seqtype(s2), s1), s2)
)
setMethod("lcsubstr", signature(s1="XString", s2="character"),
    function(s1, s2)
        XString.lcsubstr(s1, XString(seqtype(s1), s2))
)
setMethod("lcsubstr", signature(s1="XString", s2="XString"),
    function(s1, s2)
    {
        if (seqtype(s1) != seqtype(s2))
            stop("'s1' and 's2' must be XString objects containing ",
                 "sequences of the same type")
        XString.lcsubstr(s1, s2)
    }
)

