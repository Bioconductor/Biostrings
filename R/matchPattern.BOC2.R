### =========================================================================
### Preprocessed Subject Strings of type BOC2
### -------------------------------------------------------------------------


### Note that 'base1_code', 'base2_code' and 'base3_code' must be distinct
setClass("BOC2_SubjectString",
    representation(
        subject="DNAString",        # TODO: support "RNAString" too
        pattern_length="integer",   # A single integer e.g. 36
        base1_code="integer",       # integer belonging to DNA_BASE_CODES
        base2_code="integer",
        base3_code="integer",
        base4_code="integer",
        buffer="XInteger",          # must be of length nchar(subject) - pattern_length + 1
        ## The "stats" slot is a named list with the following elements:
        ##   means: vector of 4 doubles
        ##   table1, table2, table3, table4: vectors of (pattern_length + 1) integers
        stats="list"
    )
)

debug_BOC2 <- function()
{
    invisible(.Call("match_BOC2_debug", PACKAGE="Biostrings"))
}

### Typical use:
###   library(BSgenome.Hsapiens.UCSC.hg18)
###   chr1 <- Hsapiens$chr1
###   chr1boc <- new("BOC2_SubjectString", chr1, 36, c("A", "C", "G")) # 3-4 seconds on lamb1
setMethod("initialize", "BOC2_SubjectString",
    function(.Object, subject, pattern_length, base_letters)
    {
        .Object@subject <- subject
        if (!is.numeric(pattern_length) || length(pattern_length) != 1 || is.na(pattern_length))
            stop("'pattern_length' must be a single integer")
        pattern_length <- as.integer(pattern_length)
        if (pattern_length < 4L || 254L < pattern_length)
            stop("'pattern_length' must be >= 4 and <= 254")
        if (pattern_length > nchar(subject))
            stop("'pattern_length' must be <= 'nchar(subject)'")
        .Object@pattern_length <- pattern_length
        if (!is.character(base_letters) || length(base_letters) != 3
         || !all(base_letters %in% names(DNA_BASE_CODES)) || any(duplicated(base_letters)))
            stop("'base_letters' must contain 3 distinct DNA base-letters")
        buf_length <- nchar(subject) - pattern_length + 1
        code1 <- DNA_BASE_CODES[base_letters[1]]
        code2 <- DNA_BASE_CODES[base_letters[2]]
        code3 <- DNA_BASE_CODES[base_letters[3]]
        code4 <- DNA_BASE_CODES[setdiff(names(DNA_BASE_CODES), base_letters)]
        buf <- XInteger(buf_length)
        stats <- .Call("match_BOC2_preprocess",
              subject@data@xp, subject@offset, subject@length,
              pattern_length,
              code1, code2, code3, code4,
              buf@xp,
              PACKAGE="Biostrings")
        .Object@base1_code <- code1
        .Object@base2_code <- code2
        .Object@base3_code <- code3
        .Object@base4_code <- code4
        .Object@buffer <- buf
        .Object@stats <- stats
        .Object
    }
)

### Typical use:
###   Biostrings:::plotBOC2(chr1boc, "Human chr1")
plotBOC2 <- function(x, main)
{
    XLAB <- "Base Occurence Count"
    TITLE <- paste(XLAB, " for the ", x@pattern_length, "-mers in ", main, sep="")
    YLAB <- paste("number of ", x@pattern_length, "-mers", sep="")
    YMAX <- max(c(x@stats$table1, x@stats$table2, x@stats$table3, x@stats$table4))
    plot.new()
    plot.window(c(0, x@pattern_length), c(0, YMAX))
    title(TITLE, xlab=XLAB, ylab=YLAB, col.main="black")
    axis(1)
    axis(2)
    axis(4)

    par(fg="red")
    lines(0:x@pattern_length, x@stats$table1, type="l")
    par(fg="blue")
    lines(0:x@pattern_length, x@stats$table2, type="l")
    par(fg="green")
    lines(0:x@pattern_length, x@stats$table3, type="p")
    par(fg="black")
    lines(0:x@pattern_length, x@stats$table4, type="p")

    LEGEND <- c(names(x@base1_code), names(x@base2_code), names(x@base3_code), names(x@base4_code))
    LEGEND <- paste(LEGEND, "-count", sep="")
    COLORS <- c("red", "blue", "green", "black")
    legend(x=x@pattern_length, y=YMAX, legend=LEGEND, col=COLORS, lty="solid", lwd=3, xjust=1.0, yjust=1.0)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "matchPattern" method for BOC2_SubjectString objects.
###
### Typical use:
###   library(BSgenome.Hsapiens.UCSC.hg18)
###   chr1 <- Hsapiens$chr1
###   chr1boc <- new("BOC2_SubjectString", chr1, 36, c("A", "C", "G"))
###   matchPattern(chr1[1:36], chr1boc)
###
### Performance (a little bit better than with the first BOC algo):
###   for (i in 41:80) matchPattern(chr1[1:36+1000000*i], chr1boc)
###   #--> takes about 12.8 seconds on lamb1
###   for (i in 41:80) matchPattern(chr1[1:36+1000000*i], chr1, algo="boyer-moore")
###   #--> takes 14.77 seconds on lamb1
###
### See next section below for detailed benchmarks.
###   

### Must return an integer vector.
.match.BOC2.exact <- function(pattern, boc_subject, count.only)
{
    .Call("match_BOC2_exact",
          pattern@data@xp, pattern@offset, pattern@length,
          boc_subject@subject@data@xp, boc_subject@subject@offset, boc_subject@subject@length,
          boc_subject@base1_code,
          boc_subject@base2_code,
          boc_subject@base3_code,
          boc_subject@base4_code,
          boc_subject@buffer@xp,
          boc_subject@stats, count.only,
          PACKAGE="Biostrings")
}

### Must return an integer vector.
.match.BOC2.inexact <- function(pattern, boc_subject, max.mismatch, count.only)
{
    stop("NOT READY YET!")
}

.matchPattern.BOC2 <- function(pattern, boc_subject, max.mismatch, fixed, count.only=FALSE)
{
    if (class(pattern) != class(boc_subject@subject))
        pattern <- new(class(boc_subject@subject), pattern)
    pattern_length <- nchar(pattern)
    if (pattern_length != boc_subject@pattern_length)
        stop("subject was preprocessed for patterns of length ", boc_subject@pattern_length)
    max.mismatch <- normalize.max.mismatch(max.mismatch)
    fixed <- normalize.fixed(fixed, class(subject))
    if (!all(fixed))
        stop("only 'fixed=TRUE' can be used with a subject of class ", class(boc_subject))
    if (max.mismatch == 0)
        matches <- .match.BOC2.exact(pattern, boc_subject, count.only)
    else
        matches <- .match.BOC2.inexact(pattern, boc_subject, max.mismatch, count.only)
    if (count.only)
        return(matches)
    new("BStringViews", subject=boc_subject@subject,
        start=matches, end=matches+nchar(pattern)-1L, check.views=FALSE)
}

### Dispatch on 'subject' (see signature of generic).
### 'algorithm' is ignored.
setMethod("matchPattern", "BOC2_SubjectString",
    function(pattern, subject, algorithm, max.mismatch, fixed)
    {
        .matchPattern.BOC2(pattern, subject, max.mismatch, fixed)
    }
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "BOC2_SubjectString",
    function(pattern, subject, algorithm, max.mismatch, fixed)
    {
        .matchPattern.BOC2(pattern, subject, max.mismatch, fixed, count.only=TRUE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Benchmarks.
###
### Running the cmp_BOC2vsBoyerMoore_exactmatching() function with different
### values of 'pattern_length' gives the following results:
###
###   pattern_length  boyer-moore       BOC2
###   --------------  -----------  ---------
###               12       11.432      6.419
###               20        8.691      6.393
###               36        7.539      6.386
###               48        6.754      6.387
###               60        5.893      6.394
###               92        5.509      6.368
###              132        4.867      6.384
###              180        5.326      6.371
###              236        4.448      6.390
###              254        4.539      6.368
###
### Note: the BOC2 algo only supports pattern lengths <= 254 (no such limit
### for Boyer-Moore).
###

### Search for 20 patterns in Human chromosome 1
cmp_BOC2vsBoyerMoore_exactmatching <- function(pattern_length)
{
    library(BSgenome.Hsapiens.UCSC.hg18)
    chr1 <- Hsapiens$chr1
    chr1boc <- new("BOC2_SubjectString", chr1, pattern_length, c("A", "C", "G"))
    dt0 <- system.time(for (i in 41:60) matchPattern(chr1[seq_len(pattern_length)+1000000*i], chr1, algo="boyer-moore"))
    dt1 <- system.time(for (i in 41:60) matchPattern(chr1[seq_len(pattern_length)+1000000*i], chr1boc))
    c('boyer-moore'=dt0[['elapsed']], 'BOC2'=dt1[['elapsed']])
}

### A note about [40-120]mers in Human genome with a surprisingly high number
### of occurences.
###
### TODO: Move this to a more appropiate place (vignette?)
###
### In Human chr1, the 50999961-51000124 region is rich in substring of 40-80
### letters that have a lot of occurences (a few dozens) in the whole
### chromosome sequence. More remarkably, these substrings also have a similar
### number of occurences in the minus strand _and_ in the plus and minus
### strands of other chromosomes (checked with chr2 and chr3 only, need to do
### more checking).
###
### For example, this:
###
###   library(BSgenome.Hsapiens.UCSC.hg18)
###   chr1 <- Hsapiens$chr1
###   chr1boc <- new("BOC2_SubjectString", chr1, 120, c("A", "C", "G"))
###   chr2 <- Hsapiens$chr2
###   chr2boc <- new("BOC2_SubjectString", chr2, 120, c("A", "C", "G"))
###   chr3 <- Hsapiens$chr3
###   chr3boc <- new("BOC2_SubjectString", chr3, 120, c("A", "C", "G"))
###   Biostrings:::scan123(chr1[seq_len(120)+51000000-17], chr1boc, chr2boc, chr3boc)
###
### gives the following:
###
###   cp1p cp1m cp2p cp2m cp3p cp3m 
###     14   21   13   16   12   18
###
### The same with a pattern length of 40 instead of 120 gives:
###
###   cp1p cp1m cp2p cp2m cp3p cp3m 
###    212  212  239  242  243  231 
###

scan123 <- function(pattern, subject1, subject2, subject3, min.sum=-1)
{
    rc_pattern <- reverseComplement(pattern)
    cp1p <- countPattern(pattern, subject1)
    cp1m <- countPattern(rc_pattern, subject1)
    cp2p <- countPattern(pattern, subject2)
    cp2m <- countPattern(rc_pattern, subject2)
    cp3p <- countPattern(pattern, subject3)
    cp3m <- countPattern(rc_pattern, subject3)
    ans <- c('cp1p'=cp1p, 'cp1m'=cp1m, 'cp2p'=cp2p, 'cp2m'=cp2m, 'cp3p'=cp3p, 'cp3m'=cp3m)
    if (min.sum >= 0 && sum(ans) >= min.sum)
        cat(paste(ans, collapse=" "), "\n", sep="")
    ans
}

