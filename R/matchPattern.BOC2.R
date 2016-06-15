### =========================================================================
### Preprocessed Subject Strings of type BOC2  -  DEFUNCT!  DEFUNCT! DEFUNCT!
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

### Typical use:
###   library(BSgenome.Hsapiens.UCSC.hg18)
###   chr1 <- Hsapiens$chr1
###   chr1boc <- new("BOC2_SubjectString", chr1, 36, c("A", "C", "G")) # 3-4 seconds on lamb1
setMethod("initialize", "BOC2_SubjectString",
    function(.Object, subject, pattern_length, base_letters)
        .Defunct(msg="BOC2_SubjectString objects are defunct")
)

### Typical use:
###   Biostrings:::plotBOC2(chr1boc, "Human chr1")
plotBOC2 <- function(x, main)
{
    XLAB <- "Base Occurrence Count"
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

.match.BOC2.exact <- function(pattern, boc_subject, count.only)
{
    .Call2("match_BOC2_exact",
          pattern@shared@xp, pattern@offset, pattern@length,
          boc_subject@subject@shared@xp, boc_subject@subject@offset, boc_subject@subject@length,
          boc_subject@base1_code,
          boc_subject@base2_code,
          boc_subject@base3_code,
          boc_subject@base4_code,
          boc_subject@buffer@shared@xp,
          boc_subject@stats, count.only,
          PACKAGE="Biostrings")
}

.match.BOC2.inexact <- function(pattern, boc_subject, max.mismatch, count.only)
{
    stop("NOT READY YET!")
}

.matchPattern.BOC2 <- function(pattern, boc_subject, max.mismatch, fixed, count.only=FALSE)
{
    pattern <- normargPattern(pattern, boc_subject@subject)
    pattern_length <- nchar(pattern)
    if (pattern_length != boc_subject@pattern_length)
        stop("subject was preprocessed for patterns of length ", boc_subject@pattern_length)
    max.mismatch <- normargMaxMismatch(max.mismatch)
    fixed <- normargFixed(fixed, boc_subject@subject)
    if (!all(fixed))
        stop("only 'fixed=TRUE' can be used with a subject of class ", class(boc_subject))
    if (max.mismatch == 0)
        C_ans <- .match.BOC2.exact(pattern, boc_subject, count.only)
    else
        C_ans <- .match.BOC2.inexact(pattern, boc_subject, max.mismatch, count.only)
    if (count.only)
        return(C_ans)
    unsafe.newXStringViews(boc_subject@subject, start(C_ans), width(C_ans))
}

### Dispatch on 'subject' (see signature of generic).
### 'algorithm' is ignored.
setMethod("matchPattern", "BOC2_SubjectString",
    function(pattern, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto")
    {
        .Defunct(msg="BOC2_SubjectString objects are defunct")
    }
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "BOC2_SubjectString",
    function(pattern, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto")
    {
        .Defunct(msg="BOC2_SubjectString objects are defunct")
    }
)

