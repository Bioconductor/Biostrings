### =========================================================================
### Preprocessed Subject Strings of type BOC
### ------------------------------------------------------------


### Note that 'base1_code', 'base2_code' and 'base3_code' must be distinct
setClass("BOC_SubjectString",
    representation(
        subject="DNAString",        # TODO: support "RNAString" too
        pattern_length="integer",   # A single integer e.g. 36
        base1_code="integer",       # integer belonging to DNA_BASE_CODES
        base2_code="integer",
        base3_code="integer",
        base4_code="integer",
        base1_OCbuffer="XRaw",      # all buffers must be of length nchar(subject) - pattern_length + 1
        base2_OCbuffer="XRaw",
        base3_OCbuffer="XRaw",
	pre4buffer="XRaw",
        ## The "stats" slot is a named list with the following elements:
        ##   means: vector of 4 doubles
        ##   table1, table2, table3, table4: vectors of (pattern_length + 1) integers
        stats="list"
    )
)

### Typical use:
###   library(BSgenome.Hsapiens.UCSC.hg18)
###   chr1 <- Hsapiens$chr1
###   chr1boc <- new("BOC_SubjectString", chr1, 36, c("A", "C", "G")) # 3-4 seconds on lamb1
setMethod("initialize", "BOC_SubjectString",
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
        buf1 <- XRaw(buf_length)
        buf2 <- XRaw(buf_length)
        buf3 <- XRaw(buf_length)
        pre4buf <- XRaw(buf_length)
        stats <- .Call("match_BOC_preprocess",
              subject@data@xp, subject@offset, subject@length,
              pattern_length,
              code1, code2, code3, code4,
              buf1@xp, buf2@xp, buf3@xp, pre4buf@xp,
              PACKAGE="Biostrings")
        .Object@base1_code <- code1
        .Object@base2_code <- code2
        .Object@base3_code <- code3
        .Object@base4_code <- code4
        .Object@base1_OCbuffer <- buf1
        .Object@base2_OCbuffer <- buf2
        .Object@base3_OCbuffer <- buf3
        .Object@pre4buffer <- pre4buf
        .Object@stats <- stats
        .Object
    }
)

### Typical use:
###   Biostrings:::plotBOC(chr1boc, "Human chr1")
plotBOC <- function(x, main)
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
### The "matchPattern" method for BOC_SubjectString objects.
###
### Typical use:
###   library(BSgenome.Hsapiens.UCSC.hg18)
###   chr1 <- Hsapiens$chr1
###   chr1boc <- new("BOC_SubjectString", chr1, 36, c("A", "C", "G"))
###   matchPattern(chr1[1:36], chr1boc)
###
### Performance (kind of disappointing so far):
###   for (i in 41:60) matchPattern(chr1[1:36+1000000*i], chr1boc)
###   #--> takes about 11 seconds on lamb1
###   for (i in 41:60) matchPattern(chr1[1:36+1000000*i], chr1, algo="boyer-moore")
###   #--> takes about 7.6 seconds on lamb1
###   for (i in 41:60) matchPattern(chr1[1:36+1000000*i], chr1, algo="naive-exact")
###   #--> takes about 111 seconds on lamb1
###   

debug_BOC <- function()
{
    invisible(.Call("match_BOC_debug", PACKAGE="Biostrings"))
}

### Must return an integer vector.
.match.BOC.exact <- function(pattern, boc_subject, count.only)
{
    .Call("match_BOC_exact",
          pattern@data@xp, pattern@offset, pattern@length,
          boc_subject@subject@data@xp, boc_subject@subject@offset, boc_subject@subject@length,
          boc_subject@base1_code,
          boc_subject@base2_code,
          boc_subject@base3_code,
          boc_subject@base4_code,
          boc_subject@base1_OCbuffer@xp,
          boc_subject@base2_OCbuffer@xp,
          boc_subject@base3_OCbuffer@xp,
          boc_subject@pre4buffer@xp,
          boc_subject@stats, count.only,
          PACKAGE="Biostrings")
}

### Must return an integer vector.
.match.BOC.inexact <- function(pattern, boc_subject, mismatch, count.only)
{
    stop("NOT READY YET!")
}

### Dispatch on 'subject' (see signature of generic).
### 'algorithm' is ignored.
setMethod("matchPattern", "BOC_SubjectString",
    function(pattern, subject, algorithm, mismatch, fixed)
    {
        if (class(pattern) != class(subject@subject))
            pattern <- new(class(subject@subject), pattern)
        pattern_length <- nchar(pattern)
        if (pattern_length != subject@pattern_length)
            stop("subject was preprocessed for patterns of length ", subject@pattern_length)
        if (!missing(fixed)) {
            fixed <- .normalize.fixed(fixed)
            if (!all(fixed))
                stop("only 'fixed=TRUE' can be used with a subject of class ", class(subject))
        }
        if (mismatch == 0)
            matches <- .match.BOC.exact(pattern, subject, count.only=FALSE)
        else
            matches <- .match.BOC.inexact(pattern, subject, mismatch, count.only=FALSE)
        new("BStringViews", subject=subject@subject,
            start=matches, end=matches+pattern_length-1L, check.views=FALSE)
    }
)

