### =========================================================================
### alphabetFrequency()
### strrev()
### mkAllStrings()
### oligonucleotideFrequency()
### dinucleotideFrequency()
### trinucleotideFrequency()
### -------------------------------------------------------------------------


.normalize.freq <- function(freq)
{
    if (!isTRUEorFALSE(freq))
        stop("'freq' must be 'TRUE' or 'FALSE'")
    freq
}

.normalize.collapse <- function(collapse)
{
    if (is.null(collapse))
        return(FALSE)
    if (!isTRUEorFALSE(collapse))
        stop("'collapse' must be 'TRUE' or 'FALSE'")
    collapse
}

.normalize.width <- function(width)
{
    if (!is.numeric(width) || length(width) != 1 || is.na(width))
        stop("'width' must be a single integer")
    width <- as.integer(width)
    if (width < 0)
        stop("'width' must be a non-negative integer")
    width
}

.normalize.fast.moving.side <- function(fast.moving.side)
{
    if (!is.character(fast.moving.side)
     || length(fast.moving.side) != 1
     || is.na(fast.moving.side))
        stop("'fast.moving.side' must be a single string")
    match.arg(fast.moving.side, c("left", "right"))
}

.normalize.as.array <- function(as.array, fast.moving.side)
{
    if (!isTRUEorFALSE(as.array))
        stop("'as.array' must be 'TRUE' or 'FALSE'")
    if (as.array && fast.moving.side != "left")
        stop("'fast.moving.side' must be \"left\" when 'as.array' is 'TRUE'")
    as.array
}

.normalize.with.labels <- function(with.labels)
{
    if (!isTRUEorFALSE(with.labels))
        stop("'with.labels' must be 'TRUE' or 'FALSE'")
    with.labels
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "alphabetFrequency" generic and methods.
###
### sum(alphabetFrequency(x)) should always be exactly nchar(x)
###

.XString.char_frequency <- function(x, freq)
{
    freq <- .normalize.freq(freq)
    ans <- .Call("XString_char_frequency", x, PACKAGE="Biostrings")
    if (freq)
        ans <- ans / nchar(x) # nchar(x) is sum(ans) but is faster
    ans
}

.XString.code_frequency <- function(x, baseOnly, freq)
{
    codes <- codes(x, baseOnly=baseOnly)
    freq <- .normalize.freq(freq)
    ans <- .Call("XString_code_frequency", x, codes, PACKAGE="Biostrings")
    names(ans) <- names(codes)
    if (baseOnly)
        ans <- c(ans, other=nchar(x)-sum(ans))
    if (freq)
        ans <- ans / nchar(x) # nchar(x) is sum(ans) but is faster
    ans
}

.XStringSet.char_frequency <- function(x, freq, collapse)
{
    freq <- .normalize.freq(freq)
    collapse <- .normalize.collapse(collapse)
    ans <- .Call("XStringSet_char_frequency", x, collapse, PACKAGE="Biostrings")
    if (freq)
        ans <- ans / nchar(x) # nchar(x) is sum(ans) but is faster
    ans
}

.XStringSet.code_frequency <- function(x, baseOnly, freq, collapse)
{
    codes <- codes(x, baseOnly=baseOnly)
    freq <- .normalize.freq(freq)
    collapse <- .normalize.collapse(collapse)
    ans <- .Call("XStringSet_code_frequency", x, collapse, codes, PACKAGE="Biostrings")
    names(ans) <- names(codes)
    if (baseOnly)
        ans <- rbind(ans, other=nchar(x)-sum(ans))
    if (freq)
        ans <- ans / nchar(x) # nchar(x) is sum(ans) but is faster
    ans
}

setGeneric("alphabetFrequency", signature="x",
    function(x, baseOnly=FALSE, freq=FALSE, ...) standardGeneric("alphabetFrequency")
)

setMethod("alphabetFrequency", "BString",
    function(x, baseOnly=FALSE, freq=FALSE, ...)
    {
        if (!missing(baseOnly))
            warning("'baseOnly' is ignored for a non DNA or RNA sequence")
        .XString.char_frequency(x, freq)
    }
)

setMethod("alphabetFrequency", "DNAString",
    function(x, baseOnly=FALSE, freq=FALSE, ...)
        .XString.code_frequency(x, baseOnly, freq)
)

setMethod("alphabetFrequency", "RNAString",
    function(x, baseOnly=FALSE, freq=FALSE, ...)
        .XString.code_frequency(x, baseOnly, freq)
)

setMethod("alphabetFrequency", "XStringSet",
    function(x, baseOnly=FALSE, freq=FALSE, ...)
    {
        if (is(super(x), "DNAString") || is(super(x), "RNAString")) {
            frequency <- function(v) alphabetFrequency(v, baseOnly=baseOnly, freq=FALSE)
        } else {
            if (!missing(baseOnly))
                warning("'baseOnly' is ignored for non DNA or RNA sequences")
            frequency <- function(v) alphabetFrequency(v, freq=FALSE)
        }
        freq <- .normalize.freq(freq)
        collapse <- .normalize.collapse(list(...)$collapse)
        ## Generate a zero-filled answer
        ans_row <- frequency(XString(class(super(x)), ""))
        if (collapse) {
            ans <- ans_row
            for (i in seq_len(length(x)))
                ans <- ans + frequency(x[[i]])
        } else {
            ans <- matrix(rep.int(ans_row, length(x)), ncol=length(ans_row), byrow=TRUE,
                                                       dimnames=list(NULL, names(ans_row)))
            for (i in seq_len(length(x)))
                ans[i, ] <- frequency(x[[i]])
            ## The "collapsed" result could also be obtained with colSums(ans)...
        }
        if (freq)
            ans <- ans / sum(ans)
        ans
    }
)

### library(drosophila2probe)
### dict0 <- drosophila2probe$sequence
### x <- BStringViews(as.character(dict0[1:2000]), subjectClass="DNAString")
### alphabetFrequency(x, baseOnly=TRUE)
### y <- DNAStringSet(x)
### alphabetFrequency(y, baseOnly=TRUE)
setMethod("alphabetFrequency", "BStringViews",
    function(x, baseOnly=FALSE, freq=FALSE, ...)
    {
        y <- BStringViewsToSet(x, use.names=FALSE, verbose=FALSE)
        alphabetFrequency(y, baseOnly=baseOnly, freq=freq, ...)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "strrev" function.
###

strrev <- function(x)
{
    if (length(x) == 0)
        return(x)
    sapply(strsplit(x, NULL, fixed=TRUE), function(xx) paste(rev(xx), collapse=""))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "mkAllStrings" function.
###

.mkAllStringsR <- function(alphabet, width)
{
    if (width == 0)
        return("")
    ansR <- .mkAllStringsR(alphabet, width - 1)
    unlist(lapply(alphabet, function(l) paste(l, ansR, sep="")))
}

.mkAllStringsL <- function(alphabet, width)
{
    if (width == 0)
        return("")
    ansL <- .mkAllStringsL(alphabet, width - 1)
    unlist(lapply(alphabet, function(l) paste(ansL, l, sep="")))
}

.mkAllStrings <- function(alphabet, width, fast.moving.side)
{
    if (fast.moving.side == "right")
        .mkAllStringsR(alphabet, width)
    else
        .mkAllStringsL(alphabet, width)
}

mkAllStrings <- function(alphabet, width, fast.moving.side="right")
{
    if (!is.character(alphabet))
        stop("'alphabet' must be a character vector")
    .mkAllStrings(alphabet,
                  .normalize.width(width),
                  .normalize.fast.moving.side(fast.moving.side))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "oligonucleotideFrequency" generic and methods.
###
### Except for the 'other' element, oligonucleotideFrequency(x, 1L)
### should be the same as alphabetFrequency(x, baseOnly=TRUE).
###
### Something else worth checking:
###   library(BSgenome.Dmelanogaster.FlyBase.r51)
###   chr3R <- Dmelanogaster[["3R"]]
###   width <- 12
###   dict0 <- mkAllStrings(names(Biostrings:::DNAcodes(TRUE)), width)
###   names(dict0) <- dict0
###   pdict <- PDict(dict0)
###   system.time(c1 <- countPDict(pdict, chr3R))
###   system.time(c2 <- oligonucleotideFrequency(chr3R, width, with.labels=FALSE))
###   identical(c1, c2) # must be TRUE
### Then try for other values of 'width' (1 <= width <= 12).
### Of course oligonucleotideFrequency() is much better: it is >10x faster, does
### not require preprocessing, and uses much less memory.
###

.formatFreqAnswer <- function(ans, alphabet, width, freq, fast.moving.side, as.array, with.labels)
{
    if (freq)
        ans <- ans / sum(ans)
    if (as.array) {
        if (with.labels)
            dimnames <- rep(list(alphabet), each=width)
        else
            dimnames <- NULL
        ans <- array(ans, dim=rep.int(4L, width), dimnames=dimnames)
    } else if (with.labels) {
        names(ans) <- .mkAllStrings(alphabet, width, fast.moving.side)
    }
    ans
}

.oligonucleotideFrequency <- function(x, width, freq, fast.moving.side, as.array, with.labels)
{
    width <- .normalize.width(width)
    freq <- .normalize.freq(freq)
    fast.moving.side <- .normalize.fast.moving.side(fast.moving.side)
    as.array <- .normalize.as.array(as.array, fast.moving.side)
    with.labels <- .normalize.with.labels(with.labels)
    base_codes <- codes(x, baseOnly=TRUE)
    ans <- .Call("oligonucleotide_frequency",
                 x@data@xp, x@offset, x@length,
                 base_codes, width, fast.moving.side,
                 PACKAGE="Biostrings")
    .formatFreqAnswer(ans, names(base_codes), width, freq, fast.moving.side, as.array, with.labels)
}

setGeneric("oligonucleotideFrequency", signature="x",
    function(x, width, freq=FALSE, fast.moving.side="right", as.array=FALSE, with.labels=TRUE, ...)
        standardGeneric("oligonucleotideFrequency")
)

setMethod("oligonucleotideFrequency", "DNAString",
    function(x, width, freq=FALSE, fast.moving.side="right", as.array=FALSE, with.labels=TRUE, ...)
    {
        if (missing(fast.moving.side) && !missing(as.array))
            fast.moving.side <- "left"
        .oligonucleotideFrequency(x, width, freq, fast.moving.side, as.array, with.labels)
    }
)

setMethod("oligonucleotideFrequency", "RNAString",
    function(x, width, freq=FALSE, fast.moving.side="right", as.array=FALSE, with.labels=TRUE, ...)
    {
        if (missing(fast.moving.side) && !missing(as.array))
            fast.moving.side <- "left"
        .oligonucleotideFrequency(x, width, freq, fast.moving.side, as.array, with.labels)
    }
)

setMethod("oligonucleotideFrequency", "XStringSet",
    function(x, width, freq=FALSE, fast.moving.side="right", as.array=FALSE, with.labels=TRUE, ...)
    {
        if (missing(fast.moving.side) && !missing(as.array))
            fast.moving.side <- "left"
        width <- .normalize.width(width)
        freq <- .normalize.freq(freq)
        fast.moving.side <- .normalize.fast.moving.side(fast.moving.side)
        as.array <- .normalize.as.array(as.array, fast.moving.side)
        with.labels <- .normalize.with.labels(with.labels)
        collapse <- .normalize.collapse(list(...)$collapse)
        base_codes <- codes(super(x), baseOnly=TRUE)
        ans <- integer(pow.int(4L, width))
        if (collapse) {
            for (i in seq_len(length(x))) {
                xx <- x[[i]]
                ans <- ans + .Call("oligonucleotide_frequency",
                                   xx@data@xp, xx@offset, xx@length,
                                   base_codes, width, fast.moving.side,
                                   PACKAGE="Biostrings")
            }
            ans <- .formatFreqAnswer(ans, names(base_codes), width, freq, fast.moving.side, as.array, with.labels)
        } else {
            ans <- rep.int(list(ans), length(x))
            for (i in seq_len(length(x))) {
                xx <- x[[i]]
                tmp <- .Call("oligonucleotide_frequency",
                             xx@data@xp, xx@offset, xx@length,
                             base_codes, width, fast.moving.side,
                             PACKAGE="Biostrings")
                ans[[i]] <- .formatFreqAnswer(tmp, names(base_codes), width, freq, fast.moving.side, as.array, with.labels)
            }
        }
        ans
    }
)

setMethod("oligonucleotideFrequency", "BStringViews",
    function(x, width, freq=FALSE, fast.moving.side="right", as.array=FALSE, with.labels=TRUE, ...)
    {
        y <- BStringViewsToSet(x, use.names=FALSE, verbose=FALSE)
        oligonucleotideFrequency(y, width, freq=freq, fast.moving.side=fast.moving.side,
                                           as.array=as.array, with.labels=with.labels, ...)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "dinucleotideFrequency" and "trinucleotideFrequency" convenience
### wrappers.
###

dinucleotideFrequency <- function(x, freq=FALSE, fast.moving.side="right", as.matrix=FALSE, with.labels=TRUE, ...)
{
    if (missing(fast.moving.side) && !missing(as.matrix))
        fast.moving.side <- "left"
    oligonucleotideFrequency(x, 2, freq=freq,
                                fast.moving.side=fast.moving.side,
                                as.array=as.matrix, with.labels=with.labels, ...)
}

trinucleotideFrequency <- function(x, freq=FALSE, fast.moving.side="right", as.array=FALSE, with.labels=TRUE, ...)
{
    if (missing(fast.moving.side) && !missing(as.array))
        fast.moving.side <- "left"
    oligonucleotideFrequency(x, 3, freq=freq,
                                fast.moving.side=fast.moving.side,
                                as.array=as.array, with.labels=with.labels, ...)
}

