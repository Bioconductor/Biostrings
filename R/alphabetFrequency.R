### =========================================================================
### alphabetFrequency()
### strrev()
### mkAllStrings()
### allOligonucleotideFrequency()
### dinucleotideFrequency()
### trinucleotideFrequency()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "alphabetFrequency" generic and methods.
###
### sum(alphabetFrequency(x)) should always be exactly nchar(x)
###

.normalize.freq <- function(freq)
{
    if (!isTRUEorFALSE(freq))
        stop("'freq' must be 'TRUE' or 'FALSE'")
    freq
}

.BString.char_frequency <- function(x, freq)
{
    freq <- .normalize.freq(freq)
    ans <- .Call("char_frequency",
                 x@data@xp, x@offset, x@length,
                 PACKAGE="Biostrings")
    if (freq)
        ans <- ans / nchar(x) # nchar(x) is sum(ans) but is faster
    ans
}

.nucleotide_frequency <- function(x, baseOnly, freq)
{
    char_freq <- .BString.char_frequency(x, freq)
    codes <- codes(x, baseOnly=baseOnly)
    ans <- char_freq[1 + codes]
    names(ans) <- names(codes)
    if (baseOnly) {
        ## 'total_freq' is 'sum(char_freq)' but we can avoid to calculate
        ## this sum.
        if (freq)
            total_freq <- 1.0
        else
            total_freq <- nchar(x)
        ans <- c(ans, other=total_freq-sum(ans))
    }
    ans
}

setGeneric("alphabetFrequency", signature="x",
    function(x, baseOnly=FALSE, freq=FALSE) standardGeneric("alphabetFrequency")
)

setMethod("alphabetFrequency", "BString",
    function(x, baseOnly=FALSE, freq=FALSE)
    {
        if (!missing(baseOnly))
            warning("'baseOnly' is ignored for a non DNA or RNA sequence")
        .BString.char_frequency(x, freq)
    }
)

setMethod("alphabetFrequency", "DNAString",
    function(x, baseOnly=FALSE, freq=FALSE)
    {
        .nucleotide_frequency(x, baseOnly, freq)
    }
)

setMethod("alphabetFrequency", "RNAString",
    function(x, baseOnly=FALSE, freq=FALSE)
    {
        .nucleotide_frequency(x, baseOnly, freq)
    }
)

### Will fail if x contains "out of limits" views.
setMethod("alphabetFrequency", "BStringViews",
    function(x, baseOnly=FALSE, freq=FALSE)
    {
        if (is(x@subject, "DNAString") || is(x@subject, "RNAString")) {
            viewFrequency <- function(v) alphabetFrequency(v, baseOnly=baseOnly, freq=FALSE)
        } else {
            if (!missing(baseOnly))
                warning("'baseOnly' is ignored for views on a non DNA or RNA sequence")
            viewFrequency <- function(v) alphabetFrequency(v, freq=FALSE)
        }
        freq <- .normalize.freq(freq)
        ## Just a trick to generate a zero-filled answer
        ans <- viewFrequency(x@subject[1])
        ans[] <- 0L
        for (i in seq_len(length(x)))
            ans <- ans + viewFrequency(x[[i]])
        if (freq)
            ans <- ans / sum(ans)
        ans
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
### The "allOligonucleotideFrequency" generic and methods.
###
### Except for the 'other' element, allOligonucleotideFrequency(x, 1L)
### should be the same as alphabetFrequency(x, baseOnly=TRUE).
###
### Something else worth checking:
###   library(BSgenome.Dmelanogaster.FlyBase.r51)
###   chr3R <- Dmelanogaster[["3R"]]
###   width <- 5
###   dict0 <- mkAllStrings(names(Biostrings:::DNAcodes(TRUE)), width)
###   names(dict0) <- dict0
###   pdict <- new("ULdna_PDict", dict0)
###   system.time(c1 <- countPDict(pdict, chr3R))
###   system.time(c2 <- allOligonucleotideFrequency(chr3R, width))
###   identical(c1, c2) # must be TRUE
### Then try for other values of 'width' (1 <= width <= 10).
###

.all_oligonucleotide_frequency <- function(x, width, freq, fast.moving.side)
{
    width <- .normalize.width(width)
    freq <- .normalize.freq(freq)
    fast.moving.side <- .normalize.fast.moving.side(fast.moving.side)
    base_codes <- codes(x, baseOnly=TRUE)
    ans <- .Call("all_oligonucleotide_frequency",
                 x@data@xp, x@offset, x@length,
                 base_codes, width, fast.moving.side,
                 PACKAGE="Biostrings")
    if (freq)
        ans <- ans / sum(ans)
    names(ans) <- .mkAllStrings(names(base_codes), width, fast.moving.side)
    ans
}

setGeneric("allOligonucleotideFrequency", signature="x",
    function(x, width, freq=FALSE, fast.moving.side="right")
        standardGeneric("allOligonucleotideFrequency")
)

setMethod("allOligonucleotideFrequency", "DNAString",
    function(x, width, freq=FALSE, fast.moving.side="right")
    {
        .all_oligonucleotide_frequency(x, width, freq, fast.moving.side)
    }
)

setMethod("allOligonucleotideFrequency", "RNAString",
    function(x, width, freq=FALSE, fast.moving.side="right")
    {
        .all_oligonucleotide_frequency(x, width, freq, fast.moving.side)
    }
)

### Will fail if x contains "out of limits" views.
setMethod("allOligonucleotideFrequency", "BStringViews",
    function(x, width, freq=FALSE, fast.moving.side="right")
    {
        width <- .normalize.width(width)
        freq <- .normalize.freq(freq)
        fast.moving.side <- .normalize.fast.moving.side(fast.moving.side)
        base_codes <- codes(subject(x), baseOnly=TRUE)
        ans_names <- .mkAllStrings(names(base_codes), width, fast.moving.side)
        ans <- integer(length(ans_names))
        for (i in seq_len(length(x))) {
            xx <- x[[i]]
            ans <- ans + .Call("all_oligonucleotide_frequency",
                               xx@data@xp, xx@offset, xx@length,
                               base_codes, width, fast.moving.side,
                               PACKAGE="Biostrings")
        }
        if (freq)
            ans <- ans / sum(ans)
        names(ans) <- ans_names
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "dinucleotideFrequency" and "trinucleotideFrequency" convenience
### wrappers.
###

dinucleotideFrequency <- function(x, freq=FALSE, fast.moving.side="right")
{
    allOligonucleotideFrequency(x, 2, freq=freq, fast.moving.side=fast.moving.side)
}

trinucleotideFrequency <- function(x, freq=FALSE, fast.moving.side="right")
{
    allOligonucleotideFrequency(x, 3, freq=freq, fast.moving.side=fast.moving.side)
}

