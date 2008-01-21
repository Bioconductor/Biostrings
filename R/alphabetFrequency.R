### =========================================================================
### alphabetFrequency()
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
        if (!is(x@subject, "DNAString") && !is(x@subject, "RNAString") && !missing(baseOnly))
            warning("'baseOnly' is ignored for views on a non DNA or RNA sequence")
        freq <- .normalize.freq(freq)
        ## Just a trick to generate a zero-filled answer
        ans <- alphabetFrequency(x@subject[1], baseOnly, freq=FALSE)
        ans[] <- 0
        for (i in seq_len(length(x)))
            ans <- ans + alphabetFrequency(x[[i]], baseOnly, freq=FALSE)
        if (freq)
            ans <- ans / sum(ans)
        ans
    }
)


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

.mkAllStrings <- function(alphabet, width)
{
    if (width == 0)
        return("")
    tail <- .mkAllStrings(alphabet, width - 1)
    unlist(lapply(alphabet, function(l) paste(l, tail, sep="")))
}

mkAllStrings <- function(alphabet, width)
{
    if (!is.character(alphabet))
        stop("'alphabet' must be a character vector")
    width <- .normalize.width(width)
    .mkAllStrings(alphabet, width)
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
###   dict0 <- mkAllStrings(Biostrings:::codes(chr3R, baseOnly=TRUE), width)
###   names(dict0) <- dict0
###   pdict <- new("ULdna_PDict", dict0)
###   system.time(c1 <- countPDict(pdict, chr3R))
###   system.time(c2 <- allOligonucleotideFrequency(chr3R, width))
###   identical(c1, c2) # must be TRUE
### Then try for other values of 'width' (1 <= width <= 10).
###

.all_oligonucleotide_frequency <- function(x, width, freq)
{
    width <- .normalize.width(width)
    freq <- .normalize.freq(freq)
    base_codes <- codes(x, baseOnly=TRUE)
    ans <- .Call("all_oligonucleotide_frequency",
                 x@data@xp, x@offset, x@length, base_codes, width,
                 PACKAGE="Biostrings")
    if (freq)
        ans <- ans / sum(ans)
    names(ans) <- .mkAllStrings(names(base_codes), width)
    ans
}

setGeneric("allOligonucleotideFrequency", signature="x",
    function(x, width, freq=FALSE) standardGeneric("allOligonucleotideFrequency")
)

setMethod("allOligonucleotideFrequency", "DNAString",
    function(x, width, freq=FALSE)
    {
        .all_oligonucleotide_frequency(x, width, freq)
    }
)

setMethod("allOligonucleotideFrequency", "RNAString",
    function(x, width, freq=FALSE)
    {
        .all_oligonucleotide_frequency(x, width, freq)
    }
)

### Will fail if x contains "out of limits" views.
setMethod("allOligonucleotideFrequency", "BStringViews",
    function(x, width, freq=FALSE)
    {
        width <- .normalize.width(width)
        freq <- .normalize.freq(freq)
        base_codes <- codes(subject(x), baseOnly=TRUE)
        ans_names <- .mkAllStrings(names(base_codes), width)
        ans <- integer(length(ans_names))
        for (i in seq_len(length(x))) {
            xx <- x[[i]]
            ans <- ans + .Call("all_oligonucleotide_frequency",
                               xx@data@xp, xx@offset, xx@length, base_codes, width,
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

dinucleotideFrequency <- function(x, freq=FALSE) allOligonucleotideFrequency(x, 2, freq=freq)
trinucleotideFrequency <- function(x, freq=FALSE) allOligonucleotideFrequency(x, 3, freq=freq)

