### =========================================================================
### alphabetFrequency(), dinucleotideFrequency() and trinucleotideFrequency()
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "alphabetFrequency" generic and methods.
###

### sum(.BString.char_frequency(x)) should always be exactly nchar(x)
.BString.char_frequency <- function(x)
{
    .Call("char_frequency",
          x@data@xp, x@offset, x@length,
          PACKAGE="Biostrings")
}

.BString.alphabet_frequency <- function(x, baseOnly)
{
    char_freq <- .BString.char_frequency(x)
    codes <- codes(x, baseOnly=baseOnly)
    ans <- char_freq[1 + codes]
    names(ans) <- names(codes)
    if (baseOnly)
        ans <- c(ans, other=nchar(x)-sum(ans))
    ans
}

setGeneric("alphabetFrequency", signature="x",
    function(x, baseOnly=FALSE) standardGeneric("alphabetFrequency")
)

setMethod("alphabetFrequency", "BString",
    function(x, baseOnly=FALSE)
    {
        .BString.char_frequency(x)
    }
)

setMethod("alphabetFrequency", "DNAString",
    function(x, baseOnly=FALSE)
    {
        .BString.alphabet_frequency(x, baseOnly)
    }
)

setMethod("alphabetFrequency", "RNAString",
    function(x, baseOnly=FALSE)
    {
        .BString.alphabet_frequency(x, baseOnly)
    }
)

### Will fail if x contains "out of limits" views.
setMethod("alphabetFrequency", "BStringViews",
    function(x, baseOnly=FALSE)
    {
        lx <- length(x)
        ## Just a trick to generate a zero-filled answer
        ans <- alphabetFrequency(x@subject[1], baseOnly)
        ans[] <- 0
        if (lx == 0)
            return(ans)
        for (i in 1:lx)
            ans <- ans + alphabetFrequency(x[[i]], baseOnly)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "mkAllStrings" function.
###
### Something worth checking:
###   library(BSgenome.Dmelanogaster.FlyBase.r51)
###   chr3R <- Dmelanogaster[["3R"]]
###   width <- 5
###   dict0 <- mkAllStrings(Biostrings:::codes(chr3R, baseOnly=TRUE), width)
###   names(dict0) <- dict0
###   pdict <- new("ULdna_PDict", dict0)
###   system.time(c1 <- countPDict(pdict, chr3R))
###   system.time(c2 <- Biostrings:::.all_oligonucleotide_frequency(chr3R, width))
###   identical(c1, c2) # must be TRUE
### Then try for other values of 'width' (1 <= width <= 10).
###

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
    if (!is.numeric(width) || length(width) != 1 || is.na(width))
        stop("'width' must be a single integer")
    width <- as.integer(width)
    if (width < 0)
        stop("'width' must be a non-negative integer")
    .mkAllStrings(alphabet, width)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "dinucleotideFrequency" generic and methods.
###

### .all_oligonucleotide_frequency(x, 1L) should be exactly the same as
### alphabetFrequency(x, baseOnly=TRUE)
.all_oligonucleotide_frequency <- function(x, width)
{
    base_codes <- codes(x, baseOnly=TRUE)
    ans <- .Call("all_oligonucleotide_frequency",
                 x@data@xp, x@offset, x@length, base_codes, width,
                 PACKAGE="Biostrings")
    names(ans) <- .mkAllStrings(names(base_codes), width)
    ans
}

setGeneric("dinucleotideFrequency", signature="x",
    function(x) standardGeneric("dinucleotideFrequency")
)

setMethod("dinucleotideFrequency", "DNAString",
    function(x)
    {
        .all_oligonucleotide_frequency(x, 2L)
    }
)

setMethod("dinucleotideFrequency", "RNAString",
    function(x)
    {
        .all_oligonucleotide_frequency(x, 2L)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "trinucleotideFrequency" generic and methods.
###

setGeneric("trinucleotideFrequency", signature="x",
    function(x) standardGeneric("trinucleotideFrequency")
)

setMethod("trinucleotideFrequency", "DNAString",
    function(x)
    {
        .all_oligonucleotide_frequency(x, 3L)
    }
)

setMethod("trinucleotideFrequency", "RNAString",
    function(x)
    {
        .all_oligonucleotide_frequency(x, 3L)
    }
)

