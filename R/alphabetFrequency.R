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
### The "dinucleotideFrequency" generic and methods.
###

.all_multinucleotides <- function(letters, width)
{
    if (width == 0)
        return("")
    unlist(lapply(letters, function(l) paste(l, .all_multinucleotides(letters, width - 1), sep="")))
}

### .multinucleotide_frequency(x, 1L) should be exactly the same as
### alphabetFrequency(x, baseOnly=TRUE)
.multinucleotide_frequency <- function(x, width)
{
    base_codes <- codes(x, baseOnly=TRUE)
    ans <- .Call("multinucleotide_frequency",
                 x@data@xp, x@offset, x@length, base_codes, width,
                 PACKAGE="Biostrings")
    names(ans) <- .all_multinucleotides(names(base_codes), width)
    ans
}

setGeneric("dinucleotideFrequency", signature="x",
    function(x) standardGeneric("dinucleotideFrequency")
)

setMethod("dinucleotideFrequency", "DNAString",
    function(x)
    {
        .multinucleotide_frequency(x, 2L)
    }
)

setMethod("dinucleotideFrequency", "RNAString",
    function(x)
    {
        .multinucleotide_frequency(x, 2L)
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
        .multinucleotide_frequency(x, 3L)
    }
)

setMethod("trinucleotideFrequency", "RNAString",
    function(x)
    {
        .multinucleotide_frequency(x, 3L)
    }
)

