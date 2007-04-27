BString.char_frequency <- function(x)
{
    .Call("CharBuffer_char_frequency",
          x@data@xp, x@offset, x@length,
          PACKAGE="Biostrings")
}

BString.alphabet_frequency <- function(x, baseOnly)
{
    ans <- BString.char_frequency(x)
    codes <- codec(x)@codes
    names <- alphabet(x)
    if (baseOnly) {
        i <- c(1:4, 16)
        codes <- codes[i]
        names <- names[i]
    }
    ans <- ans[1 + codes]
    if (baseOnly) {
        names <- c(names, "other")
        ans <- c(ans, nchar(x) - sum(ans))
    }
    names(ans) <- names
    ans
}

setGeneric(
    "alphabetFrequency",
    function(x, baseOnly=FALSE) standardGeneric("alphabetFrequency")
)

setMethod("alphabetFrequency", "BString",
    function(x, baseOnly=FALSE)
    {
        BString.char_frequency(x)
    }
)

setMethod("alphabetFrequency", "DNAString",
    function(x, baseOnly=FALSE)
    {
        BString.alphabet_frequency(x, baseOnly)
    }
)

setMethod("alphabetFrequency", "RNAString",
    function(x, baseOnly=FALSE)
    {
        BString.alphabet_frequency(x, baseOnly)
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

