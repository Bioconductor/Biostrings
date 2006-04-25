.alphabetFrequency <- function(x)
{
    .Call("alphabetFrequency",
          x@data@xp, x@offset, x@length,
          PACKAGE="Biostrings")
}

.DNAorRNA_alphabetFrequency <- function(x, baseOnly, codec)
{
    af <- .alphabetFrequency(x)
    if (baseOnly)
        i <- c(1:4,16)
    else
        i <- 1:length(codec@letters)
    ans <- af[codec@codes[i] + 1]
    names <- codec@letters[i]
    if (baseOnly) {
        names <- c(names, "other")
        ans <- c(ans, nchar(x) - sum(ans))
    }
    names(ans) <- names
    ans
}

setGeneric(
    "alphabetFrequency",
    function(x, baseOnly=TRUE) standardGeneric("alphabetFrequency")
)

setMethod("alphabetFrequency", "BString",
    function(x, baseOnly)
    {
        .alphabetFrequency(x)
    }
)

setMethod("alphabetFrequency", "DNAString",
    function(x, baseOnly)
    {
        .DNAorRNA_alphabetFrequency(x, baseOnly, DNA_STRING_CODEC)
    }
)

setMethod("alphabetFrequency", "RNAString",
    function(x, baseOnly)
    {
        .DNAorRNA_alphabetFrequency(x, baseOnly, RNA_STRING_CODEC)
    }
)

# Will fail if x contains "out of limits" views.
setMethod("alphabetFrequency", "BStringViews",
    function(x, baseOnly)
    {
        lx <- length(x)
        # Just a trick to generate a zero-filled answer
        ans <- alphabetFrequency(x@subject[1], baseOnly)
        ans[] <- 0
        if (lx == 0)
            return(ans)
        for (i in 1:lx)
            ans <- ans + alphabetFrequency(x[[i]], baseOnly)
        ans
    }
)

