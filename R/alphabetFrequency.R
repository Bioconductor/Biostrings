setGeneric(
    "alphabetFrequency",
    function(x, baseOnly=TRUE) standardGeneric("alphabetFrequency")
)

# TODO: Implement a FAST alphabetFrequency() in C
# TODO: Make a separate function that is called by the 2 methods below
setMethod("alphabetFrequency", "DNAString",
    function(x, baseOnly)
    {
        if (baseOnly)
            letters <- DNA_STRING_CODEC@letters[c(1:4,16)]
        else
            letters <- DNA_STRING_CODEC@letters[]
        ans <- sapply(letters, function(let) countPattern(let, x, fixed=TRUE))
        if (baseOnly) {
            letters <- c(letters, "other")
            ans <- c(ans, nchar(x) - sum(ans))
        }
        names(ans) <- letters
        ans
    }
)

# TODO: Implement a FAST alphabetFrequency() in C
setMethod("alphabetFrequency", "RNAString",
    function(x, baseOnly)
    {
        if (baseOnly)
            letters <- RNA_STRING_CODEC@letters[c(1:4,16)]
        else
            letters <- RNA_STRING_CODEC@letters[]
        ans <- sapply(letters, function(let) countPattern(let, x, fixed=TRUE))
        if (baseOnly) {
            letters <- c(letters, "other")
            ans <- c(ans, nchar(x) - sum(ans))
        }
        names(ans) <- letters
        ans
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

