# Accepts an extra argument 'baseValues' that must be a named complex vector
# containing the values associated to each base e.g.
#   baseValues=c("A"=1+0i, "G"=0+1i, "T"=-1+0i, "C"=0-1i)
setMethod("as.complex", "DNAString",
    function(x, ...)
    {
        args <- list(...)
        if ("baseValues" %in% names(args))
            baseValues <- args$baseValues
        else
            baseValues <- c("A"=1+0i, "G"=0+1i, "T"=-1+0i, "C"=0-1i)
        z <- complex(nchar(x))
        for (i in 1:nchar(x)) {
            z[i] <- baseValues[as.character(x[i])]
        }
        z
    }
)

