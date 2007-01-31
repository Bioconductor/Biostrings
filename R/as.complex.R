### Accepts an extra argument 'baseValues' that must be a named complex vector
### containing the values associated to each base e.g.
###   baseValues=c("A"=1+0i, "G"=0+1i, "T"=-1+0i, "C"=0-1i)
setMethod("as.complex", "DNAString",
    function(x, ...)
    {
        args <- list(...)
        baseValues <- args[["baseValues"]]
        if (is.null(baseValues))
            baseValues <- c(A=1+0i, G=0+1i, T=-1+0i, C=0-1i)
        lookup_table <- complex(length(DNA_ALPHABET))
        for (letter in names(baseValues)) {
            code <- DNAString(letter)@data[] # dirty trick, need to find something better
            lookup_table[1 + code] <- baseValues[letter]
        }
        CharBuffer.readComplexes(x@data, x@offset + 1, x@offset + x@length, lookup_table)
    }
)

