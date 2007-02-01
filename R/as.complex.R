### Accepts an extra argument 'baseValues' that must be a named complex vector
### containing the values associated to each base e.g.
###   baseValues=c(A=1+0i, G=0+1i, T=-1+0i, C=0-1i)
setMethod("as.complex", "DNAString",
    function(x, baseValues)
    {
        letters <- names(baseValues)
        # dirty trick, need to find something better
        codes <- DNAString(paste(letters, collapse=""))@data[]
        lkup <- buildLookupTable(codes, baseValues)
        CharBuffer.readComplexes(x@data, x@offset + 1, x@offset + x@length, lkup)
    }
)

