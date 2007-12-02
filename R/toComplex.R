### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "toComplex" new generic.
###

setGeneric("toComplex",
    function(x, baseValues)
        standardGeneric("toComplex")
)

### 'baseValues' must be a named complex vector containing the values
### associated to each base e.g.
###   baseValues=c(A=1+0i, G=0+1i, T=-1+0i, C=0-1i)
setMethod("toComplex", "DNAString",
    function(x, baseValues)
    {
        letters <- names(baseValues)
        # dirty trick, need to find something better
        codes <- DNAString(paste(letters, collapse=""))@data[]
        lkup <- buildLookupTable(codes, baseValues)
        XRaw.readComplexes(x@data, x@offset + 1, x@offset + x@length, lkup)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "as.complex" generic is DEFUNCT in favor of "toComplex".
###

setMethod("as.complex", "DNAString", function(x, ...) .Defunct("toComplex"))

