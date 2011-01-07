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
        if (is.null(names(baseValues)))
            stop("'baseValues' must have names")
        if (any(duplicated(names(baseValues))))
            stop("'baseValues' must have unique names")
        base_codes <- xscodes(x)
        if (!all(names(baseValues) %in% names(base_codes)))
            stop("'baseValues' names must be valid DNA letters")
        if (!is.complex(baseValues))
            class(baseValues) <- "complex" # as.complex() would drop the names!
        lkup <- buildLookupTable(base_codes[names(baseValues)], baseValues)
        SharedRaw.readComplexes(x@shared, x@offset + 1L, x@offset + x@length, lkup)
    }
)

