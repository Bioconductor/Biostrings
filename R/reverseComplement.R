### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "reverse" generic function and methods.
###

setGeneric("reverse", signature="x",
    function(x, ...) standardGeneric("reverse")
)

setMethod("reverse", "BString",
    function(x, ...)
    {
        lx <- length(x)
        data <- XRaw(lx)
        XRaw.reverseCopy(data, x@offset + 1, x@offset + lx, src=x@data)
        ## class(x) can be "BString", "DNAString", "RNAString" or "AAString"
        new(class(x), data, 0L, length(data), check=FALSE)
    }
)

setMethod("reverse", "BStringViews",
    function(x, ...)
    {
        subject <- reverse(x@subject)
        ls <- subject@length
        start <- ls - x@views$end + 1
        end <- ls - x@views$start + 1
        ans <- views(subject, start, end)
        desc(ans) <- desc(x)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "complement" generic function and methods.
###

setGeneric("complement", signature="x",
    function(x, ...) standardGeneric("complement")
)

setMethod("complement", "DNAString",
    function(x, ...)
    {
        lx <- length(x)
        data <- XRaw(lx)
        lkup <- getDNAComplementLookup()
        XRaw.copy(data, x@offset + 1, x@offset + lx, src=x@data, lkup=lkup)
        new("DNAString", data, 0L, length(data), check=FALSE)
    }
)

setMethod("complement", "RNAString",
    function(x, ...)
    {
        lx <- length(x)
        data <- XRaw(lx)
        lkup <- getRNAComplementLookup()
        XRaw.copy(data, x@offset + 1, x@offset + lx, src=x@data, lkup=lkup)
        new("RNAString", data, 0L, length(data), check=FALSE)
    }
)

setMethod("complement", "BStringViews",
    function(x, ...)
    {
        subject <- complement(x@subject)
        ans <- views(subject, x@views$start, x@views$end)
        desc(ans) <- desc(x)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some convenience functions for transforming a DNA (or RNA) sequence into
### an RNA (or DNA) sequence.
###

transcribe <- function(x)
{
    if( !(is(x, "DNAString")) ) stop("transcribe only works on DNA input")
    RNAString(x)
}

cDNA <- function(x)
{
    if( !(is(x, "RNAString")) ) stop("cDNA only works on RNA input")
    DNAString(x)
}

dna2rna <- function(x) RNAString(complement(x))
 
rna2dna <- function(x) DNAString(complement(x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "reverseComplement" generic function and methods.
###

setGeneric("reverseComplement", signature="x",
    function(x, ...) standardGeneric("reverseComplement")
)

setMethod("reverseComplement", "DNAString",
    function(x, ...)
    {
        lx <- length(x)
        data <- XRaw(lx)
        lkup <- getDNAComplementLookup()
        XRaw.reverseCopy(data, x@offset + 1, x@offset + lx, src=x@data, lkup=lkup)
        new("DNAString", data, 0L, length(data), check=FALSE)
    }
)

setMethod("reverseComplement", "RNAString",
    function(x, ...)
    {
        lx <- length(x)
        data <- XRaw(lx)
        lkup <- getRNAComplementLookup()
        XRaw.reverseCopy(data, x@offset + 1, x@offset + lx, src=x@data, lkup=lkup)
        new("RNAString", data, 0L, length(data), check=FALSE)
    }
)

setMethod("reverseComplement", "BStringViews",
    function(x, ...)
    {
        subject <- reverseComplement(x@subject)
        ls <- subject@length
        start <- ls - x@views$end + 1
        end <- ls - x@views$start + 1
        ans <- views(subject, start, end)
        desc(ans) <- desc(x)
        ans
    }
)

