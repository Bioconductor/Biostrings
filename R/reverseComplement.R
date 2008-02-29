### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "reverse" generic function and methods.
###

setGeneric("reverse", signature="x",
    function(x, ...) standardGeneric("reverse")
)

setMethod("reverse", "BString",
    function(x, ...)
        BString.tr(x, reverse=TRUE)
)

setMethod("reverse", "BStringViews",
    function(x, ...)
    {
        subject <- reverse(x@subject)
        ls <- subject@length
        start <- ls - end(x) + 1
        end <- ls - start(x) + 1
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
        BString.tr(x, lkup=getDNAComplementLookup())
)

setMethod("complement", "RNAString",
    function(x, ...)
        BString.tr(x, lkup=getRNAComplementLookup())
)

setMethod("complement", "BStringViews",
    function(x, ...)
    {
        subject <- complement(x@subject)
        ans <- views(subject, start(x), end(x))
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
    if (!is(x, "DNAString")) stop("transcribe() only works on DNA input")
    RNAString(complement(x))
}

cDNA <- function(x)
{
    if (!is(x, "RNAString")) stop("cDNA() only works on RNA input")
    DNAString(complement(x))
}

dna2rna <- function(x)
{
    if (!is(x, "DNAString")) stop("dna2rna() only works on DNA input")
    RNAString(x)
}
 
rna2dna <- function(x)
{
    if (!is(x, "RNAString")) stop("rna2dna() only works on RNA input")
    DNAString(x)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "reverseComplement" generic function and methods.
###

setGeneric("reverseComplement", signature="x",
    function(x, ...) standardGeneric("reverseComplement")
)

setMethod("reverseComplement", "DNAString",
    function(x, ...)
        BString.tr(x, lkup=getDNAComplementLookup(), reverse=TRUE)
)

setMethod("reverseComplement", "RNAString",
    function(x, ...)
        BString.tr(x, lkup=getRNAComplementLookup(), reverse=TRUE)
)

setMethod("reverseComplement", "BStringViews",
    function(x, ...)
    {
        subject <- reverseComplement(x@subject)
        ls <- subject@length
        start <- ls - end(x) + 1L
        end <- ls - start(x) + 1L
        ans <- views(subject, start, end)
        desc(ans) <- desc(x)
        ans
    }
)

