### =========================================================================
### Mask objects
### -------------------------------------------------------------------------

setClass("Mask", contains="NormalIRanges")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###
### This is a trick to have the Mask validity method being the same as the
### NormalIRanges validity method. It's ugly and we should not need to do
### this because the Mask class inherits the validity method from its
### parent class (i.e. calling validObject() on a Mask object will call the
### validity method for NormalIRanges objects). But unfortunately, in the
### absence of a validity method for Mask objects, recursive (aka deep)
### validation of objects that have slots of type Mask is not working properly
### as reported here:
###   https://stat.ethz.ch/pipermail/r-devel/2008-April/049120.html
###

setValidity("Mask", getValidity(getClassDef(extends("Mask")[2])))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization and coercion.
###
### We need to explicitly define the "coerce,IRanges,Mask" method below,
### because the implicitly defined one doesn't check validity.
###

newEmptyMask <- function() as(newEmptyNormalIRanges(), "Mask")

setAs("IRanges", "Mask",
    function(from)
    {
        ## This calls our "coerce,IRanges,NormalIRanges" method which does
        ## check validity.
        from <- as(from, "NormalIRanges")
        ## This calls the implicitly defined "coerce,NormalIRanges,mask"
        ## which does not check anything but this is OK because a valid
        ## NormalIRanges object is always a valid Mask object.
        as(from, "Mask")
    }
)


### =========================================================================
### MaskedXString objects
### -------------------------------------------------------------------------
###
### Note that thinking of MaskedXString objects as particular XString objects
### might seem natural and therefore making the MaskedXString class an
### extension of the XString class is tempting... But in the end there would
### be no benefit in doing so because inheritance would not give us
### anything good out-of-the-box. This is because almost all XString methods
### are doing the wrong thing for MaskedXString objects and hence would need
### to be overwritten so it would not save us any work. Furthermore, some of
### the XString methods should be disabled for MaskedXString objects. But if
### MaskedXString objects are considered XString objects, then these methods
### need to be overwritten and call stop() so that they fail on MaskedXString
### objects. By not making the MaskedXString class an extension of the XString
### class, we don't have to do this.
###

### Not good (MaskedXString extends XString).
#setClass("MaskedXString",
#    contains="XString",
#    representation(
#        "VIRTUAL",
#        mask="Mask"
#    )
#)
### 4 direct "MaskedXString" derivations (no additional slot)
#setClass("MaskedBString", contains=c("MaskedXString", "BString"))
#setClass("MaskedDNAString", contains=c("MaskedXString", "DNAString"))
#setClass("MaskedRNAString", contains=c("MaskedXString", "RNAString"))
#setClass("MaskedAAString", contains=c("MaskedXString", "AAString"))

### Better (MaskedXString does NOT extend XString).
setClass("MaskedXString",
    representation(
        "VIRTUAL",
        unmasked="XString",
        mask="Mask"
    )
)

setClass("MaskedBString",
    contains="MaskedXString",
    representation(
        unmasked="BString"
    )
)
setClass("MaskedDNAString",
    contains="MaskedXString",
    representation(
        unmasked="DNAString"
    )
)
setClass("MaskedRNAString",
    contains="MaskedXString",
    representation(
        unmasked="RNAString"
    )
)
setClass("MaskedAAString",
    contains="MaskedXString",
    representation(
        unmasked="AAString"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setGeneric("unmasked", function(x) standardGeneric("unmasked"))
setMethod("unmasked", "MaskedXString", function(x) x@unmasked)

setGeneric("mask", function(x) standardGeneric("mask"))
setMethod("mask", "XString", function(x) NULL)
setMethod("mask", "MaskedXString", function(x) x@mask)

setMethod("alphabet", "MaskedXString",
    function(x) alphabet(unmasked(x))
)

setMethod("length", "MaskedXString",
    function(x) length(unmasked(x))
)

setMethod("nchar", "MaskedXString",
    function(x, type="chars", allowNA=FALSE)
    {
        length(x) - sum(width(mask(x)))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other non exported accessor methods.
###

setMethod("baseXStringSubtype", "MaskedXString",
    function(x) baseXStringSubtype(unmasked(x))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.MaskedXString.mask <- function(object)
{
    mask <- mask(object)
    if (length(mask) == 0)
        return(NULL)
    if (start(mask)[1] < 1 || length(object) < end(mask)[length(mask)])
        return("the mask contains \"out of limits\" ranges")
    NULL
}

setValidity("MaskedXString",
    function(object)
    {
        problems <- .valid.MaskedXString.mask(object)
        if (is.null(problems)) TRUE else problems
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

### From XString objects to MaskedXString objects.
setAs("BString", "MaskedBString",
    function(from) new("MaskedBString", unmasked=from, mask=newEmptyMask())
)
setAs("DNAString", "MaskedDNAString",
    function(from) new("MaskedDNAString", unmasked=from, mask=newEmptyMask())
)
setAs("RNAString", "MaskedRNAString",
    function(from) new("MaskedRNAString", unmasked=from, mask=newEmptyMask())
)
setAs("AAString", "MaskedAAString",
    function(from) new("MaskedAAString", unmasked=from, mask=newEmptyMask())
)

### From MaskedXString objects to XString objects.
setAs("MaskedBString", "BString",
    function(from) unmasked(from)
)
setAs("MaskedDNAString", "DNAString",
    function(from) unmasked(from)
)
setAs("MaskedRNAString", "RNAString",
    function(from) unmasked(from)
)
setAs("MaskedAAString", "AAString",
    function(from) unmasked(from)
)

### Others.
setMethod("as.character", "MaskedXString",
    function(x)
    {
        ans <- as.character(unmasked(x))
        mask <- mask(x)
        for (i in seq_len(length(mask))) {
            strip <- paste(rep.int("#", width(mask)[i]), collapse="")
            substr(ans,  start(mask)[i], end(mask)[i]) <- strip
        }
        ans
    }
)

setMethod("toString", "MaskedXString", function(x, ...) as.character(x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "mask<-" replacement methods.
###

setGeneric("mask<-", function(x, value) standardGeneric("mask<-"))

### Setting the mask of a MaskedXString object...
setReplaceMethod("mask", signature(x="MaskedXString", value="NULL"),
    function(x, value) unmasked(x)
)
setReplaceMethod("mask", signature(x="MaskedXString", value="Mask"),
    function(x, value)
    {
        x@mask <- restrict(value, 1L, length(x), use.names=FALSE)
        x
    }
)
setReplaceMethod("mask", signature(x="MaskedXString", value="NormalIRanges"),
    function(x, value)
    {
        mask(x) <- as(value, "Mask")
        x
    }
)
setReplaceMethod("mask", signature(x="MaskedXString", value="IRanges"),
    function(x, value)
    {
        mask(x) <- toNormalIRanges(value)
        x
    }
)
setReplaceMethod("mask", signature(x="MaskedXString", value="XString"),
    function(x, value)
    {
        ## We need to remove the current mask before calling matchPattern()
        mask(x) <- NULL
        mask(x) <- matchPattern(value, x)
        x
    }
)
setReplaceMethod("mask", signature(x="MaskedXString", value="character"),
    function(x, value)
    {
        mask(x) <- XString(baseXStringSubtype(x), value)
        x
    }
)

### Setting the mask of an XString object...
setReplaceMethod("mask", signature(x="XString", value="NULL"),
    function(x, value) x
)
setReplaceMethod("mask", signature(x="XString", value="ANY"),
    function(x, value)
    {
        x <- as(x, paste("Masked", baseXStringSubtype(x), sep=""))
        mask(x) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "MaskedXString.substr" function (NOT exported).
###
### Return a MaskedXString object (not vectorized).
### Like for XString.substr(), 'start' and 'end' must be single integers
### verifying:
###   1 <= start AND end <= length(x) AND end >= start - 1
### WARNING: This function is voluntarly unsafe (it doesn't check its
### arguments) because we want it to be the fastest possible!
###

MaskedXString.substr <- function(x, start, end)
{
    x@unmasked <- XString.substr(unmasked(x), start, end)
    mask <- mask(x)
    shift <- start - 1L
    .start(mask) <- start(mask) - shift
    mask(x) <- mask # will restrict the mask
    x
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "subseq" method for MaskedXString objects.
###

setMethod("subseq", "MaskedXString",
    function(x, start=NA, end=NA, width=NA)
    {
        limits <- new("IRanges", 1L, length(x))
        limits <- narrow(limits, start=start, end=end, width=width)
        MaskedXString.substr(x, start(limits), end(limits))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

setMethod("show", "MaskedXString",
    function(object)
    {
        lo <- length(object)
        nmasked <- sum(width(mask(object)))
        cat("  ", lo, "-letter (", nmasked, " masked with #) \"",
                  class(object), "\" instance\n", sep="")
        cat("seq:", toSeqSnippet(object, getOption("width") - 5))
        cat("\n")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "gaps" method.
###
### Allows mask inversion in a convenient way: mask(x) <- gaps(x)
###

### 'start' and 'end' are ignored.
setMethod("gaps", "MaskedXString",
    function(x, start=NA, end=NA)
    {
        gaps(mask(x), start=1, end=length(x))
    }
)

