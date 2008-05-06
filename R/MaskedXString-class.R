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
#        masks="MaskCollection"
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
        masks="MaskCollection"
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

setGeneric("masks", function(x) standardGeneric("masks"))
setMethod("masks", "XString", function(x) NULL)
setMethod("masks", "MaskedXString", function(x) x@masks)

setMethod("alphabet", "MaskedXString",
    function(x) alphabet(unmasked(x))
)

setMethod("length", "MaskedXString",
    function(x) length(unmasked(x))
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

.valid.MaskedXString.unmasked <- function(object)
{
    if (!is(unmasked(object), "XString"))
        return("the 'unmasked' slot must contain an XString object")
    if (!is(object, paste("Masked", baseXStringSubtype(object), sep="")))
        return("bad XString subtype for the unmasked sequence")
    if (length(object) != width(masks(object)))
        return("the length of the object and the width of the mask collection differ")
    NULL
}

.valid.MaskedXString.masks <- function(object)
{
    masks <- masks(object)
    if (!is(masks, "MaskCollection"))
        return("the 'masks' slot must contain a MaskCollection object")
    if (width(masks) != length(object))
        return("the length of the object and the width of the mask collection differ")
    NULL
}

.valid.MaskedXString <- function(object)
{
    c(.valid.MaskedXString.unmasked(object),
      .valid.MaskedXString.masks(object))
}

setValidity("MaskedXString",
    function(object)
    {
        problems <- .valid.MaskedXString(object)
        if (is.null(problems)) TRUE else problems
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

### From XString objects to MaskedXString objects.
setAs("BString", "MaskedBString",
    function(from)
    {
        masks <- new("MaskCollection", width=length(from))
        new("MaskedBString", unmasked=from, masks=masks)
    }
)
setAs("DNAString", "MaskedDNAString",
    function(from)
    {
        masks <- new("MaskCollection", width=length(from))
        new("MaskedDNAString", unmasked=from, masks=masks)
    }
)
setAs("RNAString", "MaskedRNAString",
    function(from)
    {
        masks <- new("MaskCollection", width=length(from))
        new("MaskedRNAString", unmasked=from, masks=masks)
    }
)
setAs("AAString", "MaskedAAString",
    function(from)
    {
        masks <- new("MaskCollection", width=length(from))
        new("MaskedAAString", unmasked=from, masks=masks)
    }
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

### From a MaskedXString object to a MaskCollection object.
setAs("MaskedXString", "MaskCollection",
    function(from) masks(from)
)

### From a MaskedXString object to a NormalIRanges object.
setAs("MaskedXString", "NormalIRanges",
    function(from) as(masks(from), "NormalIRanges")
)

### From a MaskedXString object to an XStringViews object.
setAs("MaskedXString", "XStringViews",
    function(from)
    {
        views <- gaps(reduce(masks(from)))[[1]]
        ans_start <- start(views)
        ans_width <- width(views)
        new("XStringViews", unmasked(from), start=ans_start, width=ans_width, check=FALSE)
    }
)

### NOT exported.
toXStringViewsOrXString <- function(x)
{
    x0 <- unmasked(x)
    mask1 <- reduce(masks(x))
    if (isEmpty(mask1))
        return(x0)
    views <- gaps(mask1)[[1]]
    new("XStringViews", x0, start=start(views), width=width(views), check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "nchar" method.
###

setMethod("nchar", "MaskedXString",
    function(x, type="chars", allowNA=FALSE)
    {
        length(x) - maskedwidth(reduce(masks(x)))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The transformation methods (endomorphisms) "reduce" and "gaps".
###

### 'with.inframe.attrib' is ignored.
setMethod("reduce", "MaskedXString",
    function(x, with.inframe.attrib=FALSE)
    {
        x@masks <- reduce(masks(x))
        x
    }
)

### 'start' and 'end' are ignored.
setMethod("gaps", "MaskedXString",
    function(x, start=NA, end=NA)
    {
        x@masks <- gaps(reduce(masks(x)))
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "MaskedXString.substr" function (NOT exported).
###
### Return a MaskedXString object (not vectorized).
### Like for XString.substr(), 'start' and 'end' must be single integers
### verifying:
###   1 <= start AND end <= length(x) AND start <= end + 1
### WARNING: This function is voluntarly unsafe (it doesn't check its
### arguments) because we want it to be the fastest possible!
###

MaskedXString.substr <- function(x, start, end)
{
    ## The "narrow" method for MaskCollection objects does actually check
    ## that 'start' and 'end' are safe.
    x@masks <- narrow(masks(x), start=start, end=end)
    x@unmasked <- XString.substr(unmasked(x), start, end)
    x
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "subseq" method for MaskedXString objects.
###

### This method is used by the toSeqSnippet() function when called
### on a MaskedXString object.
setMethod("subseq", "MaskedXString",
    function(x, start=NA, end=NA, width=NA)
    {
        limits <- new("IRanges", start=1L, width=length(x), check=FALSE)
        limits <- narrow(limits, start=start, end=end, width=width)
        MaskedXString.substr(x, start(limits), end(limits))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "as.character" and "toString" methods.
###

### This method is used by the toSeqSnippet() function when called
### on a MaskedXString object.
setMethod("as.character", "MaskedXString",
    function(x)
    {
        ans <- as.character(unmasked(x))
        nir0 <- as(x, "NormalIRanges")
        for (i in seq_len(length(nir0))) {
            strip <- paste(rep.int("#", width(nir0)[i]), collapse="")
            substr(ans,  start(nir0)[i], end(nir0)[i]) <- strip
        }
        ans
    }
)

setMethod("toString", "MaskedXString", function(x, ...) as.character(x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

setMethod("show", "MaskedXString",
    function(object)
    {
        lo <- length(object)
        cat("  ", lo, "-letter \"", class(object), "\" instance (# for masking)\n", sep="")
        cat("seq:", toSeqSnippet(object, getOption("width") - 5))
        cat("\n")
        MaskCollection.show_frame(masks(object))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "masks<-" replacement methods.
###

setGeneric("masks<-", function(x, value) standardGeneric("masks<-"))

### Setting the masks of a MaskedXString object...
setReplaceMethod("masks", signature(x="MaskedXString", value="NULL"),
    function(x, value) unmasked(x)
)
setReplaceMethod("masks", signature(x="MaskedXString", value="MaskCollection"),
    function(x, value)
    {
        if (width(value) != length(x))
            stop("the width of the mask collection must be equal ",
                 "to the length of the sequence")
        x@masks <- value
        x
    }
)
setReplaceMethod("masks", signature(x="MaskedXString", value="XString"),
    function(x, value)
    {
        ## We need to remove the current masks before calling matchPattern()
        masks(x) <- NULL
        nir1 <- toNormalIRanges(matchPattern(value, x))
        name1 <- paste(as.character(value), "-blocks", sep="")
        masks(x) <- new("MaskCollection",
                        nir_list=list(nir1),
                        width=length(x),
                        active=TRUE,
                        NAMES=name1)
        x
    }
)
setReplaceMethod("masks", signature(x="MaskedXString", value="character"),
    function(x, value)
    {
        masks(x) <- XString(baseXStringSubtype(x), value)
        x
    }
)

### Setting the masks of an XString object...
setReplaceMethod("masks", signature(x="XString", value="NULL"),
    function(x, value) x
)
setReplaceMethod("masks", signature(x="XString", value="ANY"),
    function(x, value)
    {
        x <- as(x, paste("Masked", baseXStringSubtype(x), sep=""))
        masks(x) <- value
        x
    }
)

