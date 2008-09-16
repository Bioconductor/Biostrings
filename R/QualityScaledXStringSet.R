### =========================================================================
### QualityScaledXStringSet objects
### -------------------------------------------------------------------------
###

setClass("QualityScaledXStringSet",
    representation(
        "VIRTUAL",
        quality = "XStringQuality"
    )
)

### QualityScaledXStringSet subtypes
setClass("QualityScaledBStringSet",
    contains=c("BStringSet", "QualityScaledXStringSet")
)
setClass("QualityScaledDNAStringSet",
    contains=c("DNAStringSet", "QualityScaledXStringSet")
)
setClass("QualityScaledRNAStringSet",
    contains=c("RNAStringSet", "QualityScaledXStringSet")
)
setClass("QualityScaledAAStringSet",
    contains=c("AAStringSet", "QualityScaledXStringSet")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setGeneric("quality", function(x) standardGeneric("quality"), useAsDefault = function(x) x@quality)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The user-friendly versatile constructors.
###

QualityScaledXStringSet <- function(x, quality) {
    if (!is(x, "XStringSet"))
        stop("'x' must be of class 'XStringSet'")
    if (!is(x, "XStringSet"))
        stop("'quality' must be of class 'XStringQuality'")
    if (!(length(quality) %in% c(1, length(x))))
        stop("'length(quality)' must equal 1 or 'length(x)'")
    output <- as(x, paste("QualityScaled", class(x), sep=""))
    slot(output, "quality", check=FALSE) <- quality
    output
}

QualityScaledBStringSet <- function(x, quality) QualityScaledXStringSet(BStringSet(x), quality)
QualityScaledDNAStringSet <- function(x, quality) QualityScaledXStringSet(DNAStringSet(x), quality)
QualityScaledRNAStringSet <- function(x, quality) QualityScaledXStringSet(RNAStringSet(x), quality)
QualityScaledAAStringSet <- function(x, quality) QualityScaledXStringSet(AAStringSet(x), quality)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

setMethod("show", "QualityScaledXStringSet",
    function(object)
    {
        cat("  A ", class(object), " instance containing:\n", sep="")
        cat("\n")
        selectMethod("show", "XStringSet")(as(object, "XStringSet"))
        cat("\n")
        show(quality(object))
        cat("\n")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Supported 'i' types: numeric vector, logical vector, NULL and missing.
setMethod("[", "QualityScaledXStringSet",
    function(x, i, j, ..., drop)
    {
        slot(x, "quality", check=FALSE) <- quality(x)[i]
        x[i]
    }
)
