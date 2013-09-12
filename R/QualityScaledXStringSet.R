### =========================================================================
### QualityScaledXStringSet objects
### -------------------------------------------------------------------------
###

setClass("QualityScaledXStringSet",
    contains="XStringSet",
    representation(
        "VIRTUAL",
        quality = "XStringQuality"
    )
)

### QualityScaledXStringSet subclasses
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
### Validity.
###

.valid.QualityScaledXStringSet <- function(object)
{
    message <- NULL
    if (!(length(object@quality) %in% c(1, length(object))))
        message <- c(message, "'length(quality)' != 1 or length of 'XStringSet'")
    if (!all(nchar(object@quality) == 1 | nchar(object@quality) == nchar(object)))
        message <- c(message, "'nchar(quality)' must equal 1 or nchar of 'XStringSet'")
    message
}

setValidity("QualityScaledXStringSet",
    function(object)
    {
        problems <- .valid.QualityScaledXStringSet(object)
        if (is.null(problems)) TRUE else problems
    }
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
    if (!is(quality, "XStringQuality"))
        stop("'quality' must be of class 'XStringQuality'")
    if (!(length(quality) %in% c(1, length(x))))
        stop("'length(quality)' must equal 1 or 'length(x)'")
    if (!all(nchar(quality) == 1 | nchar(quality) == nchar(x)))
        stop("'nchar(quality)' must equal 1 or 'nchar(x)'")
    output <- as(x, paste("QualityScaled", class(x), sep=""))
    slot(output, "quality", check=FALSE) <- quality
    output
}

QualityScaledBStringSet <- function(x, quality) QualityScaledXStringSet(BStringSet(x), quality)
QualityScaledDNAStringSet <- function(x, quality) QualityScaledXStringSet(DNAStringSet(x), quality)
QualityScaledRNAStringSet <- function(x, quality) QualityScaledXStringSet(RNAStringSet(x), quality)
QualityScaledAAStringSet <- function(x, quality) QualityScaledXStringSet(AAStringSet(x), quality)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Inherited methods.
###

setMethod("narrow", "QualityScaledXStringSet",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        x@ranges <-
          narrow(x@ranges, start=start, end=end, width=width,
                 use.names=use.names)
        if (all(width(x@quality@ranges) > 1)) {
            x@quality@ranges <-
              narrow(x@quality@ranges, start=start, end=end, width=width,
                     use.names=use.names)
        }
        x
    }
)


setMethod("append", c("QualityScaledXStringSet", "QualityScaledXStringSet"),
    function(x, values, after=length(x))
    {
        QualityScaledXStringSet(append(as(x, "XStringSet"), as(values, "XStringSet"), after=after),
                                as(append(x@quality, values@quality, after=after), class(x@quality)))
    }
)


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

.safe.subset.XStringSet <- function(x, i)
{
    if (length(x) == 1) x else x[i]
}

setMethod("[", "QualityScaledXStringSet",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i) || (is.logical(i) && all(i)))
            return(x)
        if (is.logical(i))
            i <- which(i)
        if (!is.numeric(i) || any(is.na(i)))
            stop("invalid subsetting")
        if (any(i < 1) || any(i > length(x)))
            stop("subscript out of bounds")
        slot(x, "quality", check=FALSE) <- .safe.subset.XStringSet(quality(x), i)
        callNextMethod(x, i)
    }
)

