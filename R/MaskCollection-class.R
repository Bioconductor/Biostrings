### =========================================================================
### MaskCollection objects
### -------------------------------------------------------------------------

setClass("MaskCollection",
    representation(
        masks="list",      # a list of NormalIRanges objects
        width="integer",
        NAMES="character"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setGeneric("masks", function(x) standardGeneric("masks"))
setMethod("masks", "MaskCollection", function(x) x@masks)

setMethod("length", "MaskCollection", function(x) length(masks(x)))

setMethod("width", "MaskCollection", function(x) x@width)

setMethod("names", "MaskCollection",
    function(x)
        if (length(x@NAMES) == 1 && is.na(x@NAMES)) NULL else x@NAMES
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "isEmpty" methods.
###

setMethod("isEmpty", "MaskCollection", function(x) sapply(masks(x), isEmpty))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "max" and "min" methods.
###

setMethod("max", "MaskCollection",
    function(x, ..., na.rm)
    {
        if (length(x) == 0)
            return(integer(0))
        sapply(masks(x), max)
    }
)

setMethod("min", "MaskCollection",
    function(x, ..., na.rm)
    {
        if (length(x) == 0)
            return(integer(0))
        sapply(masks(x), min)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.MaskCollection.width <- function(object)
{
    if (!isSingleInteger(width(object)) || width(object) < 0)
        return("the width of the collection must be a single non-negative integer")
    NULL
}

.valid.MaskCollection.masks <- function(object)
{
    if (!is.list(masks(object))
     || !all(sapply(masks(object), function(x) is(x, "NormalIRanges"))))
        return("the 'masks' slot must contain a list of NormalIRanges objects")
    if (!all(1 <= min(object)) || !all(max(object) <= width(object)))
        return("the min and max of the masks must be >= 1 and <= width of the collection")
    NULL
}

.valid.MaskCollection.names <- function(object)
{
    if (!is.character(object@NAMES))
        return("the 'NAMES' slot must contain a character vector")
    if (is.null(names(object)))
        return(NULL)
    if (any(is.na(names(object))))
        return("the names must be non-NA strings")
    if (length(names(object)) != length(object))
        return("number of names and number of elements differ")
    NULL
}

.valid.MaskCollection <- function(object)
{
    c(.valid.MaskCollection.width(object),
      .valid.MaskCollection.masks(object),
      .valid.MaskCollection.names(object))
}

setValidity("MaskCollection",
    function(object)
    {
        problems <- .valid.MaskCollection(object)
        if (is.null(problems)) TRUE else problems
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "maskedwidth" and "maskedratio" generics and methods.
###

setGeneric("maskedwidth", function(x) standardGeneric("maskedwidth"))
setMethod("maskedwidth", "MaskCollection",
    function(x)
    {
        masks <- masks(x)
        if (length(masks) == 0)
            integer(0)
        else
            sapply(masks, function(mask) sum(width(mask)))
    }
)

setGeneric("maskedratio", function(x) standardGeneric("maskedratio"))
setMethod("maskedratio", "MaskCollection", function(x) maskedwidth(x) / width(x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "reduce" method.
###

setMethod("reduce", "MaskCollection",
    function(x, with.inframe.attrib=FALSE)
    {
        masks <- masks(x)
        if (length(masks) == 0) {
            mask1 <- newEmptyNormalIRanges()
        } else {
            start1 <- unlist(lapply(masks, start))
            width1 <- unlist(lapply(masks, width))
            mask1 <- toNormalIRanges(new("IRanges", start=start1, width=width1, check=FALSE))
        }
        x@masks <- list(mask1)
        x@NAMES <- as.character(NA)
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

setMethod("show", "MaskCollection",
    function(object)
    {
        lo <- length(object)
        cat("  A ", class(object), " instance of length ", lo,
            " and width ", width(object), "\n", sep="")
        cat("masks:")
        if (lo == 0) {
            cat(" NONE\n")
        } else {
            cat("\n")
            frame <- data.frame(maskedwidth=maskedwidth(object),
                                maskedratio=maskedratio(object),
                                check.names=FALSE)
            frame$names <- names(object)
            show(frame)
        }
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Supported 'i' types: numeric vector, logical vector, NULL and missing.
setMethod("[", "MaskCollection",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            return(x)
        if (is.character(i))
            stop("cannot subset a ", class(x), " object by names")
        lx <- length(x)
        if (is.numeric(i)) {
            if (any(i < -lx) || any(i > lx))
                stop("subscript out of bounds")
        } else if (is.logical(i)) {
            if (length(i) > lx)
                stop("subscript out of bounds")
        } else if (!is.null(i)) {
            stop("invalid subscript type")
        }
        slot(x, "masks", check=FALSE) <- masks(x)[i]
        if (!is.null(names(x)))
            slot(x, "NAMES", check=FALSE) <- names(x)[i]
        x
    }
)

setReplaceMethod("[", "MaskCollection",
    function(x, i, j,..., value)
        stop("attempt to modify the value of a ", class(x), " instance")
)

