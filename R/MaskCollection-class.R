### =========================================================================
### MaskCollection objects
### -------------------------------------------------------------------------

setClass("MaskCollection",
    representation(
        nirlist="list",    # a list of NormalIRanges objects
        width="integer",
        NAMES="character"
    ),
    prototype(
        nirlist=list(),
        width=0L,
        NAMES=as.character(NA)
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setGeneric("nirlist", function(x) standardGeneric("nirlist"))
setMethod("nirlist", "MaskCollection", function(x) x@nirlist)

setMethod("length", "MaskCollection", function(x) length(nirlist(x)))

setMethod("width", "MaskCollection", function(x) x@width)

setMethod("names", "MaskCollection",
    function(x)
        if (length(x@NAMES) == 1 && is.na(x@NAMES)) NULL else x@NAMES
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

.valid.MaskCollection.nirlist <- function(object)
{
    if (!is.list(nirlist(object))
     || !all(sapply(nirlist(object), function(x) is(x, "NormalIRanges"))))
        return("the 'nirlist' slot must contain a list of NormalIRanges objects")
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
      .valid.MaskCollection.nirlist(object),
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
### The "isEmpty" methods.
###

setMethod("isEmpty", "MaskCollection", function(x) sapply(nirlist(x), isEmpty))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "max" and "min" methods.
###

setMethod("max", "MaskCollection",
    function(x, ..., na.rm)
    {
        if (length(x) == 0)
            return(integer(0))
        sapply(nirlist(x), max)
    }
)

setMethod("min", "MaskCollection",
    function(x, ..., na.rm)
    {
        if (length(x) == 0)
            return(integer(0))
        sapply(nirlist(x), min)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "maskedwidth" and "maskedratio" generics and methods.
###

setGeneric("maskedwidth", function(x) standardGeneric("maskedwidth"))
setMethod("maskedwidth", "MaskCollection",
    function(x)
    {
        nirlist <- nirlist(x)
        if (length(nirlist) == 0)
            integer(0)
        else
            sapply(nirlist, function(mask) sum(width(mask)))
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
        nirlist <- nirlist(x)
        if (length(nirlist) == 0) {
            mask1 <- newEmptyNormalIRanges()
        } else {
            start1 <- unlist(lapply(nirlist, start))
            width1 <- unlist(lapply(nirlist, width))
            mask1 <- toNormalIRanges(new("IRanges", start=start1, width=width1, check=FALSE))
        }
        x@nirlist <- list(mask1)
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
            if (lo >= 2) {
                cat("reduced mask (obtained with the 'reduce' method):\n")
                mask0 <- reduce(object)
                frame <- data.frame(maskedwidth=maskedwidth(mask0),
                                    maskedratio=maskedratio(mask0),
                                    check.names=FALSE)
                show(frame)
            }
        }
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Extract the i-th element of a MaskCollection object as a NormalIRanges
### object.
### Supported 'i' types: numeric vector of length 1.
setMethod("[[", "MaskCollection",
    function(x, i, j, ...)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            stop("subscript is missing")
        if (is.character(i))
            stop("cannot subset a ", class(x), " object by names")
        if (!is.numeric(i))
            stop("invalid subscript type")
        if (length(i) < 1L)
            stop("attempt to select less than one element")
        if (length(i) > 1L)
            stop("attempt to select more than one element")
        if (is.na(i))
            stop("subscript cannot be NA")
        if (i < 1L || i > length(x))
            stop("subscript out of bounds")
        nirlist(x)[[i]]
    }
)

setReplaceMethod("[[", "MaskCollection",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a ", class(x), " instance")
    }
)

### Supported 'i' types: numeric vector, logical vector, NULL and missing.
setMethod("[", "MaskCollection",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            return(x)
        if (!is.atomic(i))
            stop("invalid subscript type")
        if (is.character(i))
            stop("cannot subset a ", class(x), " object by names")
        lx <- length(x)
        if (is.numeric(i)) {
            if (any(is.na(i)))
                stop("subscript contains NAs")
            if (any(i < -lx) || any(i > lx))
                stop("subscript out of bounds")
        } else if (is.logical(i)) {
            if (any(is.na(i)))
                stop("subscript contains NAs")
            if (length(i) > lx)
                stop("subscript out of bounds")
        } else if (!is.null(i)) {
            stop("invalid subscript type")
        }
        slot(x, "nirlist", check=FALSE) <- nirlist(x)[i]
        if (!is.null(names(x)))
            slot(x, "NAMES", check=FALSE) <- names(x)[i]
        x
    }
)

setReplaceMethod("[", "MaskCollection",
    function(x, i, j,..., value)
        stop("attempt to modify the value of a ", class(x), " instance")
)

