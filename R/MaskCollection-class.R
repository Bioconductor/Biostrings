### =========================================================================
### MaskCollection objects
### -------------------------------------------------------------------------

setClass("MaskCollection",
    representation(
        nirlist="list",    # a list of NormalIRanges objects
        width="integer",
        active="logical",
        NAMES="character"  # R doesn't like @names !!
    ),
    prototype(
        nirlist=list(),
        width=0L,
        active=logical(0),
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

setGeneric("active", function(x) standardGeneric("active"))
setMethod("active", "MaskCollection", function(x) x@active)

setGeneric("active<-", signature="x",
    function(x, value) standardGeneric("active<-")
)
setReplaceMethod("active", "MaskCollection",
    function(x, value)
    {
        if (!is.logical(value) || any(is.na(value)))
            stop("'value' must be a logical vector with no NAs")
        x@active[] <- value
        x
    }
)

setMethod("names", "MaskCollection",
    function(x)
        if (length(x@NAMES) == 1 && is.na(x@NAMES)) NULL else x@NAMES
)

### The only replacement method for MaskCollection objects!
setReplaceMethod("names", "MaskCollection",
    function(x, value)
    {
        if (is.character(value)) {
            ii <- is.na(value)
            if (any(ii))
                value[ii] <- ""
        } else if (!is.null(value)) {
            stop("'value' must be NULL or a character vector")
        }
        unsafe.names(x) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.MaskCollection.width <- function(object)
{
    if (!isSingleInteger(width(object)) || width(object) < 1)
        return("the width of the collection must be a single non-negative integer")
    NULL
}

.valid.MaskCollection.nirlist <- function(object)
{
    if (!is.list(nirlist(object))
     || !all(sapply(nirlist(object), function(nir) is(nir, "NormalIRanges"))))
        return("the 'nirlist' slot must contain a list of NormalIRanges objects")
    if (!all(1 <= min(object)) || !all(max(object) <= width(object)))
        return("the min and max of the masks must be >= 1 and <= width of the collection")
    NULL
}

.valid.MaskCollection.active <- function(object)
{
    if (!is.logical(active(object)) || any(is.na(active(object))))
        return("the 'active' slot must be a logical vector with no NAs")
    if (length(active(object)) != length(object))
        return("the length of the 'active' slot differs from the length of the object")
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
      .valid.MaskCollection.active(object),
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
### The safe and user-friendly "Mask" constructor.
###

Mask <- function(mask.width, start=NULL, end=NULL, width=NULL)
{
    nir <- asNormalIRanges(IRanges(start=start, end=end, width=width))
    new("MaskCollection", nirlist=list(nir),
                          width=numeric2integer(mask.width),
                          active=TRUE)
}


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
            sapply(nirlist, function(nir) sum(width(nir)))
    }
)

setGeneric("maskedratio", function(x) standardGeneric("maskedratio"))
setMethod("maskedratio", "MaskCollection", function(x) maskedwidth(x) / width(x))


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
        slot(x, "active", check=FALSE) <- active(x)[i]
        if (!is.null(names(x)))
            slot(x, "NAMES", check=FALSE) <- names(x)[i]
        x
    }
)

setReplaceMethod("[", "MaskCollection",
    function(x, i, j,..., value)
        stop("attempt to modify the value of a ", class(x), " instance")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "append" generic and method.
###

#setGeneric("append", signature="x",
#    function(x, values, after=length(x)) standardGeneric("append")
#)

setMethod("append", "MaskCollection",
    function(x, values, after=length(x))
    {
        if (!is(values, "MaskCollection"))
            stop("'values' must be a MaskCollection object")
        if (width(values) != width(x))
            stop("'x' and 'values' must have the same width")
        if (!isSingleNumber(after))
            stop("'after' must be a single number")
        x@nirlist <- append(nirlist(x), nirlist(values), after=after)
        x@active <- append(active(x), active(values), after=after)
        nm1 <- names(x)
        nm2 <- names(values)
        if (is.null(nm1) && is.null(nm2))
            return(x)
        if (is.null(nm1))
            nm1 <- rep.int("", length(x))
        if (is.null(nm2))
            nm2 <- rep.int("", length(values))
        x@NAMES <- append(nm1, nm2, after=after)
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The transformation methods (endomorphisms) "narrow", "reduce" and "gaps".
###

setMethod("narrow", "MaskCollection",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        limits <- new("IRanges", start=1L, width=width(x), check=FALSE)
        limits <- narrow(limits, start=start, end=end, width=width)
        start <- start(limits)
        end <- end(limits)
        width <- width(limits)
        x@nirlist <- lapply(nirlist(x),
                            function(nir)
                                shift(restrict(nir, start=start, end=end), 1L - start(limits))
                     )
        x@width <- width
        if (!normalize.use.names(use.names))
            x@NAMES <- as.character(NA)
        x
    }
)

### 'with.inframe.attrib' is ignored.
setMethod("reduce", "MaskCollection",
    function(x, with.inframe.attrib=FALSE)
    {
        nirlist <- nirlist(x)
        if (length(nirlist) == 0) {
            nir1 <- newEmptyNormalIRanges()
        } else if (length(nirlist) == 1) {
            return(x)
        } else {
            start1 <- unlist(lapply(nirlist, start))
            width1 <- unlist(lapply(nirlist, width))
            nir1 <- toNormalIRanges(new("IRanges", start=start1, width=width1, check=FALSE))
        }
        x@nirlist <- list(nir1)
        x@NAMES <- as.character(NA)
        x
    }
)

### 'start' and 'end' are ignored.
setMethod("gaps", "MaskCollection",
    function(x, start=NA, end=NA)
    {
        start <- 1L
        end <- width(x)
        x@nirlist <- lapply(nirlist(x), function(nir) gaps(nir, start=start, end=end))
        x@NAMES <- as.character(NA)
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

### From a MaskCollection object to a NormalIRanges object.
setAs("MaskCollection", "NormalIRanges",
    function(from) reduce(from)[[1]]
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

MaskCollection.show_frame <- function(x)
{
    lx <- length(x)
    cat("masks:")
    if (lx == 0) {
        cat(" NONE\n")
    } else {
        cat("\n")
        frame <- data.frame(maskedwidth=maskedwidth(x),
                            maskedratio=maskedratio(x),
                            active=active(x),
                            check.names=FALSE)
        frame$names <- names(x)
        show(frame)
        if (lx >= 2) {
            cat("all masks together:\n")
            mask0 <- reduce(x)
            frame <- data.frame(maskedwidth=maskedwidth(mask0),
                                maskedratio=maskedratio(mask0),
                                check.names=FALSE)
            show(frame)
        }
    }
}

setMethod("show", "MaskCollection",
    function(object)
    {
        lo <- length(object)
        cat("  A ", class(object), " instance of length ", lo,
            " and width ", width(object), "\n", sep="")
        MaskCollection.show_frame(object)
    }
)

