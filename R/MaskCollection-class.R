### =========================================================================
### MaskCollection objects
### -------------------------------------------------------------------------

setClass("MaskCollection",
    representation(
        nir_list="list",   # a list of NormalIRanges objects
        width="integer",
        active="logical",
        NAMES="character"  # R doesn't like @names !!
    ),
    prototype(
        nir_list=list(),
        width=0L,
        active=logical(0),
        NAMES=as.character(NA)
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setGeneric("nir_list", function(x) standardGeneric("nir_list"))
setMethod("nir_list", "MaskCollection", function(x) x@nir_list)

setMethod("length", "MaskCollection", function(x) length(nir_list(x)))

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
        if (is.null(value)) {
            x@NAMES <- as.character(NA)
            return(x)
        }
        if (!is.character(value))
            stop("'value' must be NULL or a character vector")
        ii <- is.na(value)
        if (any(ii))
            value[ii] <- ""
        if (length(value) > length(x))
            stop("too many names")
        if (length(value) < length(x))
            value <- c(value, character(length(x) - length(value)))
        x@NAMES <- value
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

.valid.MaskCollection.nir_list <- function(object)
{
    if (!is.list(nir_list(object))
     || !all(sapply(nir_list(object), function(nir) is(nir, "NormalIRanges"))))
        return("the 'nir_list' slot must contain a list of NormalIRanges objects")
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
      .valid.MaskCollection.nir_list(object),
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
    new("MaskCollection", nir_list=list(nir),
                          width=numeric2integer(mask.width),
                          active=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "isEmpty" methods.
###

setMethod("isEmpty", "MaskCollection", function(x) sapply(nir_list(x), isEmpty))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "max" and "min" methods.
###

setMethod("max", "MaskCollection",
    function(x, ..., na.rm)
    {
        if (length(x) == 0)
            return(integer(0))
        sapply(nir_list(x), max)
    }
)

setMethod("min", "MaskCollection",
    function(x, ..., na.rm)
    {
        if (length(x) == 0)
            return(integer(0))
        sapply(nir_list(x), min)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "maskedwidth" and "maskedratio" generics and methods.
###

setGeneric("maskedwidth", function(x) standardGeneric("maskedwidth"))
setMethod("maskedwidth", "MaskCollection",
    function(x)
    {
        nir_list <- nir_list(x)
        if (length(nir_list) == 0)
            integer(0)
        else
            sapply(nir_list, function(nir) sum(width(nir)))
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
        nir_list(x)[[i]]
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
        slot(x, "nir_list", check=FALSE) <- nir_list(x)[i]
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

setMethod("append", "MaskCollection",
    function(x, values, after=length(x))
    {
        if (!is(values, "MaskCollection"))
            stop("'values' must be a MaskCollection object")
        if (width(values) != width(x))
            stop("'x' and 'values' must have the same width")
        if (!isSingleNumber(after))
            stop("'after' must be a single number")
        ans_nir_list <- append(nir_list(x), nir_list(values), after=after)
        ans_active <- append(active(x), active(values), after=after)
        nm1 <- names(x)
        nm2 <- names(values)
        if (is.null(nm1) && is.null(nm2)) {
            ans_NAMES <- as.character(NA)
        } else {
            if (is.null(nm1))
                nm1 <- rep.int("", length(x))
            if (is.null(nm2))
                nm2 <- rep.int("", length(values))
            ans_NAMES <- append(nm1, nm2, after=after)
        }
        ## This transformation must be atomic.
        x@nir_list <- ans_nir_list
        x@active <- ans_active
        x@NAMES <- ans_NAMES
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some useful endomorphisms: "narrow", "reduce" and "gaps".
###

setMethod("narrow", "MaskCollection",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        limits <- new("IRanges", start=1L, width=width(x), check=FALSE)
        limits <- narrow(limits, start=start, end=end, width=width)
        start <- start(limits)
        end <- end(limits)
        width <- width(limits)
        x@nir_list <- lapply(nir_list(x),
            function(nir) shift(restrict(nir, start=start, end=end), 1L - start(limits))
        )
        x@width <- width
        if (!normargUseNames(use.names))
            x@NAMES <- as.character(NA)
        x
    }
)

### 'with.inframe.attrib' is ignored.
setMethod("reduce", "MaskCollection",
    function(x, with.inframe.attrib=FALSE)
    {
        keep_it <- active(x)
        if (!all(keep_it))
            x <- x[keep_it]
        if (length(x) == 1)
            return(x)
        nir_list <- nir_list(x)
        if (length(nir_list) == 0) {
            nir1 <- newEmptyNormalIRanges()
        } else {
            start1 <- unlist(lapply(nir_list, start))
            width1 <- unlist(lapply(nir_list, width))
            ranges <- new("IRanges", start=start1, width=width1, check=FALSE)
            nir1 <- toNormalIRanges(ranges)
        }
        ## This transformation must be atomic.
        x@nir_list <- list(nir1)
        x@active <- TRUE
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
        x@nir_list <- lapply(nir_list(x),
            function(nir) gaps(nir, start=start, end=end)
        )
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
            margin <- format("", width=nchar(as.character(lx)))
            cat("all masks together:\n")
            mask0 <- reduce(`active<-`(x, TRUE))
            frame <- data.frame(maskedwidth=maskedwidth(mask0),
                                maskedratio=maskedratio(mask0),
                                check.names=FALSE)
            row.names(frame) <- margin
            show(frame)
            if (sum(active(x)) < lx) {
                cat("all active masks together:\n")
                mask1 <- reduce(x)
                frame <- data.frame(maskedwidth=maskedwidth(mask1),
                                    maskedratio=maskedratio(mask1),
                                    check.names=FALSE)
                row.names(frame) <- margin
                show(frame)
            }
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

