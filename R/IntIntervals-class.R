### =========================================================================
### IntIntervals objects
### -------------------------------------------------------------------------
###
### The IntIntervals class is a simple container for storing a set of integer
### intervals.
###

setClass("IntIntervals",
    representation(
        ## See the "initialize" method below for more details.
        inters="data.frame"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
###
### Both 'start' and 'width' must be integer vectors of equal length
### (eventually 0) with no NAs and such that all(width >= 0) is TRUE.
### 'names' must be NULL or a character vector of the same length as 'start'
### (or 'width').
###

setMethod("initialize", "IntIntervals",
    function(.Object, start=integer(0), width=integer(0), names=NULL, check=TRUE)
    {
        if (check) {
            if (!is.integer(start) || any(is.na(start)))
                stop("'start' must be an integer vector with no NAs")
            if (!is.integer(width) || any(is.na(width)))
                stop("'width' must be an integer vector with no NAs")
            if (length(start) != length(width))
                stop("'start' and 'width' must have the same length")
            if (!all(width >= 0L))
                stop("all values in 'width' must be >= 0")
        }
        if (is.null(names)) {
            slot(.Object, "inters", check=FALSE) <- data.frame(start=start, width=width,
                                                             check.names=FALSE,
                                                             stringsAsFactors=FALSE)
            return(.Object)
        }
        if (check) {
            if (!is.character(names))
                stop("'names' must be a character vector (or NULL)")
            # Disabled for now. Forbidding NAs or empty strings would not be
            # consistent with the "names<-" method that currently allows the
            # user to stick this kind of values into the 'names' slot!
            #if (any(names %in% c(NA, "")))
            #    stop("'names' cannot contain NAs or empty strings")
            if (length(names) != length(start))
                stop("'names' must have the same length as 'start' and 'width'")
        }
        inters <- data.frame(start=start, width=width, names=names,
                             check.names=FALSE, stringsAsFactors=FALSE)
        slot(.Object, "inters", check=FALSE) <- inters
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("length", "IntIntervals", function(x) nrow(x@inters))

### The substr() function uses 'start' and 'stop'.
### The substring() function uses 'first' and 'last'.
### We use 'start' and 'end'.
### Note that the "start" and "end" generic are defined in the stats package.
setMethod("start", "IntIntervals", function(x, ...) x@inters$start)

setGeneric("width", function(x) standardGeneric("width"))
setMethod("width", "IntIntervals", function(x) x@inters$width)

### Note that when width(x)[i] is 0, then end(x)[i] is start(x)[i] - 1
setMethod("end", "IntIntervals", function(x, ...) {start(x) + width(x) - 1L})

setMethod("names", "IntIntervals", function(x) x@inters$names)

setReplaceMethod("names", "IntIntervals",
    function(x, value)
    {
        if (is.null(value)) {
            x@inters$names <- NULL
            return(x)
        }
        if (!is.character(value))
            stop("'value' must be a character vector (or NULL)")
        if (length(value) > length(x))
            stop("too many names")
        length(value) <- length(x)
        x@inters$names <- value
        x
    }
)

### "desc" is an alias for "names". It may be deprecated soon...
setGeneric("desc", function(x) standardGeneric("desc"))
setMethod("desc", "ANY", function(x) names(x))
setGeneric("desc<-", signature="x", function(x, value) standardGeneric("desc<-"))
setReplaceMethod("desc", "ANY", function(x, value) `names<-`(x, value))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

setMethod("show", "IntIntervals",
    function(object)
        show(object@inters)
        #show(data.frame(object@inters, end=end(object), check.names=FALSE))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Supported 'i' types: numeric vector, logical vector, NULL and missing.
setMethod("[", "IntIntervals",
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
        x@inters <- x@inters[i, , drop=FALSE]
        x
    }
)

setReplaceMethod("[", "IntIntervals",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a ", sQuote(class(x)), " object")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other methods.
###

setMethod("as.data.frame", "IntIntervals",
    function(x, row.names=NULL, optional=FALSE, ...)
        as.data.frame(x@inters, row.names=row.names, optional=optional, ...)
)

setMethod("as.matrix", "IntIntervals",
    function(x, ...)
    {
        ans <- as.matrix(x@inters[ , c("start", "width")], ...)
        rownames(ans) <- x@inters$names
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Deprecated methods.
###

setGeneric("first", function(x) standardGeneric("first"))
setMethod("first", "IntIntervals", function(x) {.Deprecated("start"); start(x)})
setGeneric("last", function(x) standardGeneric("last"))
setMethod("last", "IntIntervals", function(x) {.Deprecated("end"); end(x)})

