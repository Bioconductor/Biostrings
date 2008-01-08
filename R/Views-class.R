### =========================================================================
### The Views class
### -------------------------------------------------------------------------
###
### The Views class is the basic container for storing a set of start/end
### locations.
###

setClass("Views",
    representation(
        ## The 'views' slot must be a data frame containing a "valid set of
        ## views" i.e. a data frame with a "start" and an "end" column, both
        ## columns being integer vectors (eventually of length 0) with no NAs
        ## and such that all(start <= end) is TRUE.
        ## Additionally this data frame can also have a "desc" column
        ## (a character vector) that gives a short description of each view.
        ## See the "initialize" method below for more details.
        views="data.frame"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
###

setMethod("initialize", "Views",
    function(.Object, start, end, desc=NULL, check.data=TRUE)
    {
        if (check.data) {
            if (!is.integer(start) || any(is.na(start)))
                stop("'start' must be an integer vector with no NAs")
            if (!is.integer(end) || any(is.na(end)))
                stop("'end' must be an integer vector with no NAs")
            if (length(start) != length(end))
                stop("'start' and 'end' must have the same length")
            if (!all(start <= end))
                stop("'start' and 'end' must verify 'all(start <= end)'")
        }
        if (is.null(desc)) {
            .Object@views <- data.frame(start=start, end=end, check.names=FALSE, stringsAsFactors=FALSE)
            return(.Object)
        }
        if (check.data) {
            if (!is.character(desc))
                stop("'desc' can only be NULL or a character vector")
            # Disabled for now. Forbidding NAs or empty strings would not be
            # consistent with the "desc<-" method that currently allows the
            # user to stick this kind of values into the desc slot!
            #if (any(desc %in% c(NA, "")))
            #    stop("'desc' cannot contain NAs or empty strings")
            if (length(desc) != length(start))
                stop("'desc' must have the same length as 'start' and 'end'")
        }
        .Object@views <- data.frame(start=start, end=end, desc=desc, check.names=FALSE, stringsAsFactors=FALSE)
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("length", "Views", function(x) nrow(x@views))

### The substr() function uses 'start' and 'stop'.
### The substring() function uses 'first' and 'last'.
### We use 'start' and 'end'. Note that the "start" and "end" generic functions
### are defined in the stats package.
setMethod("start", "Views", function(x, ...) x@views$start)
setMethod("end", "Views", function(x, ...) x@views$end)

setGeneric("first", function(x) standardGeneric("first"))
setMethod("first", "Views", function(x) {.Deprecated("start"); start(x)})
setGeneric("last", function(x) standardGeneric("last"))
setMethod("last", "Views", function(x) {.Deprecated("end"); end(x)})

setGeneric("width", function(x) standardGeneric("width"))
setMethod("width", "Views", function(x) end(x) - start(x) + 1L)

setGeneric("desc", function(x) standardGeneric("desc"))
setMethod("desc", "Views", function(x) x@views$desc)

setGeneric("desc<-", signature="x", function(x, value) standardGeneric("desc<-"))
setReplaceMethod("desc", "Views",
    function(x, value)
    {
        if (is.null(value)) {
            x@views$desc <- NULL
            return(x)
        }
        if (!is.character(value))
            stop("'value' must be NULL or a character vector")
        if (length(value) > length(x))
            stop("new 'desc' vector has more elements than the number of views")
        length(value) <- length(x)
        x@views$desc <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

setMethod("show", "Views", function(object) show(object@views))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Supported 'i' types: numeric vector, logical vector, NULL and missing.
setMethod("[", "Views",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            return(x)
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
        x@views <- x@views[i, , drop=FALSE]
        x
    }
)

setReplaceMethod("[", "Views",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a ", sQuote(class(x)), " object")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other methods.
###

setMethod("as.data.frame", "Views",
    function(x, row.names=NULL, optional=FALSE, ...)
        as.data.frame(x@views, row.names=row.names, optional=optional, ...)
)

setMethod("as.matrix", "Views",
    function(x, ...)
    {
        ans <- as.matrix(x@views[ , c("start", "end")])
        rownames(ans) <- x@views$desc
        ans
    }
)

