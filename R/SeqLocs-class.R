### =========================================================================
### The SeqLocs class
### -------------------------------------------------------------------------
###
### The SeqLocs class is the basic container for storing a set of sequence
### locations defined by their start/nchar.
###

setClass("SeqLocs",
    representation(
        ## The 'locs' slot must be a data frame containing a "valid set of
        ## start/nchar locations" i.e. a data frame with a "start" and
        ## a "nchar" column, both columns being integer vectors (eventually
        ## of length 0) with no NAs and such that all(start >= 1) and
        ## all(nchar >= 0) are TRUE.
        ## Additionally this data frame can have an extra "names" column
        ## (a character vector).
        ## See the "initialize" method below for more details.
        locs="data.frame"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
###

setMethod("initialize", "SeqLocs",
    function(.Object, start=integer(0), nchar=integer(0), names=NULL, check=TRUE)
    {
        if (check) {
            if (!is.integer(start) || any(is.na(start)))
                stop("'start' must be an integer vector with no NAs")
            if (!is.integer(nchar) || any(is.na(nchar)))
                stop("'nchar' must be an integer vector with no NAs")
            if (length(start) != length(nchar))
                stop("'start' and 'nchar' must have the same length")
            if (!all(start >= 1L))
                stop("all values in 'start' must be >= 1")
            if (!all(nchar >= 0L))
                stop("all values in 'nchar' must be >= 0")
        }
        if (is.null(names)) {
            slot(.Object, "locs", check=FALSE) <- data.frame(start=start, nchar=nchar,
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
                stop("'names' must have the same length as 'start' and 'nchar'")
        }
        slot(.Object, "locs", check=FALSE) <- data.frame(start=start, nchar=nchar, names=names,
                                                         check.names=FALSE,
                                                         stringsAsFactors=FALSE)
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("length", "SeqLocs", function(x) nrow(x@locs))

### The "start" and "end" generics are defined in the stats package.
setMethod("start", "SeqLocs", function(x, ...) x@locs$start)
setMethod("nchar", "SeqLocs", function(x, type="chars", allowNA=FALSE) x@locs$nchar)
### Note that when nchar(x)[i] is 0, the end(x)[i] is start(x)[i] - 1
setMethod("end", "SeqLocs", function(x, ...) (x@locs$start + x@locs$nchar -1L))

setMethod("names", "SeqLocs", function(x) x@locs$names)

setReplaceMethod("names", "SeqLocs",
    function(x, value)
    {
        if (is.null(value)) {
            x@locs$names <- NULL
            return(x)
        }
        if (!is.character(value))
            stop("'value' must be a character vector (or NULL)")
        if (length(value) > length(x))
            stop("new 'names' vector has more elements than 'x'")
        length(value) <- length(x)
        x@locs$names <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

setMethod("show", "SeqLocs", function(object) show(object@locs))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Supported 'i' types: numeric vector, logical vector, NULL and missing.
setMethod("[", "SeqLocs",
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
        x@locs <- x@locs[i, , drop=FALSE]
        x
    }
)

setReplaceMethod("[", "SeqLocs",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a ", sQuote(class(x)), " object")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other methods.
###

setMethod("as.data.frame", "SeqLocs",
    function(x, row.names=NULL, optional=FALSE, ...)
        as.data.frame(x@locs, row.names=row.names, optional=optional, ...)
)

setMethod("as.matrix", "SeqLocs",
    function(x, ...)
    {
        ans <- as.matrix(x@locs[ , c("start", "nchar")], ...)
        rownames(ans) <- x@locs$names
        ans
    }
)

