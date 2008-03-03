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

### "desc" is an alias for "names". It might be deprecated soon...
setGeneric("desc", function(x) standardGeneric("desc"))
setMethod("desc", "ANY", function(x) names(x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###
### Both 'start' and 'width' must be integer vectors of equal length
### (eventually 0) with no NAs and such that all(width >= 0) is TRUE.
### 'names' must be NULL or a character vector of the same length as 'start'
### (or 'width').
###

.valid.IntIntervals.start <- function(object)
{
    if (!is.integer(start(object)) || any(is.na(start(object))))
        return("the starts must be non-NA integers")
    ## This is already enforced by the fact that 'start(object)' and
    ## 'width(object)' are currently stored in the same data frame but we
    ## check their lengths anyway because there is no guarantee that this will
    ## not change in the future (e.g. they could be moved from this data.frame
    ## to separate slots).
    if (length(start(object)) != length(width(object)))
        return("number of starts and number of widths differ")
    NULL
}

.valid.IntIntervals.width <- function(object)
{
    if (!is.integer(width(object)) || any(is.na(width(object))))
        return("the widths must be non-NA integers")
    ## See comment in .valid.IntIntervals.start() in above...
    if (length(start(object)) != length(width(object)))
        return("number of starts and number of widths differ")
    if (length(width(object)) != 0 && min(width(object)) < 0L)
        return("negative widths are not allowed")
    NULL
}

.valid.IntIntervals.names <- function(object)
{
    if (is.null(names(object)))
        return(NULL)
    if (!is.character(names(object)))
        return("the names must a character vector (or the NULL value)")
    # Disabled for now. Forbidding NAs or empty strings would not be
    # consistent with the "names<-" method that currently allows the
    # user to stick this kind of values into the 'names' slot!
    #if (any(names(object) %in% c(NA, "")))
    #    return("names cannot be NAs or empty strings")
    if (length(names(object)) != length(object))
        return("number of names and number of elements differ")
    NULL
}

.valid.IntIntervals <- function(object)
{
    #cat("validating IntIntervals object of length", length(object), "...\n")
    c(.valid.IntIntervals.start(object),
      .valid.IntIntervals.width(object),
      .valid.IntIntervals.names(object))
}

setValidity("IntIntervals",
    function(object)
    {
        problems <- .valid.IntIntervals(object)
        if (is.null(problems)) TRUE else problems
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
###

.make.inters <- function(start, width, names)
{
    ## We could put the starts and the widths in the data frame and then add
    ## the "names" column to it but it might be slower (maybe it will copy the
    ## original data frame?). I've not tested this though...
    if (is.null(names))
        inters <- data.frame(start=start, width=width,
                             check.names=FALSE, stringsAsFactors=FALSE)
    else
        inters <- data.frame(start=start, width=width, names=names,
                             check.names=FALSE, stringsAsFactors=FALSE)
    inters
}

setMethod("initialize", "IntIntervals",
    function(.Object, start=integer(0), width=integer(0), names=NULL, check=TRUE)
    {
        inters <- .make.inters(start, width, names)
        slot(.Object, "inters", check=FALSE) <- inters
        if (check) {
            ## I found that using validObject() in "initialize" doesn't work
            ## properly (validation is called too many times and not in an
            ## order that makes sense to me...)
            #validObject(.Object)
            problems <- .valid.IntIntervals(.Object)
            if (!is.null(problems)) stop(paste(problems, collapse="\n       "))
        }
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Replacement methods.
###
### The rules are:
###   (1) changing the start preserves the width (so it changes the end)
###   (2) changing the width preserves the start (so it changes the end)
###   (3) changing the end preserves the width (so it changes the start)
###
### Note that we don't call validObject(x) after 'x' has been modified because
### we don't need to revalidate the entire object: validating the bits that
### have been touched is enough (and faster). However, because of this, if
### instances of derived classes must satisfy additional constraints, then some
### of the replacement methods below need to be overridden. See for example the
### "width<-" method for BStringViews objects (BStringViews-class.R file).
###

setGeneric("start<-", signature="x",
    function(x, check=TRUE, value) standardGeneric("start<-")
)

setReplaceMethod("start", "IntIntervals",
    function(x, check=TRUE, value)
    {
        x@inters$start <- value
        if (check) {
            problem <- .valid.IntIntervals.start(x)
            if (!is.null(problem)) stop(problem)
        }
        x
    }
)

setGeneric("width<-", signature="x",
    function(x, check=TRUE, value) standardGeneric("width<-")
)

setReplaceMethod("width", "IntIntervals",
    function(x, check=TRUE, value)
    {
        x@inters$width <- value
        if (check) {
            problem <- .valid.IntIntervals.width(x)
            if (!is.null(problem)) stop(problem)
        }
        x
    }
)

setGeneric("end<-", signature="x",
    function(x, check=TRUE, value) standardGeneric("end<-")
)

setReplaceMethod("end", "IntIntervals",
    function(x, check=TRUE, value)
    {
        start(x, check=check) <- value - width(x) + 1L
        x
    }
)

setReplaceMethod("names", "IntIntervals",
    function(x, value)
    {
        if (is.null(value)) {
            x@inters$names <- NULL
            return(x)
        }
        if (!is.character(value))
            stop("'value' must be NULL or a character vector")
        if (length(value) > length(x))
            stop("too many names")
        length(value) <- length(x)
        x@inters$names <- value
        x
    }
)

setGeneric("desc<-", signature="x", function(x, value) standardGeneric("desc<-"))
setReplaceMethod("desc", "ANY", function(x, value) `names<-`(x, value))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "update" method.
###
### This is a convenience method for combining multiple modifications in one
### single call.
###
### It must verify 2 important properties:
###   (1) update(x) must be identical to x (doesn't touch x at all)
###   (2) do.call("update", c(x, x@inters)) must be identical to x (it updates
###       x with its own content)
###

setMethod("update", "IntIntervals",
    function(object, ...)
    {
        args <- list(...)
        argnames <- names(args)
        if (length(args) != 0
            && (is.null(argnames) || any(argnames %in% c("", NA))))
            stop("all extra arguments must be named")
        valid_argnames <- c("start", "end", "width", "names", "check")
        if (!all(argnames %in% valid_argnames))
            stop("valid extra argument names are ",
                 paste("'", valid_argnames, "'", sep="", collapse=", "))
        if (any(duplicated(argnames)))
            stop("argument names must be unique")
        check <- args$check
        if (is.null(check)) check <- TRUE
        swe <- c("start", "end", "width")
        narg_in_swe <- sum(swe %in% argnames)
        if (narg_in_swe == 3)
            stop("only two of the ",
                 paste("'", swe, "'", sep="", collapse=", "),
                 " arguments can be specified")
        do_atomic_update <- narg_in_swe == 2 && (("names" %in% argnames)
                                                 || is.null(names(object)))
        if (do_atomic_update) {
            if ("end" %in% argnames) {
                if ("width" %in% argnames) {
                    width <- args$width
                    start <- args$end - width + 1L
                } else {
                    start <- args$start
                    width <- args$end - start + 1L
                }
            } else {
                start <- args$start
                width <- args$width
            }
            inters <- .make.inters(start, width, args$names)
            slot(object, "inters", check=FALSE) <- inters
            if (check) validObject(object)
            return(object)
        }
        if ("start" %in% argnames)
            start(object, check=check) <- args$start
        if ("end" %in% argnames)
            end(object, check=check) <- args$end
        if ("width" %in% argnames)
            width(object, check=check) <- args$width
        if ("names" %in% argnames) {
            ## "names<-" has no 'check' argument
            if (check)
                names(object) <- args$names
            else
                object@inters$names <- args$names
        }
        object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

setMethod("as.data.frame", "IntIntervals",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        ans <- data.frame(start=start(x),
                          end=end(x),
                          width=width(x),
                          check.names=FALSE,
                          stringsAsFactors=FALSE)
        ans$names <- names(x)
        ans
    }
)

setMethod("show", "IntIntervals",
    function(object) show(as.data.frame(object))
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
### The "as.list" method.
###

#Not sure we want this!
#setMethod("as.list", "IntIntervals", function(x, ...) x@inters)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other methods.
###

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

