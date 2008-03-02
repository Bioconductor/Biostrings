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
            inters <- data.frame(start=start, width=width,
                                 check.names=FALSE, stringsAsFactors=FALSE)
            slot(.Object, "inters", check=FALSE) <- inters
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
### The "intToIntervals" function.
###

intToIntervals <- function(x, use.names=TRUE)
{
    if (!is.numeric(x))
        stop("'x' must be an integer vector")
    if (!is.integer(x))
        x <- as.integer(x)
    use.names <- normalize.use.names(use.names)
    if (use.names) ans_names <- names(x) else ans_names <- NULL
    new("IntIntervals", start=rep.int(1L, length(x)), width=x,
        names=ans_names, check=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "intToAdjacentIntervals" function.
###

intToAdjacentIntervals <- function(x, use.names=TRUE)
{
    if (!is.numeric(x))
        stop("'x' must be an integer vector")
    if (!is.integer(x))
        x <- as.integer(x)
    use.names <- normalize.use.names(use.names)
    ans_start <- .Call("int_to_adjacent_intervals", x, PACKAGE="Biostrings")
    if (use.names) ans_names <- names(x) else ans_names <- NULL
    new("IntIntervals", start=ans_start, width=x,
        names=ans_names, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "restrict" function.
###

restrict <- function(x, start, end, use.names=TRUE)
{
    if (!is(x, "IntIntervals"))
        stop("'x' must be an IntIntervals object")
    if (!isSingleNumber(start))
        stop("'start' must be a single integer")
    if (!is.integer(start))
        start <- as.integer(start)
    if (!isSingleNumber(end))
        stop("'end' must be a single integer")
    if (!is.integer(end))
        end <- as.integer(end)
    if (start > end + 1L)
        stop("'start' must be <= 'end + 1'")
    use.names <- normalize.use.names(use.names)

    ans_start <- start(x)
    ans_end <- end(x)

    ## "fix" ans_start
    far_too_right <- end < ans_start
    ans_start[far_too_right] <- end + 1L
    too_left <- ans_start < start
    ans_start[too_left] <- start

    ## "fix" ans_end
    far_too_left <- ans_end < start
    ans_end[far_too_left] <- start - 1L
    too_right <- end < ans_end
    ans_end[too_right] <- end

    ans_width <- ans_end - ans_start + 1L
    if (use.names) ans_names <- names(x) else ans_names <- NULL
    new("IntIntervals", start=ans_start, width=ans_width,
        names=ans_names, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "narrow" function.
###

narrow <- function(x, start=NA, end=NA, width=NA, use.names=TRUE)
{
    if (is.numeric(x)) {
        x <- intToIntervals(x, use.names=use.names)
        ans_names <- names(x)
    } else {
        if (!is(x, "IntIntervals"))
            stop("'x' must be an IntIntervals object (or a numeric vector)")
        use.names <- normalize.use.names(use.names)
        if (use.names) ans_names <- names(x) else ans_names <- NULL
    }
    if (!isSingleNumberOrNA(start))
        stop("'start' must be a single integer or NA")
    if (!is.integer(start))
        start <- as.integer(start)
    if (!isSingleNumberOrNA(end))
        stop("'end' must be a single integer or NA")
    if (!is.integer(end))
        end <- as.integer(end)
    if (!isSingleNumberOrNA(width))
        stop("'width' must be a single integer or NA")
    if (!is.integer(width))
        width <- as.integer(width)
    C_ans <- .Call("narrow_IntIntervals", x, start, end, width, PACKAGE="Biostrings")
    new("IntIntervals", start=C_ans$start, width=C_ans$width,
        names=ans_names, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Normalization.
###

normalize <- function(x, with.inframe.attrib=FALSE)
{
    if (!is(x, "IntIntervals"))
        stop("'x' must be an IntIntervals object")
    if (!isTRUEorFALSE(with.inframe.attrib))
        stop("'with.inframe.attrib' must be 'TRUE' or 'FALSE'")
    C_ans <- .Call("normalize_IntIntervals", x, with.inframe.attrib, PACKAGE="Biostrings")
    ans <- new("IntIntervals", start=C_ans$start, width=C_ans$width, check=FALSE)
    if (with.inframe.attrib) {
        inframe <- new("IntIntervals", start=C_ans$inframe.start, width=width(x), check=FALSE)
        attr(ans, "inframe") <- inframe
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Deprecated methods.
###

setGeneric("first", function(x) standardGeneric("first"))
setMethod("first", "IntIntervals", function(x) {.Deprecated("start"); start(x)})
setGeneric("last", function(x) standardGeneric("last"))
setMethod("last", "IntIntervals", function(x) {.Deprecated("end"); end(x)})

