### =========================================================================
### IRanges objects
### -------------------------------------------------------------------------
###
### The IRanges class is a simple container for storing a set of integer
### ranges.
###

setClass(".IRanges",
    representation(
        start="integer",
        width="integer",
        NAMES="character" # R doesn't like @names !!
    )
)

setClass("IRanges", contains=".IRanges")

### A NormalIRanges object is an IRanges object where the ranges are:
###   (a) of non-null width;
###   (b) not overlapping;
###   (c) not even adjacent (there must be a non-null gap between 2
###       consecutive ranges);
###   (d) ordered from left to right.
### If 'x' is an IRanges object of length >= 2, then 'x' is normal iff:
###   start(x)[i] <= end(x)[i] < start(x)[i+1] <= end(x)[i+1]
### for every 1 <= i < length(x).
### If length(x) == 1, then 'x' is normal iff width(x)[1] >= 1.
### If length(x) == 0, then 'x' is normal.
### Subsetting 'x' preserves normality.

setClass("NormalIRanges", contains="IRanges")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("length", ".IRanges", function(x) length(x@start))

### The substr() function uses 'start' and 'stop'.
### The substring() function uses 'first' and 'last'.
### We use 'start' and 'end'.
### Note that the "start" and "end" generic are defined in the stats package.
setMethod("start", ".IRanges", function(x, ...) x@start)

setGeneric("width", function(x) standardGeneric("width"))
setMethod("width", ".IRanges", function(x) x@width)

### Note that when width(x)[i] is 0, then end(x)[i] is start(x)[i] - 1
setMethod("end", ".IRanges", function(x, ...) {start(x) + width(x) - 1L})

setMethod("names", ".IRanges",
    function(x)
        if (length(x@NAMES) == 1 && is.na(x@NAMES)) NULL else x@NAMES
)

### "desc" is an alias for "names". It might be deprecated soon...
setGeneric("desc", function(x) standardGeneric("desc"))
setMethod("desc", "ANY", function(x) names(x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "is.normal" generic and method.
###

setGeneric("is.normal", function(x) standardGeneric("is.normal"))

setMethod("is.normal", ".IRanges",
    function(x)
    {
        all(width(x) >= 1) && (length(x) <= 1 || all(end(x)[-length(x)] < start(x)[-1]))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###
### Both 'start' and 'width' must be integer vectors of equal length
### (eventually 0) with no NAs and such that all(width >= 0) is TRUE.
### 'names' must be NULL or a character vector of the same length as 'start'
### (or 'width').
###

.valid.IRanges.start <- function(object)
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

.valid.IRanges.width <- function(object)
{
    if (!is.integer(width(object)) || any(is.na(width(object))))
        return("the widths must be non-NA integers")
    ## See comment in .valid.IRanges.start() above...
    if (length(start(object)) != length(width(object)))
        return("number of starts and number of widths differ")
    if (length(width(object)) != 0 && min(width(object)) < 0L)
        return("negative widths are not allowed")
    NULL
}

.valid.IRanges.names <- function(object)
{
    if (is.null(names(object)))
        return(NULL)
    if (!is.character(names(object)) || any(is.na(names(object))))
        return("the names must be non-NA strings")
    # Disabled for now. Forbidding NAs or empty strings would not be
    # consistent with the "names<-" method that currently allows the
    # user to stick this kind of values into the 'names' slot!
    #if (any(names(object) %in% c(NA, "")))
    #    return("names cannot be NAs or empty strings")
    if (length(names(object)) != length(object))
        return("number of names and number of elements differ")
    NULL
}

.valid.IRanges <- function(object)
{
    #cat("validating IRanges object of length", length(object), "...\n")
    c(.valid.IRanges.start(object),
      .valid.IRanges.width(object),
      .valid.IRanges.names(object))
}

setValidity(".IRanges",
    function(object)
    {
        problems <- .valid.IRanges(object)
        if (is.null(problems)) TRUE else problems
    }
)

.valid.NormalIRanges <- function(object)
{
    if (!is.normal(object))
        return("object is not normal")
    NULL
}

setValidity("NormalIRanges",
    function(object)
    {
        problems <- .valid.NormalIRanges(object)
        if (is.null(problems)) TRUE else problems
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization and coercion.
###

.numeric2integer <- function(x)
{
    if (is.numeric(x) && !is.integer(x)) as.integer(x) else x
}

.set.IRanges.slots <- function(object, start, width, names)
{
    slot(object, "start", check=FALSE) <- .numeric2integer(start)
    slot(object, "width", check=FALSE) <- .numeric2integer(width)
    slot(object, "NAMES", check=FALSE) <- if (is.null(names)) as.character(NA) else names
    object
}

setMethod("initialize", ".IRanges",
    function(.Object, start, width, names)
        .set.IRanges.slots(.Object, start, width, names)
)

setMethod("initialize", "IRanges",
    function(.Object, start=integer(0), width=integer(0), names=NULL, check=TRUE)
    {
        .Object <- callNextMethod(.Object, start, width, names)
        if (check) {
            ## I found that using validObject() in "initialize" doesn't work
            ## properly (validation is called too many times and not in an
            ## order that makes sense to me...)
            #validObject(.Object)
            problems <- .valid.IRanges(.Object)
            if (!is.null(problems)) stop(paste(problems, collapse="\n  "))
        }
        .Object
    }
)

setMethod("initialize", "NormalIRanges",
    function(.Object, start=integer(0), width=integer(0), names=NULL, check=TRUE)
    {
        .Object <- callNextMethod(.Object, start=start, width=width, names=names, check=check)
        if (check) {
            problems <- .valid.NormalIRanges(.Object)
            if (!is.null(problems)) stop(paste(problems, collapse="\n  "))
        }
        .Object
    }
)

### By default as(x, "NormalIRanges") would not check that the returned object
### is a valid NormalIRanges object.
setAs(".IRanges", "NormalIRanges",
    function(from)
    {
        class(from) <- "NormalIRanges"
        validObject(from)
        from
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

setMethod("as.data.frame", ".IRanges",
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

setMethod("show", ".IRanges",
    function(object) show(as.data.frame(object))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Supported 'i' types: numeric vector, logical vector, NULL and missing.
setMethod("[", ".IRanges",
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
        slot(x, "start", check=FALSE) <- start(x)[i]
        slot(x, "width", check=FALSE) <- width(x)[i]
        if (!is.null(names(x)))
            slot(x, "NAMES", check=FALSE) <- names(x)[i]
        x
    }
)

setReplaceMethod("[", ".IRanges",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a ", sQuote(class(x)), " object")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "duplicated" method.
###
### TODO: current implementation is very inefficient and needs some C help!
###

setMethod("duplicated", ".IRanges",
    function(x, incomparables=FALSE, ...)
    {
        duplicated(data.frame(start=start(x),
                              width=width(x),
                              check.names=FALSE,
                              stringsAsFactors=FALSE))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other methods.
###

setMethod("as.matrix", ".IRanges",
    function(x, ...)
        matrix(data=c(start(x), width(x)), ncol=2, dimnames=list(names(x), NULL))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Non-exported replacement functions for .IRanges objects.
###
### IMPORTANT: They do NOT check their argument 'x' and 'value' at all, not
### even whether 'value' is of type integer or not!
###
### The rules are:
###   (1) changing the start preserves the width (so it changes the end)
###   (2) changing the width preserves the start (so it changes the end)
###   (3) changing the end preserves the width (so it changes the start)
###

`.start<-` <- function(x, value)
{
    ## Use 'x@start[]' instead of just 'x@start' so 'value' is recycled
    x@start[] <- value
    x
}

`.width<-` <- function(x, value)
{
    ## Use 'x@width[]' instead of just 'x@width' so 'value' is recycled
    x@width[] <- value
    x
}

`.end<-` <- function(x, value)
{
    .start(x) <- value - width(x) + 1L
    x
}

`.names<-` <- function(x, value)
{
    if (is.null(value))
        x@NAMES <- as.character(NA)
    else {
        if (is.null(names(x)))
            x@NAMES <- character(length(x))
        ## Use 'x@NAMES[]' instead of just 'x@NAMES' so 'value' is recycled
        x@NAMES[] <- value
    }
    x
}

.update <- function(object, ...)
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
    sew <- c("start", "end", "width")
    narg_in_sew <- sum(sew %in% argnames)
    if (narg_in_sew == 3)
        stop("only two of the ",
             paste("'", sew, "'", sep="", collapse=", "),
             " arguments can be specified")
    do_atomic_update <- narg_in_sew == 2 && (is.null(names(object))
                                             || ("names" %in% argnames))
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
        return(.set.IRanges.slots(object, start, width, args$names))
    }
    if ("start" %in% argnames)
        .start(object) <- args$start
    if ("end" %in% argnames)
        .end(object) <- args$end
    if ("width" %in% argnames)
        .width(object) <- args$width
    if ("names" %in% argnames)
        .names(object) <- args$names
    object
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Exported replacement methods.
###
### See .start<-, .width<- and .end<- above for the rules.
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

### No method for .IRanges objects!
setReplaceMethod("start", "IRanges",
    function(x, check=TRUE, value)
    {
        .start(x) <- .numeric2integer(value)
        if (check) {
            problem <- .valid.IRanges.start(x)
            if (!is.null(problem)) stop(problem)
        }
        x
    }
)

setGeneric("width<-", signature="x",
    function(x, check=TRUE, value) standardGeneric("width<-")
)

### No method for .IRanges objects!
setReplaceMethod("width", "IRanges",
    function(x, check=TRUE, value)
    {
        .width(x) <- .numeric2integer(value)
        if (check) {
            problem <- .valid.IRanges.width(x)
            if (!is.null(problem)) stop(problem)
        }
        x
    }
)

setGeneric("end<-", signature="x",
    function(x, check=TRUE, value) standardGeneric("end<-")
)

### No method for .IRanges objects!
setReplaceMethod("end", "IRanges",
    function(x, check=TRUE, value)
    {
        .end(x) <- .numeric2integer(value)
        if (check) {
            problem <- .valid.IRanges.start(x)
            if (!is.null(problem)) stop(problem)
        }
        x
    }
)

### Yes, for .IRanges objects!
setReplaceMethod("names", ".IRanges",
    function(x, value)
    {
        if (!(is.null(value) || is.character(value)))
            stop("'value' must be NULL or a character vector")
        if (length(value) > length(x))
            stop("too many names")
        length(value) <- length(x)
        .names(x) <- value
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
###   (2) update(x, start=start(x), width=width(x), names=names(x))
###       must be identical to x too (but this time it updates x with its own content)
###

### No method for .IRanges objects!
setMethod("update", "IRanges",
    function(object, ...)
    {
        object <- .update(object, ...)
        check <- list(...)$check
        if (is.null(check))
            check <- TRUE
        if (check)
            validObject(object)
        object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Deprecated methods.
###

setGeneric("first", function(x) standardGeneric("first"))
setMethod("first", ".IRanges", function(x) {.Deprecated("start"); start(x)})
setGeneric("last", function(x) standardGeneric("last"))
setMethod("last", ".IRanges", function(x) {.Deprecated("end"); end(x)})

