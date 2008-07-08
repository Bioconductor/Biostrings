### =========================================================================
### IRanges objects
### -------------------------------------------------------------------------
###
### The IRanges class is a simple container for storing a set of integer
### ranges.
###

setClass("IRanges",
    representation(
        start="integer",
        width="integer",
        NAMES="character"  # R doesn't like @names !!
    ),
    prototype(
        start=integer(0),
        width=integer(0),
        NAMES=as.character(NA)
    )
)

setClass("UnlockedIRanges", contains="IRanges")

### A Views object is an UnlockedIRanges object with no null widths.
setClass("Views", contains="UnlockedIRanges")

setClass("LockedIRanges", contains="IRanges")

### A NormalIRanges object is a LockedIRanges object where the ranges are:
###   (a) not empty (i.e. they have a non-null width);
###   (b) not overlapping;
###   (c) ordered from left to right;
###   (d) not even adjacent (i.e. there must be a non empty gap between 2
###       consecutive ranges).
### If 'x' is an IRanges object of length >= 2, then 'x' is normal iff:
###   start(x)[i] <= end(x)[i] < start(x)[i+1] <= end(x)[i+1]
### for every 1 <= i < length(x).
### If length(x) == 1, then 'x' is normal iff width(x)[1] >= 1.
### If length(x) == 0, then 'x' is normal.
setClass("NormalIRanges", contains="LockedIRanges")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("length", "IRanges", function(x) length(x@start))

### The substr() function uses 'start' and 'stop'.
### The substring() function uses 'first' and 'last'.
### We use 'start' and 'end'.
### Note that the "start" and "end" generics are defined in the stats package.
setMethod("start", "IRanges", function(x, ...) x@start)

setGeneric("width", function(x) standardGeneric("width"))
setMethod("width", "IRanges", function(x) x@width)

### Note that when width(x)[i] is 0, then end(x)[i] is start(x)[i] - 1
setMethod("end", "IRanges", function(x, ...) {start(x) + width(x) - 1L})

setMethod("names", "IRanges",
    function(x)
        if (length(x@NAMES) == 1 && is.na(x@NAMES)) NULL else x@NAMES
)

### "desc" is an alias for "names". It might be deprecated soon...
desc <- function(x) names(x)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "isNormal" and "whichFirstNotNormal" generics and methods.
###

setGeneric("isNormal", function(x) standardGeneric("isNormal"))

setMethod("isNormal", "IRanges",
    function(x)
    {
        all_ok <- all(width(x) >= 1)
        if (length(x) >= 2)
            all_ok <- all_ok && all(start(x)[-1] - end(x)[-length(x)] >= 2)
        all_ok
    }
)

setGeneric("whichFirstNotNormal", function(x) standardGeneric("whichFirstNotNormal"))

setMethod("whichFirstNotNormal", "IRanges",
    function(x)
    {
        is_ok <- width(x) >= 1
        if (length(x) >= 2)
            is_ok <- is_ok & c(TRUE, start(x)[-1] - end(x)[-length(x)] >= 2)
        which(!is_ok)[1]
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "isEmpty" generic and methods.
###
### An IRanges object is considered empty iff all its ranges are empty.
### 

setGeneric("isEmpty", function(x) standardGeneric("isEmpty"))

setMethod("isEmpty", "IRanges", function(x) all(width(x) == 0))

### The "isEmpty" method for IRanges objects would work fine on
### NormalIRanges objects but it can be made faster.
setMethod("isEmpty", "NormalIRanges", function(x) length(x) == 0)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "max" and "min" methods.
###
### Note: defined for NormalIRanges objects only.
### For an ordinary IRanges object 'x', it's not clear what the semantic
### should. In particular, should empty ranges be ignored or not? If not then
### we could end up with 'min(x)' > 'max(x)' (e.g. when 'x' is made of 1 empty
### range) which is not nice. Another (and more pragmatic) reason for not
### defining these methods for IRanges objects is that I don't need them at
### the moment.
###

setMethod("max", "NormalIRanges",
    function(x, ..., na.rm)
    {
        if (isEmpty(x)) {
            warning("empty ", class(x), " object; returning -Inf")
            -Inf
        } else {
            end(x)[length(x)]
        }
    }
)

setMethod("min", "NormalIRanges",
    function(x, ..., na.rm)
    {
        if (isEmpty(x)) {
            warning("empty ", class(x), " object; returning Inf")
            Inf
        } else {
            start(x)[1]
        }
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

### IRanges objects
.valid.IRanges.start <- function(object)
{
    if (!is.integer(start(object)) || any(is.na(start(object))))
        return("the starts must be non-NA integers")
    if (length(start(object)) != length(width(object)))
        return("number of starts and number of widths differ")
    NULL
}
.valid.IRanges.width <- function(object)
{
    if (!is.integer(width(object)) || any(is.na(width(object))))
        return("the widths must be non-NA integers")
    if (length(start(object)) != length(width(object)))
        return("number of starts and number of widths differ")
    if (length(width(object)) != 0 && min(width(object)) < 0L)
        return("negative widths are not allowed")
    NULL
}
.valid.IRanges.names <- function(object)
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
.valid.IRanges <- function(object)
{
    #cat("validating IRanges object of length", length(object), "...\n")
    c(.valid.IRanges.start(object),
      .valid.IRanges.width(object),
      .valid.IRanges.names(object))
}
setValidity("IRanges",
    function(object)
    {
        problems <- .valid.IRanges(object)
        if (is.null(problems)) TRUE else problems
    }
)

### Views objects
.valid.Views.width <- function(object)
{
    if (length(object) != 0 && any(width(object) == 0))
        return("null widths are not allowed")
    NULL
}
setValidity("Views",
    function(object)
    {
        problems <- .valid.Views.width(object)
        if (is.null(problems)) TRUE else problems
    }
)

### NormalIRanges objects
.valid.NormalIRanges <- function(object)
{
    if (!isNormal(object))
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
### Initialization.
###
### I found that using validObject() in "initialize" doesn't work properly
### (validation is called too many times and not in an order that makes sense
### to me...)
###

.set.IRanges.slots <- function(object, start, width, names, check=TRUE)
{
    slot(object, "start", check=FALSE) <- numeric2integer(start)
    slot(object, "width", check=FALSE) <- numeric2integer(width)
    slot(object, "NAMES", check=FALSE) <- if (is.null(names)) as.character(NA) else names
    if (check)
        stopIfProblems(.valid.IRanges(object))
    object
}

setMethod("initialize", "IRanges",
    function(.Object, start=integer(0), width=integer(0),
                      names=NULL, check=TRUE)
        .set.IRanges.slots(.Object, start, width, names, check=check)
)

setMethod("initialize", "Views",
    function(.Object, start=integer(0), width=integer(0),
                      names=NULL, check=TRUE)
    {
        .Object <- callNextMethod(.Object, start=start, width=width,
                                           names=names, check=check)
        if (check)
            stopIfProblems(.valid.Views.width(.Object))
        .Object
    }
)

setMethod("initialize", "NormalIRanges",
    function(.Object, start=integer(0), width=integer(0),
                      names=NULL, check=TRUE)
    {
        .Object <- callNextMethod(.Object, start=start, width=width,
                                           names=names, check=check)
        if (check)
            stopIfProblems(.valid.NormalIRanges(.Object))
        .Object
    }
)

### For internal use only. No need to export.
newEmptyNormalIRanges <- function() new("NormalIRanges", check=FALSE)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The safe and user-friendly "IRanges" constructor.
###

.IRanges.normargStartEndWidth <- function(start_end_width, argname)
{
    if (is.null(start_end_width))
        return(start_end_width)
    if (!is.numeric(start_end_width))
        stop("'", argname, "' must be a numeric vector (or NULL)")
#    if (length(start_end_width) == 0)
#        stop("'", argname, "' must contain at least one integer value")
    if (!is.integer(start_end_width))
        start_end_width <- as.integer(start_end_width)
    if (any(is.na(start_end_width)))
        stop("'", argname,, "' cannot contain NAs")
    start_end_width
}

IRanges <- function(start=NULL, end=NULL, width=NULL)
{
    start <- .IRanges.normargStartEndWidth(start, "start")
    end <- .IRanges.normargStartEndWidth(end, "end")
    width <- .IRanges.normargStartEndWidth(width, "width")
    if (sum(!c(is.null(start), is.null(end), is.null(width))) != 2)
        stop("exactly two of the 'start', 'end' and 'width' arguments must be specified")
    if (is.null(width)) { 
        if (length(start) != length(end))
            stop("'start' and 'end' must have the same length")
        if (length(start) != 0)
            width <- end - start + 1L
        else
            width <- integer(0)
    } else if (is.null(start)) {
        if (length(width) > length(end))
            stop("'width' has more elements than 'end'")
        if (length(width) != 0)
            start <- end - width + 1L
        else if (length(end) == 0)
            start <- integer(0)
        else
            stop("cannot recycle zero length 'width'")
    } else {
        if (length(width) > length(start))
            stop("'width' has more elements than 'start'")
        if (length(width) < length(start)) {
            if (length(width) != 0)
                width <- rep(width, length.out=length(start))
            else
                stop("cannot recycle zero length 'width'")
        }
    }
    new("IRanges", start=start, width=width)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###
### We cannot rely on the implicit "coerce" methods for coercing an arbitrary
### IRanges object into a NormalIRanges object because they do NOT check that
### the returned object is valid! Yes, implicit "coerce" methods were supposed
### to be a nice S4 "feature"...
###

### NOT exported and unsafe: 'from' MUST be an IRanges object.
asNormalIRanges <- function(from, check=TRUE)
{
    new("NormalIRanges", start=start(from), width=width(from),
                         names=names(from), check=check)
}

.asNormalIRanges <- function(from) asNormalIRanges(from, check=TRUE)

### No, defining the IRanges->NormalIRanges "coerce" method is not enough and
### we also need to define the other methods! Otherwise a silly implicit
### method would be called when calling as(x, "NormalIRanges") on an
### UnlockedIRanges, Views or LockedIRanges object. Yes, this is another S4
### "feature":
###   https://stat.ethz.ch/pipermail/r-devel/2008-April/049027.html
setAs("IRanges", "NormalIRanges", .asNormalIRanges)
setAs("UnlockedIRanges", "NormalIRanges", .asNormalIRanges)
setAs("Views", "NormalIRanges", .asNormalIRanges)
setAs("LockedIRanges", "NormalIRanges", .asNormalIRanges)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

setMethod("as.data.frame", "IRanges",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        if (!(is.null(row.names) || is.character(row.names)))
            stop("'row.names'  must be NULL or a character vector")
        ans <- data.frame(start=start(x),
                          end=end(x),
                          width=width(x),
                          row.names=row.names,
                          check.rows=TRUE,
                          check.names=FALSE,
                          stringsAsFactors=FALSE)
        ans$names <- names(x)
        ans
    }
)

setMethod("show", "IRanges",
    function(object)
    {
        cat(class(object), " object:\n", sep="")
        show(as.data.frame(object))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Non-exported (and unsafe) replacement functions for IRanges objects.
###
### IMPORTANT: They do NOT check their arguments ('x' and 'value'). In
### particular they do not check that 'value' is of the expected type (integer
### for "unsafe.start<-", "unsafe.width<-", "unsafe.end<-", and character for
### "unsafe.names<-").
###
### The "sliding rules" are:
###   (1) changing the start preserves the width (so it changes the end)
###   (2) changing the width preserves the start (so it changes the end)
###   (3) changing the end preserves the width (so it changes the start)
###

### 'value' is recycled.
`unsafe.start<-` <- function(x, value)
{
    ## Use 'x@start[]' instead of just 'x@start' so the right value is recycled
    x@start[] <- numeric2integer(value)
    x
}

### 'value' is recycled.
`unsafe.width<-` <- function(x, value)
{
    ## Use 'x@width[]' instead of just 'x@width' so the right value is recycled
    x@width[] <- numeric2integer(value)
    x
}

### 'value' is recycled.
`unsafe.end<-` <- function(x, value)
{
    unsafe.start(x) <- numeric2integer(value) - width(x) + 1L
    x
}

### 'value' is NOT recycled so we stay close to what the standard R "names<-"
### methods generally do (except that if 'value' is shorter than 'x', we extend
### it by empty strings, not by character NAs).
`unsafe.names<-` <- function(x, value)
{
    if (is.null(value))
        x@NAMES <- as.character(NA)
    else {
        if (length(value) > length(x))
            stop("too many names")
        if (length(value) < length(x))
            value <- c(value, character(length(x) - length(value)))
        x@NAMES <- value
    }
    x
}

unsafe.update <- function(object, ...)
{
    valid_argnames <- c("start", "end", "width", "names")
    args <- extraArgsAsList(valid_argnames, ...)
    argnames <- names(args)
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
        return(.set.IRanges.slots(object, start, width, args$names, check=FALSE))
    }
    if ("start" %in% argnames)
        unsafe.start(object) <- args$start
    if ("end" %in% argnames)
        unsafe.end(object) <- args$end
    if ("width" %in% argnames)
        unsafe.width(object) <- args$width
    if ("names" %in% argnames)
        unsafe.names(object) <- args$names
    object
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Exported (and safe) replacement methods.
###
### See the unsafe replacement functions for IRanges objects above for the
### "sliding rules".
###
### Note that we don't call validObject(x) after 'x' has been modified because
### we don't need to revalidate the entire object: validating the bits that
### have been touched is enough (and faster). However, because of this, when
### defining a new class that extends the UnlockedIRanges class, if objects of
### the new class must satisfy additional constraints, then some of the
### replacement methods below need to be overridden for this new class.
###

setGeneric("start<-", signature="x",
    function(x, check=TRUE, value) standardGeneric("start<-")
)

setReplaceMethod("start", "UnlockedIRanges",
    function(x, check=TRUE, value)
    {
        unsafe.start(x) <- value
        if (check)
            stopIfProblems(.valid.IRanges.start(x))
        x
    }
)

setGeneric("width<-", signature="x",
    function(x, check=TRUE, value) standardGeneric("width<-")
)

setReplaceMethod("width", "UnlockedIRanges",
    function(x, check=TRUE, value)
    {
        unsafe.width(x) <- value
        if (check) {
            stopIfProblems(.valid.IRanges.width(x))
            if (is(x, "Views"))
                stopIfProblems(.valid.Views.width(x))
        }
        x
    }
)

setGeneric("end<-", signature="x",
    function(x, check=TRUE, value) standardGeneric("end<-")
)

setReplaceMethod("end", "UnlockedIRanges",
    function(x, check=TRUE, value)
    {
        unsafe.end(x) <- value
        if (check)
            stopIfProblems(.valid.IRanges.start(x))
        x
    }
)

### Yes, for IRanges objects!
setReplaceMethod("names", "IRanges",
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


### "desc<-" is an alias for "names<-". It might be deprecated soon...
`desc<-` <- function(x, value) `names<-`(x, value)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "update" method.
###
### This is a convenience method for combining multiple modifications in one
### single call.
###
### It must verify 2 important properties:
###   (1) update(x) must be identical to x (doesn't touch x at all)
###   (2) update(x, start=start(x), width=width(x), names=names(x))
###       must be identical to x too (but this time it updates x with its own
###       content)
###

setMethod("update", "UnlockedIRanges",
    function(object, ...)
    {
        object <- unsafe.update(object, ...)
        check <- list(...)$check
        if (is.null(check))
            check <- TRUE
        if (check)
            validObject(object)
        object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Supported 'i' types: numeric vector, logical vector, NULL and missing.
setMethod("[", "IRanges",
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
            if (is(x, "NormalIRanges") && all(i >= 0)) {
                i <- i[i != 0]
                if (isNotStrictlySorted(i))
                    stop("positive numeric subscript must be strictly increasing ",
                         "for NormalIRanges objects")
            }
        } else if (is.logical(i)) {
            if (any(is.na(i)))
                stop("subscript contains NAs")
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

setReplaceMethod("[", "IRanges",
    function(x, i, j,..., value)
        stop("attempt to modify the value of a ", class(x), " instance")
)

setMethod("rep", "IRanges",
    function(x, times)
        x[rep.int(seq_len(length(x)), times)]
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "duplicated" method.
###
### TODO: current implementation is very inefficient and needs some C help!
###

setMethod("duplicated", "IRanges",
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

setMethod("as.matrix", "IRanges",
    function(x, ...)
        matrix(data=c(start(x), width(x)), ncol=2, dimnames=list(names(x), NULL))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Deprecated generics and methods.
###

setGeneric("first", function(x) standardGeneric("first"))
setMethod("first", "IRanges", function(x) {.Deprecated("start"); start(x)})
setGeneric("last", function(x) standardGeneric("last"))
setMethod("last", "IRanges", function(x) {.Deprecated("end"); end(x)})

