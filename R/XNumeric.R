### =========================================================================
### External numeric vectors: the "XNumeric" class
### -------------------------------------------------------------------------
###
### The "XNumeric" class implements the concept of "XRaw" objects but
### for numerics instead of bytes.
### Some differences between "numeric" and "XNumeric":
###
###   1. XNumeric(10) does not initialize its values (numeric(10) does)
###
###   2. XNumeric(10)[i] produces an error if i is out of bounds
###
###   3. XNumeric objects are faster:
###
###        > a <- numeric(100000000)
###        > system.time(tmp <- a[])
###        [1] 0.65 0.30 0.95 0.00 0.00
###        > system.time(a[] <- 100:1)
###        [1] 3.08 0.52 3.61 0.00 0.00
###
###        > b <- XNumeric(100000000)
###        > system.time(tmp <- b[])
###        [1] 0.39 0.52 0.91 0.00 0.00
###        > system.time(b[] <- 100:1)
###        [1] 0.56 0.00 0.56 0.00 0.00

setClass("XNumeric", representation(xp="externalptr"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
###
### Note that unlike numeric vectors, XNumeric objects are not initialized
### with 0's.
###

### This:
###   xn <- XNumeric(30)
### will call this "initialize" method.
setMethod("initialize", "XNumeric",
    function(.Object, length, verbose=FALSE)
    {
        if (!isSingleNumber(length))
            stop("'length' must be a single integer")
        length <- as.integer(length)
        if (length < 0)
            stop("'length' must be a non-negative integer")
        xp <- .Call("Biostrings_xp_new", PACKAGE="Biostrings")
        if (verbose)
            cat("Allocating memory for new", class(.Object), "object...")
        .Call("XNumeric_alloc", xp, length, PACKAGE="Biostrings")
        if (verbose) {
            cat(" OK\n")
            show_string <- .Call("XNumeric_get_show_string", xp, PACKAGE="Biostrings")
            cat("New", show_string, "successfully created\n")
        }
        .Object@xp <- xp
        .Object
    }
)

XNumeric <- function(...)
{
    new("XNumeric", ...)
}

setMethod("show", "XNumeric",
    function(object)
    {
        show_string <- .Call("XNumeric_get_show_string", object@xp, PACKAGE="Biostrings")
        cat(show_string, "\n", sep="")
        ## What is correct here? The documentation (?show) says that 'show'
        ## should return an invisible 'NULL' but, on the other hand, the 'show'
        ## method for numerics returns its 'object' argument...
        invisible(object)
    }
)

setMethod("length", "XNumeric",
    function(x)
    {
        .Call("XNumeric_length", x@xp, PACKAGE="Biostrings")
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Read/write functions.
### These are almost safe wrappers to unsafe C functions ("almost" because
### they don't check for NAs in their arguments).
### If length(i) == 0 then the read functions return an empty vector
### and the write functions don't do anything.

XNumeric.read <- function(x, i, imax=integer(0))
{
    if (!is.integer(i))
        i <- as.integer(i)
    if (length(i) == 1) {
        if (length(imax) == 0)
            imax <- i
        else
            imax <- as.integer(imax)
        .Call("XNumeric_read_nums_from_i1i2", x@xp, i, imax, PACKAGE="Biostrings")
    } else {
        .Call("XNumeric_read_nums_from_subset", x@xp, i, PACKAGE="Biostrings")
    }
}

XNumeric.write <- function(x, i, imax=integer(0), value)
{
    if (!is.numeric(value))
        stop("'value' is not a numeric vector")
    if (!is.integer(i))
        i <- as.integer(i)
    if (length(i) == 1) {
        if (length(imax) == 0)
            imax <- i
        else
            imax <- as.integer(imax)
        .Call("XNumeric_write_nums_to_i1i2", x@xp, i, imax, value, PACKAGE="Biostrings")
    } else {
        .Call("XNumeric_write_nums_to_subset", x@xp, i, value, PACKAGE="Biostrings")
    }
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### length(as.numeric(b)) is equivalent to length(b)
### but the latter is MUCH faster!
setMethod("as.numeric", "XNumeric",
    function(x)
    {
        XNumeric.read(x, 1, length(x))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting

setMethod("[", "XNumeric",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            return(as.numeric(x))
        XNumeric.read(x, i)
    }
)

setReplaceMethod("[", "XNumeric",
    function(x, i, j,..., value)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")

        ## We want to allow this: b[3] <- 4, even if storage.mode(value)
        ## is not "numeric"
        if (!is.numeric(value)) {
            if (length(value) >= 2)
                stop("'storage.mode(value)' must be \"numeric\"")
            tmp <- value
            value <- as.numeric(value)
            if (value != tmp)
                stop("'value' is not numeric")
        }
        ## Now 'value' is a numeric vector
        if (missing(i))
            return(XNumeric.write(x, 1, length(x), value=value))
        XNumeric.write(x, i, value=value)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Equality

setMethod("==", signature(e1="XNumeric", e2="XNumeric"),
    function(e1, e2)
    {
        address(e1@xp) == address(e2@xp)
    }
)
setMethod("!=", signature(e1="XNumeric", e2="XNumeric"),
    function(e1, e2)
    {
        address(e1@xp) != address(e2@xp)
    }
)

### A wrapper to the very fast memcmp() C-function.
### Arguments MUST be the following or it will crash R:
###   x1, x2: "XNumeric" objects
###   start1, start2, width: single integers
### In addition: 1 <= start1 <= start1+width-1 <= length(x1)
###              1 <= start2 <= start2+width-1 <= length(x2)
### WARNING: This function is voluntarly unsafe (it doesn't check its
### arguments) because we want it to be the fastest possible!
XNumeric.compare <- function(x1, start1, x2, start2, width)
{
    .Call("XNumeric_memcmp", x1@xp, start1, x2@xp, start2, width, PACKAGE="Biostrings")
}
