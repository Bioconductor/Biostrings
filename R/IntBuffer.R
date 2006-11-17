### =========================================================================
### INTEGER BUFFER objects
### -------------------------------------------------------------------------
###
### The "IntBuffer" class implements the concept of "CharBuffer" objects but
### for integers instead of chars.
### Some differences between "integer" and "IntBuffer":
###   1. an "IntBuffer" object can't be of length 0 (IntBuffer(0) produces an error)
###   2. IntBuffer(10) does not initialize its values (integer(10) does)
###   3. IntBuffer(10)[i] produces an error if i is out of bounds
###   4. "IntBuffer" objects are fast:
###        > a <- integer(100000000)
###        > system.time(tmp <- a[])
###        [1] 0.65 0.30 0.95 0.00 0.00
###        > system.time(a[] <- 100:1)
###        [1] 3.08 0.52 3.61 0.00 0.00
###
###        > ib <- IntBuffer(100000000)
###        > system.time(tmp <- ib[])
###        [1] 0.39 0.52 0.91 0.00 0.00
###        > system.time(ib[] <- 100:1)
###        [1] 0.56 0.00 0.56 0.00 0.00

setClass("IntBuffer", representation(xp="externalptr"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization

setMethod("initialize", "IntBuffer",
    function(.Object, length)
    {
        if (missing(length))
            stop("argument 'length' is missing")
        if (!is.numeric(length) || is.nan(length))
            stop("'length' is not a number")
        if (length(length) != 1)
            stop("'length' must have only one element")
        length <- as.integer(length)
        if (length < 1)
            stop("buffer length must be >= 1")
        xp <- .Call("xp_new", PACKAGE="Biostrings")
        .Call("IntBuffer_alloc", xp, length, PACKAGE="Biostrings")
        .Object@xp <- xp
        .Object
    }
)

IntBuffer <- function(...)
{
    new("IntBuffer", ...)
}

setMethod("show", "IntBuffer",
    function(object)
    {
        .Call("IntBuffer_show", object@xp, PACKAGE="Biostrings")
    }
)

setMethod("length", "IntBuffer",
    function(x)
    {
        .Call("IntBuffer_length", x@xp, PACKAGE="Biostrings")
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Read/write functions.
### These are safe wrappers to unsafe C functions.
### If length(i) == 0 then the read functions return an empty vector
### and the write functions don't do anything.

IntBuffer.read <- function(x, i, imax=integer(0))
{
    if (!is.integer(i))
        i <- as.integer(i)
    if (length(i) == 1) {
        if (length(imax) == 0)
            imax <- i
        else
            imax <- as.integer(imax)
        .Call("IntBuffer_read_ints_from_i1i2", x@xp, i, imax, PACKAGE="Biostrings")
    } else {
        .Call("IntBuffer_read_ints_from_subset", x@xp, i, PACKAGE="Biostrings")
    }
}

IntBuffer.write <- function(x, i, imax=integer(0), value)
{
    if (!is.integer(value))
        stop("'value' is not an integer vector")
    if (!is.integer(i))
        i <- as.integer(i)
    if (length(i) == 1) {
        if (length(imax) == 0)
            imax <- i
        else
            imax <- as.integer(imax)
        .Call("IntBuffer_write_ints_to_i1i2", x@xp, i, imax, value, PACKAGE="Biostrings")
    } else {
        .Call("IntBuffer_write_ints_to_subset", x@xp, i, value, PACKAGE="Biostrings")
    }
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### length(as.integer(ib)) is equivalent to length(ib)
### but the latter is MUCH faster!
setMethod("as.integer", "IntBuffer",
    function(x)
    {
        IntBuffer.read(x, 1, length(x))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting

setMethod("[", "IntBuffer",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            return(as.integer(x))
        IntBuffer.read(x, i)
    }
)

setReplaceMethod("[", "IntBuffer",
    function(x, i, j,..., value)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")

        ## We want to allow this: ib[3] <- 4, even if storage.mode(value)
        ## is not "integer"
        if (!is.integer(value)) {
            if (length(value) >= 2)
                stop("'storage.mode(value)' must be \"integer\"")
            tmp <- value
            value <- as.integer(value)
            if (value != tmp)
                stop("'value' is not an integer")
        }
        ## Now 'value' is an integer vector
        if (missing(i))
            return(IntBuffer.write(x, 1, length(x), value=value))
        IntBuffer.write(x, i, value=value)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Equality

setMethod("==", signature(e1="IntBuffer", e2="IntBuffer"),
    function(e1, e2)
    {
        address(e1@xp) == address(e2@xp)
    }
)
setMethod("!=", signature(e1="IntBuffer", e2="IntBuffer"),
    function(e1, e2)
    {
        address(e1@xp) != address(e2@xp)
    }
)

### A wrapper to the very fast memcmp() C-function.
### Arguments MUST be the following or it will crash R:
###   x1, x2: "IntBuffer" objects
###   start1, start2, width: single integers
### In addition: 1 <= start1 <= start1+width-1 <= length(x1)
###              1 <= start2 <= start2+width-1 <= length(x2)
### WARNING: This function is voluntarly unsafe (it doesn't check its
### arguments) because we want it to be the fastest possible!
IntBuffer.compare <- function(x1, start1, x2, start2, width)
{
    .Call("IntBuffer_memcmp", x1@xp, start1, x2@xp, start2, width, PACKAGE="Biostrings")
}
