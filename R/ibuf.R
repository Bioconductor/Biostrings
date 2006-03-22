# ===========================================================================
# INTEGER BUFFER objects
# ---------------------------------------------------------------------------
#
# The "ibuf" class implements the concept of "bbuf" objects but
# for integers instead of bytes.
# Some differences between "integer" and "ibuf":
#   1. an "ibuf" object can't be of length 0 (ibuf(0) produces an error)
#   2. ibuf(10) does not initialize its values (integer(10) does)
#   3. ibuf(10)[i] produces an error if i is out of bounds
#   4. "ibuf"objects are fast:
#        > a <- integer(100000000)
#        > system.time(tmp <- a[])
#        [1] 0.65 0.30 0.95 0.00 0.00
#        > system.time(a[] <- 100:1)
#        [1] 3.08 0.52 3.61 0.00 0.00
#
#        > ib <- ibuf(100000000)
#        > system.time(tmp <- ib[])
#        [1] 0.39 0.52 0.91 0.00 0.00
#        > system.time(ib[] <- 100:1)
#        [1] 0.56 0.00 0.56 0.00 0.00

setClass("ibuf", representation(xp="externalptr"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Initialization

setMethod("initialize", "ibuf",
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
        .Call("ibuf_alloc", xp, length, PACKAGE="Biostrings")
        .Object@xp <- xp
        .Object
    }
)

ibuf <- function(...)
{
    new("ibuf", ...)
}

setMethod("show", "ibuf",
    function(object)
    {
        .Call("ibuf_show", object@xp, PACKAGE="Biostrings")
    }
)

setMethod("length", "ibuf",
    function(x)
    {
        .Call("ibuf_length", x@xp, PACKAGE="Biostrings")
    }
)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Read/write functions.
# These are safe wrappers to unsafe C functions.
# If length(i) then read functions return an empty vector and write functions
# don't do anything.

ibReadInts <- function(x, i, imax=integer(0))
{
    if (!is.integer(i))
        i <- as.integer(i)
    if (length(i) == 1) {
        if (length(imax) == 0)
            imax <- i
        else
            imax <- as.integer(imax)
        .Call("ibuf_read_ints", x@xp, i, imax, PACKAGE="Biostrings")
    } else {
        .Call("ibuf_readii_ints", x@xp, i, PACKAGE="Biostrings")
    }
}

ibWriteInts <- function(x, i, imax=integer(0), value)
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
        .Call("ibuf_write_ints", x@xp, i, imax, value, PACKAGE="Biostrings")
    } else {
        .Call("ibuf_writeii_ints", x@xp, i, value, PACKAGE="Biostrings")
    }
    x
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# length(as.integer(bb)) is equivalent to length(bb)
# but the latter is MUCH faster!
setMethod("as.integer", "ibuf",
    function(x)
    {
        ibReadInts(x, 1, length(x))
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Subsetting

setMethod("[", "ibuf",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            return(as.integer(x))
        ibReadInts(x, i)
    }
)

setReplaceMethod("[", "ibuf",
    function(x, i, j,..., value)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")

        # We want to allow this: bb[3] <- 4, even if storage.mode(value)
        # is not "integer"
        if (!is.integer(value)) {
            if (length(value) >= 2)
                stop("'storage.mode(value)' must be \"integer\"")
            tmp <- value
            value <- as.integer(value)
            if (value != tmp)
                stop("'value' is not an integer")
        }
        # Now 'value' is an integer vector
        if (missing(i))
            return(ibWriteInts(x, 1, length(x), value=value))
        ibWriteInts(x, i, value=value)
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Equality

setMethod("==", signature(e1="ibuf", e2="ibuf"),
    function(e1, e2)
    {
        address(e1@xp) == address(e2@xp)
    }
)

# A wrapper to the very fast memcmp() C-function.
# Arguments MUST be the following or it will crash R:
#   x1, x2: "ibuf" objects
#   first1, first2, width: single integers
# In addition: 1 <= first1 <= first1+width-1 <= length(x1)
#              1 <= first2 <= first2+width-1 <= length(x2)
# WARNING: This function is voluntarly unsafe (it doesn't check its
# arguments) because we want it to be the fastest possible!
ibCompare <- function(x1, first1, x2, first2, width)
{
    .Call("ibuf_memcmp", x1@xp, first1, x2@xp, first2, width, PACKAGE="Biostrings")
}
