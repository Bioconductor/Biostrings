# ===========================================================================
# BYTE BUFFER objects
# ---------------------------------------------------------------------------
#
# The "bbuf" class describes "byte buffer" objects.
# A "byte buffer" object is a chunk of memory that:
#   a. Contains bytes (characters).
#   b. Is readable and writable.
#   c. Is not copied on object duplication i.e. when doing
#        bb2 <- bb1
#      both bb1 and bb2 point to the same place in memory.
#      This is achieved by using R predefined type "externalptr".
#   d. Is not 0-terminated (so it can contain zeros). This is achieved by
#      having the length of the buffer stored in the "bbuf" object.


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Return the hexadecimal address of any R object in a string.
address <- function(x)
{
    .Call("sexp_address", x, PACKAGE="Biostrings")
}

# Helper function (for debugging purpose).
# Print some obscure info about an "externalptr" object.
# Typical use:
#   show(new("externalptr"))
setMethod("show", "externalptr",
    function(object)
    {
        .Call("xp_show", object, PACKAGE="Biostrings")
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The "byte buffer" class.
# Note: instead of defining the "bbuf" class with just one slot of type
# "externalptr" (it HAS an "externalptr", and nothing else), an alternative
# would be to simply extend the "externalptr" type.
# After all, a "bbuf" object IS an "externalptr" object.
# However, I tried this but was not able to implement the "initialize" method
# in such a way that it returns a new instance of the "bbuf" class (the
# returned object was ALWAYS the same instance everytime the method was
# called, I found no workaround).
setClass("bbuf", representation(xp="externalptr"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Initialization

# This:
#   bb <- bbuf(30)
# will call the "initialize" method.
setMethod("initialize", "bbuf",
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
        .Call("bbuf_alloc", xp, length, PACKAGE="Biostrings")
        .Object@xp <- xp
        .Object
    }
)

bbuf <- function(...)
{
    new("bbuf", ...)
}

setMethod("show", "bbuf",
    function(object)
    {
        .Call("bbuf_show", object@xp, PACKAGE="Biostrings")
    }
)

setMethod("length", "bbuf",
    function(x)
    {
        .Call("bbuf_length", x@xp, PACKAGE="Biostrings")
    }
)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Read/write functions.
# These are safe wrappers to unsafe C functions.
# If length(i) == 0 then read functions return an empty vector and
# write functions don't do anything.

bbReadInts <- function(x, i, imax=integer(0))
{
    if (!is.integer(i))
        i <- as.integer(i)
    if (length(i) == 1) {
        if (length(imax) == 0)
            imax <- i
        else
            imax <- as.integer(imax)
        .Call("bbuf_read_ints", x@xp, i, imax, PACKAGE="Biostrings")
    } else {
        .Call("bbuf_readii_ints", x@xp, i, PACKAGE="Biostrings")
    }
}

bbWriteInts <- function(x, i, imax=integer(0), value)
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
        .Call("bbuf_write_ints", x@xp, i, imax, value, PACKAGE="Biostrings")
    } else {
        .Call("bbuf_writeii_ints", x@xp, i, value, PACKAGE="Biostrings")
    }
    x
}

bbReadChars <- function(x, i, imax=integer(0), dec_hash=NULL)
{
    if (!is.integer(i))
        i <- as.integer(i)
    if (length(i) == 1) {
        if (length(imax) == 0)
            imax <- i
        else
            imax <- as.integer(imax)
        if (is.null(dec_hash))
            .Call("bbuf_read_chars",
                  x@xp, i, imax, PACKAGE="Biostrings")
        else
            .Call("bbuf_read_enc_chars",
                  x@xp, i, imax, dec_hash@xp, PACKAGE="Biostrings")
    } else {
        if (is.null(dec_hash))
            .Call("bbuf_readii_chars",
                  x@xp, i, PACKAGE="Biostrings")
        else
            .Call("bbuf_readii_enc_chars",
                  x@xp, i, dec_hash@xp, PACKAGE="Biostrings")
    }
}

bbWriteChars <- function(x, i, imax=integer(0), value, enc_hash=NULL)
{
    if (!is.character(value))
        stop("'value' is not a character vector")
    if (!is.integer(i))
        i <- as.integer(i)
    if (length(i) == 1) {
        if (length(imax) == 0)
            imax <- i
        else
            imax <- as.integer(imax)
        if (is.null(enc_hash))
            .Call("bbuf_write_chars",
                  x@xp, i, imax, value, PACKAGE="Biostrings")
        else
            .Call("bbuf_write_enc_chars",
                  x@xp, i, imax, value, enc_hash@xp, PACKAGE="Biostrings")
    } else {
        if (is.null(enc_hash))
            .Call("bbuf_writeii_chars",
                  x@xp, i, value, PACKAGE="Biostrings")
        else
            .Call("bbuf_writeii_enc_chars",
                  x@xp, i, value, enc_hash@xp, PACKAGE="Biostrings")
    }
    x
}

bbCopy <- function(dest, i, imax=integer(0), src)
{
    if (class(src) != "bbuf")
        stop("'src' is not a byte buffer")
    if (!is.integer(i))
        i <- as.integer(i)
    if (length(i) == 1) {
        if (length(imax) == 0)
            imax <- i
        else
            imax <- as.integer(imax)
        .Call("bbuf_copy", dest@xp, i, imax, src@xp, PACKAGE="Biostrings")
    } else {
        .Call("bbuf_copyii", dest@xp, i, src@xp, PACKAGE="Biostrings")
    }
    dest
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# length(as.integer(bb)) is equivalent to length(bb)
# but the latter is MUCH faster!
setMethod("as.integer", "bbuf",
    function(x)
    {
        bbReadInts(x, 1, length(x))
    }
)

# Typical use:
#   bb <- bbuf(15)
#   bb[] <- 65
#   toString(bb)
#   bb[] <- "Hello"
#   toString(bb)
# So this should always rewrite the content of a "bbuf" object
# to itself, without any modification:
#   bb[] <- toString(bb)
# whetever the content of bb is!
setMethod("toString", "bbuf",
    function(x)
    {
        bbReadChars(x, 1, length(x))
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Subsetting

# Select bytes from a "byte buffer".
# Typical use:
#   bb <- bbuf(30)
#   bb[25:20]
#   bb[25:31] # subscript out of bounds
# Note: bb[] can be used as a shortcut for as.integer(bb)
setMethod("[", "bbuf",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            return(as.integer(x))
        bbReadInts(x, i)
    }
)

# Replace bytes in a "byte buffer".
# Typical use:
#   bb <- bbuf(30)
#   bb[] <- 12 # fill with 12
#   bb[3:10] <- 1:-2
#   bb[3:10] <- "Ab"
#   bb[0] <- 4 # subscript out of bounds
#   bb[31] <- 4 # subscript out of bounds
#   bb[3] <- -12 # subscript out of bounds
setReplaceMethod("[", "bbuf",
    function(x, i, j,..., value)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")

        # 'value' is a string
        if (is.character(value)) {
            if (length(value) >= 2)
                stop("character vector 'value' has more than one string")
            if (missing(i))
                return(bbWriteChars(x, 1, length(x), value=value))
            return(bbWriteChars(x, i, value=value))
        }

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
            return(bbWriteInts(x, 1, length(x), value=value))
        bbWriteInts(x, i, value=value)
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Equality

# Be careful to the semantic of the "==" operator:
#   2 "bbuf" objects are equals if their @xp slot is the
#   same "externalptr" instance (then they obviously have
#   the same length and contain the same data).
# With this definition, bb1 and bb2 can be 2 different "bbuf" objects
# (bb1 != bb2) and contain the same data.
setMethod("==", signature(e1="bbuf", e2="bbuf"),
    function(e1, e2)
    {
        address(e1@xp) == address(e2@xp)
    }
)

# A wrapper to the very fast memcmp() C-function.
# Arguments MUST be the following or it will crash R:
#   x1, x2: "bbuf" objects
#   first1, first2, width: single integers
# In addition: 1 <= first1 <= first1+width-1 <= length(x1)
#              1 <= first2 <= first2+width-1 <= length(x2)
# WARNING: This function is voluntarly unsafe (it doesn't check its
# arguments) because we want it to be the fastest possible!
bbCompare <- function(x1, first1, x2, first2, width)
{
    .Call("bbuf_memcmp", x1@xp, first1, x2@xp, first2, width, PACKAGE="Biostrings")
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The "cbuf" class is a simple extention of the "bbuf" class (no additional
# slots)

setClass("cbuf", representation("bbuf"))

setMethod("show", "cbuf",
    function(object)
    {
        print(toString(object))
    }
)

# Safe alternative to 'strsplit(x, NULL)'.
safeExplode <- function(x)
{
    if (!is.character(x) || length(x) != 1)
        stop("'x' must be a single string")
    .Call("safe_explode", x, PACKAGE="Biostrings")
}

setMethod("[", "cbuf",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            s <- toString(x)
        else
            s <- bbReadChars(x, i)
        safeExplode(s)
    }
)
