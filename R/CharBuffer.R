### =========================================================================
### "CharBuffer" objects
### -------------------------------------------------------------------------
###
### A "CharBuffer" object is a chunk of memory that:
###   a. Contains bytes (characters).
###   b. Is readable and writable.
###   c. Is not copied on object duplication i.e. when doing
###        cb2 <- cb1
###      both cb1 and cb2 point to the same place in memory.
###      This is achieved by using R predefined type "externalptr".
###   d. Is not 0-terminated (so it can contain zeros). This is achieved by
###      having the length of the buffer stored in the "CharBuffer" object.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

debug_utils <- function()
{
    invisible(.Call("utils_debug", PACKAGE="Biostrings"))
}

debug_CharBuffer <- function()
{
    invisible(.Call("CharBuffer_debug", PACKAGE="Biostrings"))
}


### Return the hexadecimal address of any R object in a string.
address <- function(x)
{
    .Call("sexp_address", x, PACKAGE="Biostrings")
}

### Helper function (for debugging purpose).
### Print some obscure info about an "externalptr" object.
### Typical use:
###   show(new("externalptr"))
setMethod("show", "externalptr",
    function(object)
    {
        .Call("xp_show", object, PACKAGE="Biostrings")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "CharBuffer" class.
### Note: instead of defining the "CharBuffer" class with just one slot of type
### "externalptr" (it HAS an "externalptr", and nothing else), an alternative
### would be to simply extend the "externalptr" type.
### After all, a "CharBuffer" object IS an "externalptr" object.
### However, I tried this but was not able to implement the "initialize" method
### in such a way that it returns a new instance of the "CharBuffer" class (the
### returned object was ALWAYS the same instance everytime the method was
### called, I found no workaround).
setClass("CharBuffer", representation(xp="externalptr"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization

### This:
###   cb <- CharBuffer(30)
### will call the "initialize" method.
setMethod("initialize", "CharBuffer",
    function(.Object, length, verbose=FALSE)
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
        .Call("CharBuffer_alloc", xp, length, PACKAGE="Biostrings")
        .Object@xp <- xp
        if (verbose) {
            show_string <- .Call("CharBuffer_get_show_string", .Object@xp, PACKAGE="Biostrings")
            cat("Allocating new ", show_string, "\n", sep="")
        }
        .Object
    }
)

CharBuffer <- function(...)
{
    new("CharBuffer", ...)
}

setMethod("show", "CharBuffer",
    function(object)
    {
        show_string <- .Call("CharBuffer_get_show_string", object@xp, PACKAGE="Biostrings")
        cat(show_string, "\n", sep="")
        ## What is correct here? The documentation (?show) says that 'show'
        ## should return an invisible 'NULL' but, on the other hand, the 'show'
        ## method for intergers returns its 'object' argument...
        invisible(object)
    }
)

setMethod("length", "CharBuffer",
    function(x)
    {
        .Call("CharBuffer_length", x@xp, PACKAGE="Biostrings")
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Read/write functions.
### These are safe wrappers to unsafe C functions.
### If length(i) == 0 then the read functions return an empty vector
### and the write functions don't do anything.

CharBuffer.readInts <- function(x, i, imax=integer(0))
{
    if (!is.integer(i))
        i <- as.integer(i)
    if (length(i) == 1) {
        if (length(imax) == 0)
            imax <- i
        else
            imax <- as.integer(imax)
        .Call("CharBuffer_read_ints_from_i1i2", x@xp, i, imax, PACKAGE="Biostrings")
    } else {
        .Call("CharBuffer_read_ints_from_subset", x@xp, i, PACKAGE="Biostrings")
    }
}

CharBuffer.writeInts <- function(x, i, imax=integer(0), value)
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
        .Call("CharBuffer_write_ints_to_i1i2", x@xp, i, imax, value, PACKAGE="Biostrings")
    } else {
        .Call("CharBuffer_write_ints_to_subset", x@xp, i, value, PACKAGE="Biostrings")
    }
    x
}

CharBuffer.read <- function(x, i, imax=integer(0), dec_hash=NULL)
{
    if (!is.integer(i))
        i <- as.integer(i)
    if (length(i) == 1) {
        if (length(imax) == 0)
            imax <- i
        else
            imax <- as.integer(imax)
        if (is.null(dec_hash))
            .Call("CharBuffer_read_chars_from_i1i2",
                  x@xp, i, imax, PACKAGE="Biostrings")
        else
            .Call("CharBuffer_read_enc_chars_from_i1i2",
                  x@xp, i, imax, dec_hash@xp, PACKAGE="Biostrings")
    } else {
        if (is.null(dec_hash))
            .Call("CharBuffer_read_chars_from_subset",
                  x@xp, i, PACKAGE="Biostrings")
        else
            .Call("CharBuffer_read_enc_chars_from_subset",
                  x@xp, i, dec_hash@xp, PACKAGE="Biostrings")
    }
}

CharBuffer.write <- function(x, i, imax=integer(0), value, enc_hash=NULL)
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
            .Call("CharBuffer_write_chars_to_i1i2",
                  x@xp, i, imax, value, PACKAGE="Biostrings")
        else
            .Call("CharBuffer_write_enc_chars_to_i1i2",
                  x@xp, i, imax, value, enc_hash@xp, PACKAGE="Biostrings")
    } else {
        if (is.null(enc_hash))
            .Call("CharBuffer_write_chars_to_subset",
                  x@xp, i, value, PACKAGE="Biostrings")
        else
            .Call("CharBuffer_write_enc_chars_to_subset",
                  x@xp, i, value, enc_hash@xp, PACKAGE="Biostrings")
    }
    x
}

CharBuffer.copy <- function(dest, i, imax=integer(0), src, hash=NULL)
{
    if (class(src) != "CharBuffer")
        stop("'src' is not a \"CharBuffer\" object")
    if (!is.integer(i))
        i <- as.integer(i)
    if (length(i) == 1) {
        if (length(imax) == 0)
            imax <- i
        else
            imax <- as.integer(imax)
        if (is.null(hash))
            .Call("CharBuffer_copy_from_i1i2", dest@xp, src@xp,
                  i, imax, PACKAGE="Biostrings")
        else
            .Call("CharBuffer_translate_copy_from_i1i2", dest@xp, src@xp,
                  i, imax, hash@xp, PACKAGE="Biostrings")
    } else {
        if (is.null(hash))
            .Call("CharBuffer_copy_from_subset", dest@xp, src@xp,
                  i, PACKAGE="Biostrings")
        else
            .Call("CharBuffer_translate_copy_from_subset", dest@xp, src@xp,
                  i, hash@xp, PACKAGE="Biostrings")
    }
    dest
}

CharBuffer.reverseCopy <- function(dest, i, imax=integer(0), src, hash=NULL)
{
    if (class(src) != "CharBuffer")
        stop("'src' is not a \"CharBuffer\" object")
    if (length(i) != 1)
        stop("'i' must be a single integer")
    if (!is.integer(i))
        i <- as.integer(i)
    if (length(imax) == 0)
        imax <- i
    else
        imax <- as.integer(imax)
    if (is.null(hash))
        .Call("CharBuffer_reverse_copy_from_i1i2", dest@xp, src@xp, i, imax, PACKAGE="Biostrings")
    else
        .Call("CharBuffer_reverse_translate_copy_from_i1i2", dest@xp, src@xp, i, imax, hash@xp, PACKAGE="Biostrings")
    dest
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### length(as.integer(cb)) is equivalent to length(cb)
### but the latter is MUCH faster!
setMethod("as.integer", "CharBuffer",
    function(x)
    {
        CharBuffer.readInts(x, 1, length(x))
    }
)

### Typical use:
###   cb <- CharBuffer(15)
###   cb[] <- 65
###   toString(cb)
###   cb[] <- "Hello"
###   toString(cb)
### So this should always rewrite the content of a "CharBuffer" object
### to itself, without any modification:
###   cb[] <- toString(cb)
### whatever the content of cb is!
setMethod("toString", "CharBuffer",
    function(x)
    {
        CharBuffer.read(x, 1, length(x))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting

### Select bytes from a "CharBuffer" object.
### Typical use:
###   cb <- CharBuffer(30)
###   cb[25:20]
###   cb[25:31] # subscript out of bounds
### Note: cb[] can be used as a shortcut for as.integer(cb)
setMethod("[", "CharBuffer",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            return(as.integer(x))
        CharBuffer.readInts(x, i)
    }
)

### Replace bytes in a "CharBuffer" object.
### Typical use:
###   cb <- CharBuffer(30)
###   cb[] <- 12 # fill with 12
###   cb[3:10] <- 1:-2
###   cb[3:10] <- "Ab"
###   cb[0] <- 4 # subscript out of bounds
###   cb[31] <- 4 # subscript out of bounds
###   cb[3] <- -12 # subscript out of bounds
setReplaceMethod("[", "CharBuffer",
    function(x, i, j,..., value)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")

        ## 'value' is a string
        if (is.character(value)) {
            if (length(value) >= 2)
                stop("character vector 'value' has more than one string")
            if (missing(i))
                return(CharBuffer.write(x, 1, length(x), value=value))
            return(CharBuffer.write(x, i, value=value))
        }

        ## We want to allow this: cb[3] <- 4, even if storage.mode(value)
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
            return(CharBuffer.writeInts(x, 1, length(x), value=value))
        CharBuffer.writeInts(x, i, value=value)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Equality

### Be careful to the semantic of the "==" operator:
###   2 "CharBuffer" objects are equals if their @xp slot is the
###   same "externalptr" instance (then they obviously have
###   the same length and contain the same data).
### With this definition, cb1 and cb2 can be 2 different "CharBuffer" objects
### (cb1 != cb2) and contain the same data.
setMethod("==", signature(e1="CharBuffer", e2="CharBuffer"),
    function(e1, e2)
    {
        address(e1@xp) == address(e2@xp)
    }
)
setMethod("!=", signature(e1="CharBuffer", e2="CharBuffer"),
    function(e1, e2)
    {
        address(e1@xp) != address(e2@xp)
    }
)

### A wrapper to the very fast memcmp() C-function.
### Arguments MUST be the following or it will crash R:
###   x1, x2: "CharBuffer" objects
###   start1, start2, width: single integers
### In addition: 1 <= start1 <= start1+width-1 <= length(x1)
###              1 <= start2 <= start2+width-1 <= length(x2)
### WARNING: This function is voluntarly unsafe (it doesn't check its
### arguments) because we want it to be the fastest possible!
CharBuffer.compare <- function(x1, start1, x2, start2, width)
{
    .Call("CharBuffer_memcmp", x1@xp, start1, x2@xp, start2, width, PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "PrintableCharBuffer" class is a simple extention of the "CharBuffer"
### class (no additional slots)

setClass("PrintableCharBuffer", representation("CharBuffer"))

setMethod("show", "PrintableCharBuffer",
    function(object)
    {
        print(toString(object))
    }
)

### Safe alternative to 'strsplit(x, NULL, fixed=TRUE)[[1]]'.
safeExplode <- function(x)
{
    if (!is.character(x) || length(x) != 1)
        stop("'x' must be a single string")
    .Call("safe_explode", x, PACKAGE="Biostrings")
}

### pcb <- new("PrintableCharBuffer", 10)
### pcb[] <- "ab-C."
### pcb[] == strsplit(toString(pcb), NULL, fixed=TRUE)[[1]]
setMethod("[", "PrintableCharBuffer",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            s <- toString(x)
        else
            s <- CharBuffer.read(x, i)
        safeExplode(s)
    }
)
