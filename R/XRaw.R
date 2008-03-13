### =========================================================================
### External raw vectors: the "XRaw" class
### -------------------------------------------------------------------------
###
### An "XRaw" object is a chunk of memory that:
###   a. Contains bytes (char at the C level).
###   b. Is readable and writable.
###   c. Has a passed by address semantic i.e. the data it contains are not
###      copied on object duplication. For example when doing
###        xr2 <- xr1
###      both xr1 and xr2 point to the same place in memory.
###      This is achieved by using R predefined type "externalptr".
###
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


### Return the hexadecimal address of any R object in a string.
address <- function(x)
{
    .Call("Biostrings_sexp_address", x, PACKAGE="Biostrings")
}

### Helper function (for debugging purpose).
### Print some obscure info about an "externalptr" object.
### Typical use:
###   show(new("externalptr"))
setMethod("show", "externalptr",
    function(object)
    {
        .Call("Biostrings_xp_show", object, PACKAGE="Biostrings")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "XRaw" class.
### Note: instead of defining the "XRaw" class with just one slot of type
### "externalptr" (it HAS an "externalptr", and nothing else), an alternative
### would be to simply extend the "externalptr" type.
### After all, a "XRaw" object IS an "externalptr" object.
### However, I tried this but was not able to implement the "initialize" method
### in such a way that it returns a new instance of the "XRaw" class (the
### returned object was ALWAYS the same instance everytime the method was
### called, I found no workaround).
###

setClass("XRaw", representation(xp="externalptr"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
###
### Note that unlike raw vectors, XRaw objects are not initialized with 0's.
###

### This:
###   xr <- XRaw(30)
### will call this "initialize" method.
setMethod("initialize", "XRaw",
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
        .Call("Biostrings_XRaw_alloc", xp, length, PACKAGE="Biostrings")
        if (verbose) {
            cat(" OK\n")
            show_string <- .Call("Biostrings_XRaw_get_show_string", xp, PACKAGE="Biostrings")
            cat("New", show_string, "successfully created\n")
        }
        .Object@xp <- xp
        .Object
    }
)

XRaw <- function(...)
{
    new("XRaw", ...)
}

setMethod("show", "XRaw",
    function(object)
    {
        show_string <- .Call("Biostrings_XRaw_get_show_string", object@xp, PACKAGE="Biostrings")
        cat(show_string, "\n", sep="")
        ## What is correct here? The documentation (?show) says that 'show'
        ## should return an invisible 'NULL' but, on the other hand, the 'show'
        ## method for intergers returns its 'object' argument...
        invisible(object)
    }
)

setMethod("length", "XRaw",
    function(x)
    {
        .Call("Biostrings_XRaw_length", x@xp, PACKAGE="Biostrings")
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Read/write functions.
### These are almost safe wrappers to unsafe C functions ("almost" because
### they don't check for NAs in their arguments).
### If length(i) == 0 then the read functions return an empty vector
### and the write functions don't do anything.

XRaw.readInts <- function(x, i, imax=integer(0))
{
    if (!is.integer(i))
        i <- as.integer(i)
    if (length(i) == 1) {
        if (length(imax) == 0)
            imax <- i
        else
            imax <- as.integer(imax)
        .Call("XRaw_read_ints_from_i1i2", x@xp, i, imax, PACKAGE="Biostrings")
    } else {
        .Call("XRaw_read_ints_from_subset", x@xp, i, PACKAGE="Biostrings")
    }
}

XRaw.writeInts <- function(x, i, imax=integer(0), value)
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
        .Call("XRaw_write_ints_to_i1i2", x@xp, i, imax, value, PACKAGE="Biostrings")
    } else {
        .Call("XRaw_write_ints_to_subset", x@xp, i, value, PACKAGE="Biostrings")
    }
    x
}

### 'dec_lkup' must be NULL or a vector of integers
XRaw.read <- function(x, i, imax=integer(0), dec_lkup=NULL)
{
    if (!is.integer(i))
        i <- as.integer(i)
    if (length(i) == 1) {
        if (length(imax) == 0)
            imax <- i
        else
            imax <- as.integer(imax)
        if (is.null(dec_lkup))
            .Call("Biostrings_XRaw_read_chars_from_i1i2",
                  x@xp, i, imax, PACKAGE="Biostrings")
        else
            .Call("XRaw_read_enc_chars_from_i1i2",
                  x@xp, i, imax, dec_lkup, PACKAGE="Biostrings")
    } else {
        if (is.null(dec_lkup))
            .Call("Biostrings_XRaw_read_chars_from_subset",
                  x@xp, i, PACKAGE="Biostrings")
        else
            .Call("XRaw_read_enc_chars_from_subset",
                  x@xp, i, dec_lkup, PACKAGE="Biostrings")
    }
}

### 'enc_lkup' must be NULL or a vector of integers
XRaw.write <- function(x, i, imax=integer(0), value, enc_lkup=NULL)
{
    if (!isSingleString(value))
        stop("'value' must be a single string")
    if (!is.integer(i))
        i <- as.integer(i)
    if (length(i) == 1) {
        if (length(imax) == 0)
            imax <- i
        else
            imax <- as.integer(imax)
        if (is.null(enc_lkup))
            .Call("Biostrings_XRaw_write_chars_to_i1i2",
                  x@xp, i, imax, value, PACKAGE="Biostrings")
        else
            .Call("XRaw_write_enc_chars_to_i1i2",
                  x@xp, i, imax, value, enc_lkup, PACKAGE="Biostrings")
    } else {
        if (is.null(enc_lkup))
            .Call("Biostrings_XRaw_write_chars_to_subset",
                  x@xp, i, value, PACKAGE="Biostrings")
        else
            .Call("XRaw_write_enc_chars_to_subset",
                  x@xp, i, value, enc_lkup, PACKAGE="Biostrings")
    }
    x
}

### 'lkup' must be NULL or a vector of integers
XRaw.copy <- function(dest, i, imax=integer(0), src, lkup=NULL)
{
    if (class(src) != "XRaw")
        stop("'src' is not a \"XRaw\" object")
    if (!is.integer(i))
        i <- as.integer(i)
    if (length(i) == 1) {
        if (length(imax) == 0)
            imax <- i
        else
            imax <- as.integer(imax)
        if (is.null(lkup))
            .Call("Biostrings_XRaw_copy_from_i1i2", dest@xp, src@xp,
                  i, imax, PACKAGE="Biostrings")
        else
            .Call("XRaw_translate_copy_from_i1i2", dest@xp, src@xp,
                  i, imax, lkup, PACKAGE="Biostrings")
    } else {
        if (is.null(lkup))
            .Call("Biostrings_XRaw_copy_from_subset", dest@xp, src@xp,
                  i, PACKAGE="Biostrings")
        else
            .Call("XRaw_translate_copy_from_subset", dest@xp, src@xp,
                  i, lkup, PACKAGE="Biostrings")
    }
    dest
}

### 'lkup' must be NULL or a vector of integers
XRaw.reverseCopy <- function(dest, i, imax=integer(0), src, lkup=NULL)
{
    if (class(src) != "XRaw")
        stop("'src' is not a \"XRaw\" object")
    if (length(i) != 1)
        stop("'i' must be a single integer")
    if (!is.integer(i))
        i <- as.integer(i)
    if (length(imax) == 0)
        imax <- i
    else
        imax <- as.integer(imax)
    if (is.null(lkup))
        .Call("XRaw_reverse_copy_from_i1i2", dest@xp, src@xp, i, imax, PACKAGE="Biostrings")
    else
        .Call("XRaw_reverse_translate_copy_from_i1i2", dest@xp, src@xp, i, imax, lkup, PACKAGE="Biostrings")
    dest
}

### 'lkup' must be a vector of complexes
XRaw.readComplexes <- function(x, i, imax=integer(0), lkup)
{
    if (!is.integer(i))
        i <- as.integer(i)
    if (length(i) == 1) {
        if (length(imax) == 0)
            imax <- i
        else
            imax <- as.integer(imax)
        .Call("XRaw_read_complexes_from_i1i2",
              x@xp, i, imax, lkup, PACKAGE="Biostrings")
    } else {
        .Call("XRaw_read_complexes_from_subset",
              x@xp, i, lkup, PACKAGE="Biostrings")
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 2 convenience wrappers
###

normalize.start <- function(start)
{
    if (!isSingleNumber(start))
        stop("'start' must be a single integer")
    if (!is.integer(start))
        start <- as.integer(start)
    if (start < 1L)
        stop("'start' must be >= 1")
    start
}

normalize.nchar <- function(start, nchar, seq_nchar)
{
    if (!isSingleNumberOrNA(nchar))
        stop("'nchar' must be a single integer or NA")
    if (is.na(nchar)) {
        nchar <- seq_nchar - start + 1L
        if (nchar < 0L)
            stop("cannot read a negative number of letters")
        return(nchar)
    }
    if (!is.integer(nchar))
        nchar <- as.integer(nchar)
    if (nchar < 0L)
        stop("cannot read a negative number of letters")
    end <- start + nchar - 1L
    if (end > seq_nchar)
        stop("cannot read beyond the end of 'seq'")
    nchar
}

charToXRaw <- function(x, start=NA, end=NA, width=NA, collapse=NULL, lkup=NULL, check=TRUE)
{
    safe_locs <- narrow(nchar(x, type="bytes"), start, end, width)
    .Call("new_XRaw_from_STRSXP",
          x, start(safe_locs), width(safe_locs), collapse, lkup,
          PACKAGE="Biostrings")
}

copySubXRaw <- function(x, start=1, nchar=NA, lkup=NULL, check=TRUE)
{
    if (check) {
        start <- normalize.start(start)
        nchar <- normalize.nchar(start, nchar, length(x))
    }
    ans <- XRaw(nchar)
    XRaw.copy(ans, start, start + nchar - 1L, src=x, lkup=lkup)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### XRaw IO
###

XRaw.saveFASTA <- function(x, filepath, dec_lkup=NULL)
{
    stop("not ready yet, sorry!")
}

### Return a list of 4 elements (see comments for XRaw_loadFASTA() in
### src/XRaw_fillread.c for the details).
### 'filepath' must a path to an uncompressed FASTA file. Note that,
### unlike with the file() function, it cannot an URL, '""', '"stdin"'
### or '"clipboard"'.
XRaw.loadFASTA <- function(x, filepath, collapse="", enc_lkup=NULL)
{
    if (!isSingleString(filepath))
        stop("'filepath' must be a single string")
    if (!isSingleString(collapse))
        stop("'collapse' must be a single string")
    if (!is.null(enc_lkup) && !is.integer(enc_lkup))
        stop("'enc_lkup' must be an integer vector")
    filepath <- path.expand(filepath)
    .Call("XRaw_loadFASTA",
          x@xp, filepath, collapse, enc_lkup, PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### length(as.integer(xr)) is equivalent to length(xr)
### but the latter is MUCH faster!
setMethod("as.integer", "XRaw",
    function(x)
    {
        XRaw.readInts(x, 1, length(x))
    }
)

### Typical use:
###   xr <- XRaw(15)
###   xr[] <- 65
###   toString(xr)
###   xr[] <- "Hello"
###   toString(xr)
### So this should always rewrite the content of a "XRaw" object
### to itself, without any modification:
###   xr[] <- toString(xr)
### whatever the content of xr is!
setMethod("toString", "XRaw",
    function(x, ...)
    {
        XRaw.read(x, 1, length(x))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

### Select bytes from a "XRaw" object.
### Typical use:
###   xr <- XRaw(30)
###   xr[25:20]
###   xr[25:31] # subscript out of bounds
### Note: xr[] can be used as a shortcut for as.integer(xr)
setMethod("[", "XRaw",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            return(as.integer(x))
        XRaw.readInts(x, i)
    }
)

### Replace bytes in a "XRaw" object.
### Typical use:
###   xr <- XRaw(30)
###   xr[] <- 12 # fill with 12
###   xr[3:10] <- 1:-2
###   xr[3:10] <- "Ab"
###   xr[0] <- 4 # subscript out of bounds
###   xr[31] <- 4 # subscript out of bounds
###   xr[3] <- -12 # subscript out of bounds
setReplaceMethod("[", "XRaw",
    function(x, i, j,..., value)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")

        ## 'value' is a string
        if (is.character(value)) {
            if (missing(i))
                return(XRaw.write(x, 1, length(x), value=value))
            return(XRaw.write(x, i, value=value))
        }

        ## We want to allow this: xr[3] <- 4, even if storage.mode(value)
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
            return(XRaw.writeInts(x, 1, length(x), value=value))
        XRaw.writeInts(x, i, value=value)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Equality.
###

### Be careful to the semantic of the "==" operator:
###   2 "XRaw" objects are equals if their @xp slot is the
###   same "externalptr" instance (then they obviously have
###   the same length and contain the same data).
### With this definition, xr1 and xr2 can be 2 different "XRaw" objects
### (xr1 != xr2) and contain the same data.
setMethod("==", signature(e1="XRaw", e2="XRaw"),
    function(e1, e2)
    {
        address(e1@xp) == address(e2@xp)
    }
)
setMethod("!=", signature(e1="XRaw", e2="XRaw"),
    function(e1, e2)
    {
        address(e1@xp) != address(e2@xp)
    }
)

### A wrapper to the very fast memcmp() C-function.
### Arguments MUST be the following or it will crash R:
###   x1, x2: "XRaw" objects
###   start1, start2, width: single integers
### In addition: 1 <= start1 <= start1+width-1 <= length(x1)
###              1 <= start2 <= start2+width-1 <= length(x2)
### WARNING: This function is voluntarly unsafe (it doesn't check its
### arguments) because we want it to be the fastest possible!
XRaw.compare <- function(x1, start1, x2, start2, width)
{
    .Call("Biostrings_XRaw_memcmp", x1@xp, start1, x2@xp, start2, width, PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "PrintableXRaw" class is a simple extention of the "XRaw"
### class (no additional slots)

setClass("PrintableXRaw", representation("XRaw"))

setMethod("show", "PrintableXRaw",
    function(object)
    {
        print(toString(object))
    }
)

### Safe alternative to 'strsplit(x, NULL, fixed=TRUE)[[1]]'.
safeExplode <- function(x)
{
    if (!isSingleString(x))
        stop("'x' must be a single string")
    .Call("Biostrings_safe_explode", x, PACKAGE="Biostrings")
}

### pxr <- new("PrintableXRaw", 10)
### pxr[] <- "ab-C."
### pxr[] == strsplit(toString(pxr), NULL, fixed=TRUE)[[1]]
setMethod("[", "PrintableXRaw",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            s <- toString(x)
        else
            s <- XRaw.read(x, i)
        safeExplode(s)
    }
)
