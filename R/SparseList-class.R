### =========================================================================
### SparseList objects
### -------------------------------------------------------------------------
###

setClass("SparseList",
    representation(
        length="integer",
        env="environment"
    )
)

### Typical use:
###   env <- new.env(hash=TRUE, parent=emptyenv())
###   key <- formatC(98L, width=10, format="d", flag="0")
###   value <- 3:-2
###   assign(key, value, envir=env)
###   x <- new("SparseList", length=100L, env=env)
###   length(x)
###   ls(x)
###   ls(x, all.names=TRUE)
###   as.list(x)
###   as.list(x, all.names=TRUE)
###   x[[1]]
###   x[[98]
###   x[[101]]]
###

setMethod("length", "SparseList", function(x) x@length)

### 'pos', 'envir' and 'pattern' args are ignored
setMethod("ls", signature(name="SparseList"),
    function(name, pos, envir, all.names=FALSE, pattern)
    {
        if (!all.names)
            return(ls(name@env, all.names=TRUE))
        seq_len(length(name))
    }
)

setMethod("as.list", "SparseList",
    function(x, all.names=FALSE, ...)
    {
        if (!all.names)
            return(as.list(x@env, all.names=TRUE))
        ans <- vector(mode="list", length=length(x))
        symbols <- ls(x)
        for (symb in symbols)
            ans[[as.integer(symb)]] <- get(symb, envir=x@env)
        ans
    }
)

### Supported 'i' types: character or numeric vector of length 1.
setMethod("[[", "SparseList",
    function(x, i, j, ...)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i))
            stop("subscript is missing")
        if (!is.character(i) && !is.numeric(i))
            stop("invalid subscript type")
        if (length(i) < 1L)
            stop("attempt to select less than one element")
        if (length(i) > 1L)
            stop("attempt to select more than one element")
        if (is.na(i))
            stop("subscript cannot be NA")
        if (is.character(i))
            return(get(i, envir=x@env))
        if (!is.integer(i))
            i <- as.integer(i)
        if (i < 1L || i > length(x))
            stop("subscript out of bounds")
        i <- formatC(i, width=10, format="d", flag="0")
        get(i, envir=x@env)
    }
)

