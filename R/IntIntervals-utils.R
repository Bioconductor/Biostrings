### =========================================================================
### Utility functions for creating or modifying IntIntervals objects
### -------------------------------------------------------------------------
###

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
    C_ans <- .Call("narrow_IntIntervals",
                   x, start, end, width,
                   PACKAGE="Biostrings")
    new("IntIntervals", start=C_ans$start, width=C_ans$width,
        names=ans_names, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "reduce" function.
###

reduce <- function(x, with.inframe.attrib=FALSE)
{
    if (!is(x, "IntIntervals"))
        stop("'x' must be an IntIntervals object")
    if (!isTRUEorFALSE(with.inframe.attrib))
        stop("'with.inframe.attrib' must be 'TRUE' or 'FALSE'")
    C_ans <- .Call("reduce_IntIntervals",
                    x, with.inframe.attrib,
                    PACKAGE="Biostrings")
    ans <- new("IntIntervals", start=C_ans$start, width=C_ans$width, check=FALSE)
    if (with.inframe.attrib) {
        inframe <- new("IntIntervals", start=C_ans$inframe.start,
                       width=width(x), check=FALSE)
        attr(ans, "inframe") <- inframe
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "mask" method.
###

setGeneric("mask", signature="x",
    function(x, start, end, ...) standardGeneric("mask")
)

.mask.IntIntervals <- function(x, start, end, ...)
{
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
    y <- restrict(reduce(x), start, end)
    y <- y[width(y) != 0]
    ## Intervals in 'y' are ordered from left to right and separated by gaps
    ans_start <- ans_end <- integer(0)
    start0 <- start
    for (i in seq_len(length(y))) {
        end0 <- start(y)[i] - 1L
        start1 <- end(y)[i] + 1L
        if (end0 >= start0) {
            ans_start <- c(ans_start, start0)
            ans_end <- c(ans_end, end0)
        }
        if (start1 > start0)
            start0 <- start1
    }
    if (start0 <= end) {
        ans_start <- c(ans_start, start0)
        ans_end <- c(ans_end, end)
    }
    ans_width <- ans_end - ans_start + 1L
    update(x, start=ans_start, width=ans_width, check=FALSE)
}

setMethod("mask", "IntIntervals", .mask.IntIntervals)

