### =========================================================================
### Utility functions for creating or modifying IRanges objects
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "intToRanges" function.
###

intToRanges <- function(x, use.names=TRUE)
{
    if (!is.numeric(x))
        stop("'x' must be an integer vector")
    if (!is.integer(x))
        x <- as.integer(x)
    use.names <- normalize.use.names(use.names)
    if (use.names) ans_names <- names(x) else ans_names <- NULL
    new("IRanges", start=rep.int(1L, length(x)), width=x,
        names=ans_names, check=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "intToAdjacentRanges" function.
###

intToAdjacentRanges <- function(x, use.names=TRUE)
{
    if (!is.numeric(x))
        stop("'x' must be an integer vector")
    if (!is.integer(x))
        x <- as.integer(x)
    use.names <- normalize.use.names(use.names)
    ans_start <- .Call("int_to_adjacent_ranges", x, PACKAGE="Biostrings")
    if (use.names) ans_names <- names(x) else ans_names <- NULL
    new("IRanges", start=ans_start, width=x,
        names=ans_names, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "restrict" generic and methods.
###

setGeneric("restrict", signature="x",
    function(x, start, end, keep.nonoverlapping=FALSE, use.names=TRUE)
        standardGeneric("restrict")
)

.restrict.IRanges <- function(x, start, end,
                                   keep.nonoverlapping=FALSE, use.names=TRUE)
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
    if (!isTRUEorFALSE(keep.nonoverlapping))
        stop("'keep.nonoverlapping' must be 'TRUE' or 'FALSE'")
    use.names <- normalize.use.names(use.names)

    ans_start <- start(x)
    ans_end <- end(x)
    if (use.names) ans_names <- names(x) else ans_names <- NULL

    far_too_right <- end < ans_start
    far_too_left <- ans_end < start
    if (keep.nonoverlapping) {
        ans_start[far_too_right] <- end + 1L
        ans_end[far_too_left] <- start - 1L
    } else {
        keep <- !(far_too_right | far_too_left)
        ans_start <- ans_start[keep]
        ans_end <- ans_end[keep]
        if (!is.null(ans_names))
            ans_names <- ans_names[keep]
    }
    ## "fix" ans_start
    too_left <- ans_start < start
    ans_start[too_left] <- start
    ## "fix" ans_end
    too_right <- end < ans_end
    ans_end[too_right] <- end

    ans_width <- ans_end - ans_start + 1L

    update(x, start=ans_start, width=ans_width, names=ans_names, check=FALSE)
}

setMethod("restrict", "IRanges", .restrict.IRanges)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "narrow" generic and methods.
###

setGeneric("narrow", signature="x",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
        standardGeneric("narrow")
)

.narrow.IRanges <- function(x, start=NA, end=NA, width=NA, use.names=TRUE)
{
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
    use.names <- normalize.use.names(use.names)

    C_ans <- .Call("narrow_IRanges",
                   x, start, end, width,
                   PACKAGE="Biostrings")
    if (use.names) ans_names <- names(x) else ans_names <- NULL

    update(x, start=C_ans$start, width=C_ans$width, names=ans_names, check=FALSE)
}

setMethod("narrow", "IRanges", .narrow.IRanges)

setMethod("narrow", "numeric",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        y <- intToRanges(x, use.names=use.names)
        narrow(y, start=start, end=end, width=width, use.names=TRUE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "reduce" generic and methods.
###

setGeneric("reduce", signature="x",
    function(x, with.inframe.attrib=FALSE) standardGeneric("reduce")
)

.reduce.IRanges <- function(x, with.inframe.attrib=FALSE)
{
    if (!isTRUEorFALSE(with.inframe.attrib))
        stop("'with.inframe.attrib' must be 'TRUE' or 'FALSE'")
    C_ans <- .Call("reduce_IRanges",
                    x, with.inframe.attrib,
                    PACKAGE="Biostrings")
    ans <- update(x, start=C_ans$start, width=C_ans$width, check=FALSE)
    if (with.inframe.attrib) {
        inframe <- new("IRanges", start=C_ans$inframe.start,
                       width=width(x), check=FALSE)
        attr(ans, "inframe") <- inframe
    }
    ans
}

setMethod("reduce", "IRanges", .reduce.IRanges)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "mask" generic and methods.
###

setGeneric("mask", signature="x",
    function(x, start, end, ...) standardGeneric("mask")
)

.mask.IRanges <- function(x, start, end, ...)
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
    ## Ranges in 'y' are ordered from left to right and separated by gaps
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

setMethod("mask", "IRanges", .mask.IRanges)

