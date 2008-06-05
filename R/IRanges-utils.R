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
    use.names <- normargUseNames(use.names)
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
    use.names <- normargUseNames(use.names)
    ans_start <- .Call("int_to_adjacent_ranges", x, PACKAGE="Biostrings")
    if (use.names) ans_names <- names(x) else ans_names <- NULL
    new("IRanges", start=ans_start, width=x,
                   names=ans_names, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "shift" generic and methods.
###
### Shifting preserves normality.
###

setGeneric("shift", signature="x",
    function(x, shift, use.names=TRUE) standardGeneric("shift")
)

setMethod("shift", "IRanges",
    function(x, shift, use.names=TRUE)
    {
        if (!isSingleNumber(shift))
            stop("'shift' must be a single integer")
        if (!is.integer(shift))
            shift <- as.integer(shift)
        if (!normargUseNames(use.names))
            names(x) <- NULL
        unsafe.start(x) <- start(x) + shift
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "reverse" generic and methods.
###

setGeneric("reverse", signature="x",
    function(x, ...) standardGeneric("reverse")
)

### This method does NOT preserve normality.
IRanges.reverse <- function(x, ...)
{
    args <- extraArgsAsList(NULL, ...)
    argnames <- names(args)
    n2p <- match(c("start", "end", "use.names"), argnames)
    if (is.na(n2p[1]))
        stop("'start' must be specified for \"reverse\" method for IRanges objects")
    start <- normargSingleStart(args[[n2p[1]]])
    if (is.na(n2p[2]))
        stop("'end' must be specified for \"reverse\" method for IRanges objects")
    end <- normargSingleEnd(args[[n2p[2]]])
    if (!is.na(n2p[3]) && !normargUseNames(args[[n2p[3]]]))
        unsafe.names(x) <- NULL
    unsafe.start(x) <- start + end - end(x)
    x
}

setMethod("reverse", "IRanges", IRanges.reverse)

setMethod("reverse", "NormalIRanges",
    function(x, ...)
    {
        ## callNextMethod() temporarily breaks 'x' as a NormalIRanges object
        ## because the returned ranges are ordered from right to left.
        x <- callNextMethod()
        unsafe.update(x, start=rev(start(x)), width=rev(width(x)), names=rev(names(x)))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "restrict" generic and methods.
###
### Note that when used with 'keep.all.ranges=FALSE', restrict() preserves
### normality.
###

setGeneric("restrict", signature="x",
    function(x, start=NA, end=NA, keep.all.ranges=FALSE, use.names=TRUE)
        standardGeneric("restrict")
)

setMethod("restrict", "IRanges",
    function(x, start=NA, end=NA, keep.all.ranges=FALSE, use.names=TRUE)
    {
        start <- normargSingleStartOrNA(start)
        end <- normargSingleEndOrNA(end)
        if (!isTRUEorFALSE(keep.all.ranges))
            stop("'keep.all.ranges' must be 'TRUE' or 'FALSE'")
        use.names <- normargUseNames(use.names)

        ans_start <- start(x)
        ans_end <- end(x)
        if (use.names) ans_names <- names(x) else ans_names <- NULL

        if (!is.na(start)) {
            far_too_left <- ans_end < start
            if (keep.all.ranges) {
                ans_end[far_too_left] <- start - 1L
            } else {
                keep_it <- !far_too_left
                ans_start <- ans_start[keep_it]
                ans_end <- ans_end[keep_it]
                if (!is.null(ans_names))
                    ans_names <- ans_names[keep_it]
            }
            ## "fix" ans_start
            too_left <- ans_start < start
            ans_start[too_left] <- start
        }
        if (!is.na(end)) {
            far_too_right <- end < ans_start
            if (keep.all.ranges) {
                ans_start[far_too_right] <- end + 1L
            } else {
                keep_it <- !far_too_right
                ans_start <- ans_start[keep_it]
                ans_end <- ans_end[keep_it]
                if (!is.null(ans_names))
                    ans_names <- ans_names[keep_it]
            }
            ## "fix" ans_end
            too_right <- end < ans_end
            ans_end[too_right] <- end
        }
        ans_width <- ans_end - ans_start + 1L

        unsafe.update(x, start=ans_start, width=ans_width, names=ans_names)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "narrow" generic and methods.
###
### Note that in general, narrow() does NOT preserve normality.
###

setGeneric("narrow", signature="x",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
        standardGeneric("narrow")
)

setMethod("narrow", "IRanges",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        start <- normargSingleStartOrNA(start)
        end <- normargSingleEndOrNA(end)
        width <- normargSingleWidthOrNA(width)
        use.names <- normargUseNames(use.names)

        C_ans <- .Call("narrow_IRanges",
                       x, start, end, width,
                       PACKAGE="Biostrings")
        if (use.names) ans_names <- names(x) else ans_names <- NULL

        unsafe.update(x, start=C_ans$start, width=C_ans$width, names=ans_names)
    }
)

setMethod("narrow", "NormalIRanges",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
        stop("narrowing a ", class(x), " instance is not supported")
)

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
### Note that reduce() preserves normality (of course).
###

setGeneric("reduce", signature="x",
    function(x, with.inframe.attrib=FALSE) standardGeneric("reduce")
)

setMethod("reduce", "IRanges",
    function(x, with.inframe.attrib=FALSE)
    {
        if (!isTRUEorFALSE(with.inframe.attrib))
            stop("'with.inframe.attrib' must be 'TRUE' or 'FALSE'")
        C_ans <- .Call("reduce_IRanges",
                        x, with.inframe.attrib,
                        PACKAGE="Biostrings")
        ans <- unsafe.update(x, start=C_ans$start, width=C_ans$width, names=NULL)
        if (with.inframe.attrib) {
            inframe <- new("IRanges", start=C_ans$inframe.start,
                                      width=width(x), check=FALSE)
            attr(ans, "inframe") <- inframe
        }
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "toNormalIRanges" function.
###

toNormalIRanges <- function(x)
{
    if (!is(x, "IRanges"))
        stop("'x' must be an IRanges object")
    x1 <- as(x, "IRanges") # downgrade
    x2 <- reduce(x1)
    x3 <- x2[width(x2) != 0]
    asNormalIRanges(x3, check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "gaps" generic and methods.
###
### Note that gaps() will always return a normal IRanges object (so, obviously,
### it preserves normality).
###

setGeneric("gaps", signature="x",
    function(x, start=NA, end=NA)
        standardGeneric("gaps")
)

setMethod("gaps", "IRanges",
    function(x, start=NA, end=NA)
    {
        start <- normargSingleStartOrNA(start)
        end <- normargSingleEndOrNA(end)
        ## No matter in what order restricting and normalizing are done, the final
        ## result should always be exactly the same.
        ## Now which order is the most efficient? It depends...
        xx <- toNormalIRanges(x)
        xx0 <- restrict(xx, start=start, end=end) # preserves normality
        ans_start <- ans_width <- integer(0)
        if (isEmpty(xx0)) {
            if (is.na(start) || is.na(end))
                stop("'x' is not overlapping with the unbounded region ",
                     "represented by 'start' and 'end'")
            if (start <= end) {
                ans_start <- start
                ans_width <- end - start + 1L
            }
        } else {
            start0 <- start(xx0)
            end0 <- end(xx0)
            if (!is.na(start) && start < min(xx0)) {
                start0 <- c(start, start0)
                end0 <- c(start - 1L, end0)
            }
            if (!is.na(end) && max(xx0) < end) {
                start0 <- c(start0, end + 1L)
                end0 <- c(end0, end)
            }
            if (length(start0) >= 2) {
                ans_start <- end0[-length(start0)] + 1L
                ans_width <- start0[-1] - ans_start
            } 
        }
        unsafe.update(x, start=ans_start, width=ans_width, names=NULL)
    }
)

