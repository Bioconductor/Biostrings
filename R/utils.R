### Some low-level (not exported) helper functions

### TODO: Clean this file (some stuff is not needed anymore).

extraArgsAsList <- function(.valid.argnames, ...)
{
    args <- list(...)
    argnames <- names(args)
    if (length(args) != 0
        && (is.null(argnames) || any(argnames %in% c("", NA))))
        stop("all extra arguments must be named")
    if (!is.null(.valid.argnames) && !all(argnames %in% .valid.argnames))
        stop("valid extra argument names are ",
             paste("'", .valid.argnames, "'", sep="", collapse=", "))
    if (any(duplicated(argnames)))
        stop("argument names must be unique")
    args
}

isTRUEorFALSE <- function(x)
{
    is.logical(x) && length(x) == 1 && !is.na(x)
}

isSingleInteger <- function(x)
{
    is.integer(x) && length(x) == 1 && !is.na(x)
}

isSingleNumber <- function(x)
{
    is.numeric(x) && length(x) == 1 && !is.na(x)
}

isSingleNumberOrNA <- function(x)
{
    is.vector(x) && is.atomic(x) && length(x) == 1 && (is.numeric(x) || is.na(x))
}

isNumericOrNAs <- function(x)
{
    is.numeric(x) || (is.atomic(x) && is.vector(x) && all(is.na(x)))
}

isSingleString <- function(x)
{
    is.character(x) && length(x) == 1 && !is.na(x)
}

isSingleStringOrNA <- function(x)
{
    is.vector(x) && is.atomic(x) && length(x) == 1 && (is.character(x) || is.na(x))
}

numeric2integer <- function(x)
{
    if (is.numeric(x) && !is.integer(x)) as.integer(x) else x
}

stopIfProblems <- function(problems)
{
    if (!is.null(problems)) stop(paste(problems, collapse="\n  "))
}

normargIntegerOrNA <- function(x, name)
{
    if (!isNumericOrNAs(x))
        stop(paste("'", name, "' must be an integer vector", sep = ""))
    if (!is.integer(x))
        x <- as.integer(x)
    x
}

normargSingleStart <- function(start)
{
    if (!isSingleNumber(start))
        stop("'start' must be a single integer")
    if (!is.integer(start))
        start <- as.integer(start)
    start
}

normargSingleEnd <- function(end)
{
    if (!isSingleNumber(end))
        stop("'end' must be a single integer")
    if (!is.integer(end))
        end <- as.integer(end)
    end
}

normargSingleStartOrNA <- function(start)
{
    if (!isSingleNumberOrNA(start))
        stop("'start' must be a single integer or NA")
    if (!is.integer(start))
        start <- as.integer(start)
    start
}

normargSingleEndOrNA <- function(end)
{
    if (!isSingleNumberOrNA(end))
        stop("'end' must be a single integer or NA")
    if (!is.integer(end)) 
        end <- as.integer(end)
    end
}

normargSingleWidthOrNA <- function(width)
{
    if (!isSingleNumberOrNA(width))
        stop("'width' must be a single integer or NA")
    if (!is.integer(width))    
        width <- as.integer(width)
    width
}

normargStart <- function(start)
{
    if (!isSingleNumber(start))
        stop("'start' must be a single integer")
    if (!is.integer(start))
        start <- as.integer(start)
    if (start < 1L)
        stop("'start' must be >= 1")
    start
}

normargNchar <- function(start, nchar, seq_nchar)
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

normargUseNames <- function(use.names)
{
    if (is.null(use.names))
        return(TRUE)
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    use.names
}

### Returns an integer vector.
pow.int <- function(x, y)
{
    if (!is.numeric(x))
        stop("'x' must be a numeric vector")
    if (!is.integer(x))
        x <- as.integer(x)
    ans <- rep.int(1L, length(x))
    for (i in seq_len(y))
        ans <- ans * x
    ans
}

