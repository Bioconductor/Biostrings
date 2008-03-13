### Some low-level (not exported) helper functions

isTRUEorFALSE <- function(x)
{
    is.logical(x) && length(x) == 1 && !is.na(x)
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

normalize.use.names <- function(use.names)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be 'TRUE' or 'FALSE'")
    use.names
}

recycleVector <- function(x, length)
{
    y <- vector(storage.mode(x), length)
    y[] <- x
    y
}

### Returns
pow.int <- function(x, y)
{
    if (!is.numeric(x))
        stop("'x' must be a numeric vector")
    x <- as.integer(x)
    ans <- rep.int(1L, length(x))
    for (i in seq_len(y))
        ans <- ans * x
    ans
}

