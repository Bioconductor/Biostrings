### Some low-level (not exported) helper functions

isTRUEorFALSE <- function(x)
{
    return(is.logical(x) && length(x) == 1 && !is.na(x))
}

isSingleNumber <- function(x)
{
    return(is.numeric(x) && length(x) == 1 && !is.na(x))
}

isNumericOrNAs <- function(x)
{
    return(is.numeric(x) || (is.atomic(x) && is.vector(x) && all(is.na(x))))
}

isSingleString <- function(x)
{
    return(is.character(x) && length(x) == 1 && !is.na(x))
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

debug_utils <- function()
{
    invisible(.Call("Biostrings_debug_utils", PACKAGE="Biostrings"))
}

debug_bufutils <- function()
{
    invisible(.Call("Biostrings_debug_bufutils", PACKAGE="Biostrings"))
}

debug_views_buffer <- function()
{
    invisible(.Call("Biostrings_debug_views_buffer", PACKAGE="Biostrings"))
}

