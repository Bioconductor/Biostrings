### Some low-level helper functions

isTRUEorFALSE <- function(x)
{
    return(is.logical(x) && length(x) == 1 && !is.na(x))
}

isNumericOrNAs <- function(x)
{
    return(is.numeric(x) || (is.atomic(x) && is.vector(x) && all(is.na(x))))
}

recycleVector <- function(x, length)
{
    y <- vector(storage.mode(x), length)
    y[] <- x
    y
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

