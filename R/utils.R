### Some low-level helper functions

isLooseNumeric <- function(x)
{
    return(is.numeric(x) || (!is.null(x) && all(is.na(x))))
}

recycleVector <- function(x, length)
{
    y <- vector(storage.mode(x), length)
    y[] <- x
    y
}

