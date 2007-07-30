### Some low-level helper functions

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

