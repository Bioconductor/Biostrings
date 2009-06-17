### Some miscellaneous stuff

N50 <- function(csizes)
{
    if (!is.numeric(csizes))
        stop("'csizes' must be a vector containing the contig sizes")
    if (!is.integer(csizes))
        csizes <- as.integer(csizes)
    sorted_csizes <- sort(csizes)
    tmp <- cumsum(sorted_csizes)
    N50 <- sorted_csizes[which(tmp >= max(tmp)/2)[1]]
    return(N50)
}

