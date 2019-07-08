### Some miscellaneous stuff

N50 <- function(csizes)
{
    if (!is.numeric(csizes))
        stop("'csizes' must be a vector containing the contig sizes")
    decreasing_csizes <- sort(csizes, decreasing=TRUE)
    tmp <- cumsum(decreasing_csizes)
    total_size <- tmp[length(tmp)]
    N50 <- decreasing_csizes[which(tmp >= total_size/2)[1]]
    return(N50)
}

