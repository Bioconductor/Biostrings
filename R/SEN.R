### =========================================================================
### The SEN (Start/End/Nchar) interface
### -------------------------------------------------------------------------

SEN2safelocs <- function(start, end, nchar, seq_nchars, check=TRUE)
{
    if (check) {
        ## Only limited checking here, more is done at the C level
        if (!isSingleNumberOrNA(start))
            stop("'start' must be a single integer or NA")
        if (!is.integer(start))
            start <- as.integer(start)
        if (!isSingleNumberOrNA(end))
            stop("'end' must be a single integer or NA")
        if (!is.integer(end))
            end <- as.integer(end)
        if (!isSingleNumberOrNA(nchar))
            stop("'nchar' must be a single integer or NA")
        if (!is.integer(nchar))
            nchar <- as.integer(nchar)
    }
    safelocs <- .Call("SEN_to_safelocs", start, end, nchar, seq_nchars, PACKAGE="Biostrings")
    new("SeqLocs", start=safelocs$start, nchar=safelocs$nchar, check=FALSE)
}

getStartForAdjacentSeqs <- function(seq_nchars)
{
    .Call("get_start_for_adjacent_seqs", seq_nchars, PACKAGE="Biostrings")
}

