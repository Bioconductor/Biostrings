### Robert's contribution
readFASTA <- function(file, checkComments=TRUE, strip.descs=TRUE)
{
    if (missing(strip.descs))
        warning("use 'strip.descs=FALSE' for compatibility with old version\n",
                "  of readFASTA(), or 'strip.descs=TRUE' to remove the \">\"\n",
                "  at the beginning of the description lines and to get\n",
                "  rid of this warning (see '?readFASTA' for more details)")
    if (is.character(file)) {
        file <- file(file, "r")
        on.exit(close(file))
    } else {
        if (!inherits(file, "connection"))
            stop("'file' must be a character string or connection")
        if (!isOpen(file)) {
            open(file, "r")
            on.exit(close(file))
        }
    }

    s1 <- scan(file=file, what="", sep="\n",
               quote="", allowEscapes=FALSE, quiet=TRUE)
    if (checkComments) {
        ##comments are supposedly lines beginning with semi-colons
        comments <- grep("^;", s1)
        if (length(comments) > 0)
            s1 <- s1[-comments]
    }
    descriptions <- which(substr(s1, 1L, 1L) == ">")
    numF <- length(descriptions)
    if (numF == 0)
        stop("no FASTA sequences found")
    dp <- descriptions + 1L
    dm <- descriptions - 1L
    end <- c(dm[-1], length(s1))
    lapply(seq_len(numF),
           function(i)
           {
               desc <- s1[descriptions[i]]
               if (strip.descs)
                   desc <- substr(desc, 2L, nchar(desc))
               seq <- paste(s1[dp[i]:end[i]], collapse="")
               list(desc=desc, seq=seq)
           }
    )
}

writeFASTA <- function(x, file="", append=FALSE, width=80)
{
    if (!isTRUEorFALSE(append))
        stop("'append' must be TRUE or FALSE")
    if (isSingleString(file)) {
        if (file == "") {
            file <- stdout()
        } else {
            file <- file(file, ifelse(append, "a", "w"))
            on.exit(close(file))
        }
    } else if (inherits(file, "connection")) {
        if (!isOpen(file)) {
            file <- file(file, ifelse(append, "a", "w"))
            on.exit(close(file))
        }
    } else {
        stop("'file' must be a single string or connection")
    }
    if (!isSingleNumber(width))
        stop("'width' must be an integer >= 1")
    if (!is.integer(width))
        width <- as.integer(width)
    if (width < 1L)
        stop("'width' must be an integer >= 1")
    for (rec in x) {
        cat(">", rec$desc, "\n", file=file, sep="")
        nlines <- nchar(rec$seq) %/% width + 1L
        for (i in seq_len(nlines)) {
            start <- (i-1L) * width + 1L
            stop <- start + width - 1L
            if (stop > nchar(rec$seq))
                stop <- nchar(rec$seq)
            line <- substr(rec$seq, start, stop)
            cat(line, "\n", file=file, sep="")
        }
    }
}

