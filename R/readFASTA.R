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
               if (end[i] >= dp[i]) {
                   seq <- paste(s1[dp[i]:end[i]], collapse="")
               } else {
                   warning("record \"", desc, "\" contains no sequence")
                   seq <- ""
               }
               list(desc=desc, seq=seq)
           }
    )
}

writeFASTA <- function(x, file="", desc=NULL, append=FALSE, width=80)
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
    writeBString <- function(bstring)
    {
        if (length(bstring) == 0L)
            return()
        nlines <- (length(bstring) - 1L) %/% width + 1L
        lineIdx <- seq_len(nlines)
        start <- (lineIdx - 1L) * width + 1L
        end <- start + width - 1L
        if (end[nlines] > length(bstring))
            end[nlines] <- length(bstring)
        bigstring <- paste(
            as.character(Views(bstring, start = start, end = end)),
            collapse="\n"
        )
        cat(bigstring, "\n", file=file, sep="")
    }
    if (is.null(desc)) {
        for (rec in x) {
            cat(">", rec$desc, "\n", file=file, sep="")
            writeBString(BString(rec$seq))
        }
    } else {
        if (!is.character(desc))
            stop("when specified, 'desc' must be a character vector")
        if (length(desc) > length(x))
            stop("'desc' is longer than the number of sequences in 'x'")
        for (i in seq_len(length(desc))) {
            cat(">", desc[i], "\n", file=file, sep="")
            xi <- x[[i]]
            if (is.list(xi))
                xi <- xi$seq
            writeBString(BString(xi))
        }
    }
}

