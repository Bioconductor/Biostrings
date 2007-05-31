### Robert contribution
readFASTA <- function(file, checkComments=TRUE)
{
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
               quote="", allowEscapes=FALSE);
    if (checkComments) {
        ##comments are supposedly lines beginning with semi-colons
        comments = grep("^;", s1)
        if (length(comments) > 0)
            s1 = s1[-comments]
    }
    descriptions <- grep("^>", s1)
    numF <- length(descriptions)
    if (numF == 0)
        stop("no FASTA sequences found")
    dp <- descriptions + 1
    dm <- descriptions - 1
    end = c(dm[-1], length(s1))
    ans = NULL
    for(i in 1:numF)
        ans[[i]] <- list(desc=s1[descriptions[i]],
                         seq=paste(s1[dp[i]:end[i]], collapse=""))
    ans
}

writeFASTA <- function(x, file, width=80)
{
    if (is.character(file)) {
        file <- file(file, "w")
        on.exit(close(file))
    } else {
        if (!inherits(file, "connection"))
            stop("'file' must be a character string or connection")
        if (!isOpen(file)) {
            open(file, "w")
            on.exit(close(file))
        }
    }

    for (rec in x) {
        cat(rec$desc, "\n", file=file, sep="")
        nlines <- nchar(rec$seq) %/% width + 1
        for (i in seq_len(nlines)) {
            start <- (i-1) * width + 1
            stop <- start + width - 1
            if (stop > nchar(rec$seq))
                stop <- nchar(rec$seq)
            line <- substr(rec$seq, start, stop)
            cat(line, "\n", file=file, sep="")
        }
    }
}

