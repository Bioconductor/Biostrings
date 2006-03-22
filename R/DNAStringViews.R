# Similar to the DNAString() constructor, but:
#   - return a BStringViews object instead of a BString object,
#   - it is vectorized.

# For now: just a naive implementation.
# Test:
#   alphabet <- strsplit("-TGCANBDHKMRSVWY", NULL)[[1]]
#   n <- 40000
#   src <- sapply(1:n, function(i) {paste(sample(alphabet, 250, replace=TRUE), collapse="")})
#   dsv <- DNAStringViews(src)
#
# Comparing DNAStringViews() speed vs "old" DNAString() speed:
#
#       n  DNAStringViews  "old" DNAString
#   -----  --------------  ---------------
#    5000          0.26 s           4.15 s
#   10000          0.51 s          16.29 s
#   20000          0.99 s          64.85 s
#   40000          1.69 s         488.43 s
#
# Well, the naive implementation is not too bad!
# (the quadratic behaviour of "old" DNAString() was first reported
# by Wolfgang)


# 'width' is the vector of view widths.
setLimitsForConsecutiveViews <- function(x, width)
{
    lw <- length(width)
    first <- integer(lw)
    last <- integer(lw)
    one <- as.integer(1)
    first[1] <- one
    last[1] <- width[1]
    if (lw >= 2) {
        for (i in 2:lw) {
            first[i] <- last[i-1] + one
            last[i] <- first[i] + width[i] - one
        }
    }
    x@first <- first
    x@last <- last
    x
}

# Works if 'src' is a FASTA file:
#   src <- file('Genomes/Celegans/chrI.fa')
#   chrI <- DNAStringViews(src)
DNAStringViews <- function(src)
{
    desc <- character(0)
    if (any(class(src) == "file")) {
        fasta <- readFASTA(src)
        src <- sapply(fasta, function(rec) toupper(rec$seq))
        desc <- sapply(fasta, function(rec) rec$desc)
    } else if (class(src) == "BStringViews") {
        if (class(subject(src)) != "RNAString")
            stop("subject of 'src' is not a \"RNAString\" object")
        x <- src
        x@subject <- DNAString(subject(src))
        return(x)
    }
    big <- paste(src, collapse="")
    x <- new("BStringViews", DNAString(big))
    x@desc <- desc
    setLimitsForConsecutiveViews(x, nchar(src))
}

RNAStringViews <- function(src)
{
    if (class(src) == "BStringViews") {
        if (class(subject(src)) != "DNAString")
            stop("subject of 'src' is not a \"DNAString\" object")
        x <- src
        x@subject <- RNAString(subject(src))
        return(x)
    }
    big <- paste(src, collapse="")
    x <- new("BStringViews", RNAString(big))
    setLimitsForConsecutiveViews(x, nchar(src))
}
