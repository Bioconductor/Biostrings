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
        ans <- src
        ans@subject <- DNAString(subject(src))
        return(ans)
    }
    big <- paste(src, collapse="")
    ans <- new("BStringViews", DNAString(big))
    ans@desc <- desc
    setLimitsForAdjacentViews(ans, nchar(src))
}

RNAStringViews <- function(src)
{
    if (class(src) == "BStringViews") {
        if (class(subject(src)) != "DNAString")
            stop("subject of 'src' is not a \"DNAString\" object")
        ans <- src
        ans@subject <- RNAString(subject(src))
        return(ans)
    }
    big <- paste(src, collapse="")
    ans <- new("BStringViews", RNAString(big))
    setLimitsForAdjacentViews(ans, nchar(src))
}
