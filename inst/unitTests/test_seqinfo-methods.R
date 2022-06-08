test_seqinfo <- function()
{
    x0 <- DNAStringSet(c("AAA", "TTTTT", "GG"))

    target <- Seqinfo(as.character(1:3), width(x0))
    checkIdentical(target, seqinfo(x0))
    checkIdentical(seqlevels(target), seqlevels(x0))
    checkIdentical(isCircular(target), isCircular(x0))
    checkIdentical(genome(target), genome(x0))

    ## no-op:
    x <- x0
    seqinfo(x) <- seqinfo(x)
    checkIdentical(x0, x)

    ## The seqnames in the supplied Seqinfo will set or replace the
    ## current names on 'x':
    value <- Seqinfo(c("chr1", "chr2", "chrM"))
    checkException(seqinfo(x) <- value)  # seqlengths don't match
    seqlengths(value) <- width(x)
    seqinfo(x) <- value
    checkIdentical(value, seqinfo(x))
    checkIdentical(seqnames(value), names(x))

    ## Set circularity flags and/or genome:
    isCircular(x) <- c(TRUE, TRUE, FALSE)
    target <- setNames(c(TRUE, TRUE, FALSE), seqlevels(x))
    checkIdentical(target, isCircular(x))
    checkIdentical(target, metadata(x)$is_circular)
    genome(x) <- "The Dude"
    target <- setNames(rep("The Dude", 3), seqlevels(x))
    checkIdentical(target, genome(x))
    checkIdentical(target, metadata(x)$genome)

    ## Change seqlevels style:
    seqlevelsStyle(x) <- "Ensembl"
    checkIdentical(c("1", "2", "MT"), names(x))

    ## Remove names, circularity flags, and genome:
    seqinfo(x) <- seqinfo(x0)
    checkIdentical(x0, x)
}

