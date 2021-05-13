test_replace_seqinfo <- function() {
    value <- Seqinfo(c("a", "b", "c"), c(9, 2, 0))
    x <- DNAStringSet(c(a="TTTCATAGG", b="AA"))
    ## names are not identical (length)
    checkException({
        seqinfo(x) <- value
    })

    x <- DNAStringSet(c(A="TTTCATAGG", B="AA", C="TCGA"))
    value <- Seqinfo(c("a", "b", "c"), c(9, 2, 0))
    ## names are not identical (case)
    checkException({
        seqinfo(x) <- value
    })

    x <- DNAStringSet(c(a="TTTCATAGG", b="AA", c="TCGA"))
    value <- Seqinfo(c("a", "b", "c"), c(9, 2, 0))
    ## widths do not match
    checkException({
        seqinfo(x) <- value
    })

    x <- DNAStringSet(c("TTTCATAGG", "AA",  "TCGA"))
    value <- Seqinfo(c("1", "2", "3"), c(9, 2, 0))
    ## names match but not widths
    checkException({
        seqinfo(x) <- value
    })

    x <- DNAStringSet(c("TTTCATAGG", "AA",  "TCGA"))
    value <- Seqinfo(c("1", "2", "3"), c(9, 2, 4))
    ## names and widths match
    seqinfo(x) <- value
    checkTrue(validObject(x))
}
