
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqinfo setter and getter
###

setMethod("seqinfo", "DNAStringSet", function(x) {
    si <- metadata(x)$seqinfo
    if (is.null(si)) {
        sn <- names(x)
        if (is.null(sn))
            sn <- as.character(seq_along(x))
        si <- Seqinfo(sn, width(x))
    }
    si
})

setReplaceMethod("seqinfo", "DNAStringSet", function(x, value) {
    metadata(x)$seqinfo <- value
    x
})


