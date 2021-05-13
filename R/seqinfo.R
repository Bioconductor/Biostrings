
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

setReplaceMethod("seqinfo", "DNAStringSet",
    function(x, new2old = NULL, pruning.mode = c("error", "coarse", "fine", "tidy"), value) {
        if (!is.null(new2old))
            stop("'new2old' not supported in 'seqinfo<-,DNAStringSet-method'",
                call. = FALSE)

        xnames <- names(x)
        same_names <-
            if (is.null(xnames))
                identical(as.character(seq_along(x)), names(value))
            else
                identical(xnames, names(value))
        if (!same_names)
            stop("'names(x)' and 'names(value)' are not identical")

        if (!identical(width(x), unname(seqlengths(value))))
            stop("'width(x)' is not the same as 'seqlengths(value)'")

        metadata(x)$seqinfo <- value
        x
    }
)


