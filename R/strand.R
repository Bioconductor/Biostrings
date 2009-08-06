### The strand generic and the strand() method, for lack of a better home.

setGeneric("strand", function(x) standardGeneric("strand"))

setMethod("strand", "missing", function(x) factor(levels=c("-","+","*")))

setMethod("strand", "character",
    function(x) {
        lvls <- levels(strand())
        if (!all(is.na(x) | (x %in% lvls)))
            stop("strand values must be in '", paste(lvls, collapse="' '"), "'")
        factor(x, levels=lvls)
    })

setMethod("strand", "DataTable",
    function(x) {
        ans <- x[["strand"]]
        if (is.null(ans))
            ans <- rep(NA_character_, nrow(x))
        if (is.character(ans))
            ans <- strand(ans)
        ans
    })
