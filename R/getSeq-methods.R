setMethod("getSeq", "XStringSet",
    function(x, names)
    {
        stopifnot(is.character(names) || is(names, "GRanges") ||
                  is(names, "GRangesList"),
                  !is.null(names(x)))
        if (is.character(names)) {
            found <- names %in% names(x)
            regexNames <- unlist(lapply(names[!found], grep, names(x),
                                        value=TRUE))
            names <- c(names[found], regexNames)
            return(x[names])
        } else if (is(names, "GRangesList")) {
            gr <- unlist(names, use.names=FALSE)
        } else {
            gr <- names
        }
        ignoringStrand <- any(strand(gr) != "*") &&
            !hasMethod(reverseComplement, class(x))
        if (ignoringStrand) {
            warning("some strand(x) != '*' but ",
                    "strand has no meaning for ", class(x))
        }
        rl <- as(gr, "IntegerRangesList")
        ans <- unsplit(extractAt(x[names(rl)], unname(rl)), seqnames(gr))
        if (!ignoringStrand) {
            minus <- strand(gr) == "-"
            ans[minus] <- reverseComplement(ans[minus])
        }
        if (is(names, "GRangesList")) {
            ans <- relist(ans, names)
        }
        ans
    }
)

