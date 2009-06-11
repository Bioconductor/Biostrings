### =========================================================================
### WCP objects
### -------------------------------------------------------------------------


setClass("WCP",
    representation(
        "VIRTUAL",
        clusters="Binning"
    )
)

setClass("B_WCP",
    contains="WCP",
    representation(
        dictList="BKeySortedDataList"
    )
)

setClass("DNA_WCP",
    contains="WCP",
    representation(
        dictList="DNAKeySortedDataList"
    )
)

setClass("RNA_WCP",
    contains="WCP",
    representation(
        dictList="RNAKeySortedDataList"
    )
)

setClass("AA_WCP",
    contains="WCP",
    representation(
        dictList="AAKeySortedDataList"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.WCP <- function(object)
{
    ncharStrings <- lapply(object@dictList, function(x) unique(nchar(dataKey(x))))
    if (!all(sapply(ncharStrings, length) == 1)) {
        message <- "within cluster strings do not have the same nchar"
    } else {
        ncharStrings <- unlist(ncharStrings)
        binSizes <- elementLengths(object@clusters)
        if (length(ncharStrings) != length(binSizes)) {
            message <- "number of clusters do not match the number expected"
        } else if (any(ncharStrings != binSizes)) {
            message <- "within cluster nchar do not match the number expected"
        } else {
            message <- NULL
        }
    }
    message
}

setValidity("WCP",
    function(object)
    {
        problems <- .valid.WCP(object)
        if (is.null(problems)) TRUE else problems
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "xsbasetype" method.
###

setMethod("xsbasetype", "WCP", function(x) xsbasetype(x@dictList))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

setMethod("show", "WCP",
    function(object)
    {
        NG <- length(object@clusters)
        NO <- nobj(object@clusters)
        cat("  A ", class(object), " instance with ", NG,
            ifelse(NG == 1, " cluster on ", " clusters on "), NO,
            ifelse(NO == 1, " letter\n", " letters\n"), sep = "")
    }
)
