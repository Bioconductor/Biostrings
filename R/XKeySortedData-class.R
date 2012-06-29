### =========================================================================
### XKeySortedData objects
### -------------------------------------------------------------------------


setClass("XKeySortedData",
    contains="DataTable",
    representation(
        "VIRTUAL",
        key="XStringSet",
        table="DataFrame"
    )
)

setClass("BKeySortedData",
    contains="XKeySortedData",
    representation(
        key="BStringSet"
    )
)

setClass("DNAKeySortedData",
    contains="XKeySortedData",
    representation(
        key="DNAStringSet"
    )
)

setClass("RNAKeySortedData",
    contains="XKeySortedData",
    representation(
        key="RNAStringSet"
    )
)

setClass("AAKeySortedData",
    contains="XKeySortedData",
    representation(
        key="AAStringSet"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.XKeySortedData <- function(object)
{
    .Defunct(msg="the XKeySortedData class and subclasses are defunct")
}

setValidity("XKeySortedData",
    function(object)
    {
        problems <- .valid.XKeySortedData(object)
        if (is.null(problems)) TRUE else problems
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setGeneric("dataKey", function(x) standardGeneric("dataKey"))

setMethod("dataKey", "XKeySortedData",
    function(x)
    {
        .Defunct(msg="the dataKey() generic is defunct")
        x@key
    }
)

setGeneric("dataTable", function(x) standardGeneric("dataTable"))

setGeneric("dataTable<-", signature="x",
    function(x, value) standardGeneric("dataTable<-")
)

setMethod("dataTable", "XKeySortedData",
    function(x)
    {
        .Defunct(msg="the dataTable() generic is defunct")
        x@table
    }
)

setReplaceMethod("dataTable", "XKeySortedData",
    function(x, value)
    {
        .Defunct(msg="the dataTable() generic is defunct")
        slot(x, "table") <- value
        x
    }
)

setMethod("length", "XKeySortedData", function(x) length(dataTable(x)))

setMethod("dim", "XKeySortedData", function(x) dim(dataTable(x)))

setMethod("dimnames", "XKeySortedData",
    function(x) list(as.character(dataKey(x)), colnames(dataTable(x)))
)

setMethod("seqtype", "XKeySortedData", function(x) seqtype(dataKey(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### User-friendly versatile constructors.
###

XKeySortedData <- function(seqtype, key=character(0), table=NULL)
{
    .Defunct(msg="the XKeySortedData class and subclasses are defunct")
}

BKeySortedData <- function(key=character(0), table=NULL)
    XKeySortedData("B", key, table)

DNAKeySortedData <- function(key=character(0), table=NULL)
    XKeySortedData("DNA", key, table)

RNAKeySortedData <- function(key=character(0), table=NULL)
    XKeySortedData("RNA", key, table)

AAKeySortedData <- function(key=character(0), table=NULL)
    XKeySortedData("AA", key, table)
