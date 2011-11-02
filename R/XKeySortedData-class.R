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
    .Deprecated(msg="the XKeySortedData class and subclasses are deprecated")
    message <- NULL
    if (!is.null(elementMetadata(object)))
        message <- c(message, "elementMetdata is not NULL")
    if (length(dataKey(object)) != nrow(dataTable(object)))
        message <-
          c(message, "length(dataKey(object)) != nrow(dataTable(object))")
    if (is.unsorted(dataKey(object), strictly = TRUE)) {
        if (is.unsorted(sort(dataKey(object)), strictly = TRUE))
            message <- c(message, "dataKey(object) contains duplicates")
        else
            message <- c(message, "dataKey(object) is unsorted")
    }
    message
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
        .Deprecated(msg="the dataKey() generic is deprecated")
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
        .Deprecated(msg="the dataTable() generic is deprecated")
        x@table
    }
)

setReplaceMethod("dataTable", "XKeySortedData",
    function(x, value)
    {
        .Deprecated(msg="the dataTable() generic is deprecated")
        slot(x, "table") <- value
        x
    }
)

setMethod("length", "XKeySortedData", function(x) length(dataTable(x)))

setMethod("dim", "XKeySortedData", function(x) dim(dataTable(x)))

setMethod("dimnames", "XKeySortedData",
    function(x) list(as.character(dataKey(x)), colnames(dataTable(x)))
)

setMethod("xsbasetype", "XKeySortedData", function(x) xsbasetype(dataKey(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

setMethod("[[", "XKeySortedData",
    function(x, i, j, ...) dataTable(x)[[i]]
)

setMethod("[", "XKeySortedData",
    function(x, i, j, ..., drop)
    {
        if (length(list(...)) > 0)
            warning("parameters in '...' not supported")
        if ((nargs() - !missing(drop)) < 3) {
            if (!missing(drop))
                warning("parameter 'drop' ignored by list-style subsetting")
            if (!missing(i))
                dataTable(x) <- dataTable(x)[i]
        } else {
            if (missing(i)) {
                if (!missing(j)) {
                    if (missing(drop)) {
                        subset <- dataTable(x)[,j]
                    } else {
                        subset <- dataTable(x)[,j,drop=drop]
                    }
                    if (is(subset, "DataFrame")) {
                        dataTable(x) <- subset
                    } else {
                        x <- subset
                    }
                }
            } else {
                if (is.character(i))
                    i <- XStringSet(xsbasetype(dataKey(x)), i)
                if (is(i, "XStringSet")) {
                    i <- match(i, dataKey(x))
                    if (any(is.na(i)))
                        stop("selecting rows: subscript out of bounds")
                }
                if (is.integer(i) && is.unsorted(i, strictly = TRUE))
                    stop("cannot extract duplicate items or items out of order")
                if (missing(drop)) {
                    subset <- dataTable(x)[i,j]
                } else {
                    subset <- dataTable(x)[i,j,drop=drop]
                }
                if (is(subset, "DataFrame")) {
                    x <- x[i]
                    dataTable(x) <- subset
                } else {
                    x <- subset
                }
            } 
        }
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### User-friendly versatile constructors.
###

XKeySortedData <- function(basetype, key=character(0), table=NULL)
{
    .Deprecated(msg="the XKeySortedData class and subclasses are deprecated")
    stringSetClass <- paste(basetype, "StringSet", sep = "")
    stringSortedDataClass <- paste(basetype, "KeySortedData", sep = "")
    if (!is(key, stringSetClass))
        key <- do.call(stringSetClass, list(key))
    if (!is(table, "DataFrame"))
        table <- as(table, "DataFrame")
    if (length(key) != nrow(table))
        stop("length(key) != nrow(table)")
    orderKey <- order(key)
    new(stringSortedDataClass,
        key = key[orderKey],
        table = table[orderKey,,drop=FALSE])
}

BKeySortedData <- function(key=character(0), table=NULL)
    XKeySortedData("B", key, table)

DNAKeySortedData <- function(key=character(0), table=NULL)
    XKeySortedData("DNA", key, table)

RNAKeySortedData <- function(key=character(0), table=NULL)
    XKeySortedData("RNA", key, table)

AAKeySortedData <- function(key=character(0), table=NULL)
    XKeySortedData("AA", key, table)
