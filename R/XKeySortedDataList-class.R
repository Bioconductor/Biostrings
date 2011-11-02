### =========================================================================
### XKeySortedDataList objects
### -------------------------------------------------------------------------

setClass("XKeySortedDataList",
    contains = "SimpleList",
    representation = representation("VIRTUAL"),
    prototype = prototype(elementType = "XKeySortedData")
)

setClass("BKeySortedDataList",
    contains = "XKeySortedDataList",
    prototype = prototype(elementType = "BKeySortedData")
)

setClass("DNAKeySortedDataList",
    contains = "XKeySortedDataList",
    prototype = prototype(elementType = "DNAKeySortedData")
)

setClass("RNAKeySortedDataList",
    contains = "XKeySortedDataList",
    prototype = prototype(elementType = "RNAKeySortedData")
)

setClass("AAKeySortedDataList",
    contains = "XKeySortedDataList",
    prototype = prototype(elementType = "AAKeySortedData")
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "xsbasetype" method.
###

setMethod("xsbasetype", "XKeySortedDataList",
    function(x) xsbasetype(new(x@elementType))
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### User-friendly versatile constructors.
###

BKeySortedDataList <- function(...)
{
  .Deprecated(msg="the XKeySortedDataList class and subclasses are deprecated")
  IRanges:::newSimpleList("BKeySortedDataList", list(...))
}

DNAKeySortedDataList <- function(...)
{
  .Deprecated(msg="the XKeySortedDataList class and subclasses are deprecated")
  IRanges:::newSimpleList("DNAKeySortedDataList", list(...))
}

RNAKeySortedDataList <- function(...)
{
  .Deprecated(msg="the XKeySortedDataList class and subclasses are deprecated")
  IRanges:::newSimpleList("RNAKeySortedDataList", list(...))
}
AAKeySortedDataList <- function(...)
{
  .Deprecated(msg="the XKeySortedDataList class and subclasses are deprecated")
  IRanges:::newSimpleList("AAKeySortedDataList", list(...))
}

