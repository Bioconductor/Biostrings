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
### The "seqtype" method.
###

setMethod("seqtype", "XKeySortedDataList",
    function(x) seqtype(new(x@elementType))
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### User-friendly versatile constructors.
###

BKeySortedDataList <- function(...)
{
  .Defunct(msg="the XKeySortedDataList class and subclasses are defunct")
}

DNAKeySortedDataList <- function(...)
{
  .Defunct(msg="the XKeySortedDataList class and subclasses are defunct")
}

RNAKeySortedDataList <- function(...)
{
  .Defunct(msg="the XKeySortedDataList class and subclasses are defunct")
}
AAKeySortedDataList <- function(...)
{
  .Defunct(msg="the XKeySortedDataList class and subclasses are defunct")
}

