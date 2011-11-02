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

setValidity("WCP",
    function(object)
        .Defunct(msg="the WCP class and subclasses are defunct")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

setMethod("show", "WCP",
    function(object)
        .Defunct(msg="the WCP class and subclasses are defunct")
)
