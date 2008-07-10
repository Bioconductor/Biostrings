### ==========================================================================
### InDel objects
### --------------------------------------------------------------------------
### An InDel object contains the insertion and deletion information.


setClass("InDel",
    representation(
        insertion="ANY",
        deletion="ANY"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setGeneric("insertion", function(x) standardGeneric("insertion"))
setMethod("insertion", "InDel", function(x) x@insertion)

setGeneric("deletion", function(x) standardGeneric("deletion"))
setMethod("deletion", "InDel", function(x) x@deletion)
