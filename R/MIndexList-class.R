### =========================================================================
### MIndexList objects
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "MIndexList" class.
###
### This class serves as a base class to hold collections of MIndex objects.
###
### As with MIndex, in normal operations, the user should never need to create
### MIndex objects directly or to modify existing ones. Those objects are
### typically returned by a sequence matching/alignment function like
### vmatchPDict(). For this reason, coercion methods are relatively barebones
### and will throw errors rather than try to correct
###

setClass("MIndexList",
    contains="CompressedList",
    representation(
#        "VIRTUAL",
        unlistData="MIndex"
    ),
    prototype(
        elementType="MIndex"
    )
)


.from_list_to_MIndexList <- function(from)
{
    ## keeping ans_class variable for potential future extensions
    ans_class <- "MIndexList"

    ## ensure elements are all MIndex objects
    if (!all(vapply(from, inherits, logical(1L), what="MIndex")))
        stop("Only lists of MIndex objects can be coerced to MIndexList")
    IRanges:::new_CompressedList_from_list(ans_class, from)
}

## barebones, but users shouldn't be calling this directly
## expecting to be used internally, list coercion provided for QoL
setAs("list", "MIndexList",
    function(from) .from_list_to_MIndexList(from)
)
setAs("List", "MIndexList",
    function(from) .from_list_to_MIndexList(as.list(from))
)

.from_MIndexList_to_compact_char <- function(object)
{
    lapply(object,
        function(m_ind) {
            paste(length(m_ind), "patterns,",
                  sum(lengths(m_ind)), "matches")
        }
    )
}

setMethod("show", "MIndexList",
    function(object)
    {
        lo <- length(object)
        k <- min(5, length(object))
        diffK <- lo - 5
        cat(classNameForDisplay(object), " of length ", lo, "\n", sep="")

        repr <- .from_MIndexList_to_compact_char(head(object, k))
        IRanges:::.showAtomicList(CharacterList(repr), minLines=10L)
        if (diffK > 0) {
            cat("...\n<", diffK,
                ifelse(diffK == 1,
                       " more element>\n", " more elements>\n"),
                sep="")
        }
    }
)

setMethod("showAsCell", "MIndexList",
    function(object)
    {
        repr <- .from_MIndexList_to_compact_char(object)
        showAsCell(CharacterList(repr))
    }
)
