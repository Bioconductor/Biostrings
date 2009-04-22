### =========================================================================
### XStringSetList objects
### -------------------------------------------------------------------------
###

setClass("XStringSetList",
    contains="ListLike",
    representation(
        "VIRTUAL",
        unlisted="XStringSet",
        cum_eltlength="integer"
    )
)

setClass("BStringSetList",
    contains="XStringSetList",
    representation(
        unlisted="BStringSet"
    )
)
setClass("DNAStringSetList",
    contains="XStringSetList",
    representation(
        unlisted="DNAStringSet"
    )
)
setClass("RNAStringSetList",
    contains="XStringSetList",
    representation(
        unlisted="RNAStringSet"
    )
)
setClass("AAStringSetList",
    contains="XStringSetList",
    representation(
        unlisted="AAStringSet"
    )
)

