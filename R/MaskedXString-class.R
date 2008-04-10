### =========================================================================
### MaskedXString objects
### -------------------------------------------------------------------------

setClass("MaskedXString",
    contains="XString",
    representation(
        "VIRTUAL",
        mask="NormalIRanges"
    )
)

### 4 direct "MaskedXString" derivations (no additional slot)
setClass("MaskedBString", contains=c("MaskedXString", "BString"))
setClass("MaskedDNAString", contains=c("MaskedXString", "DNAString"))
setClass("MaskedRNAString", contains=c("MaskedXString", "RNAString"))
setClass("MaskedAAString", contains=c("MaskedXString", "AAString"))

