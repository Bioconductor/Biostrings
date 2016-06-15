### =========================================================================
### Preprocessed Subject Strings of type BOC  -- DEFUNCT!  DEFUNCT!  DEFUNCT!
### -------------------------------------------------------------------------


### Note that 'base1_code', 'base2_code' and 'base3_code' must be distinct
setClass("BOC_SubjectString",
    representation(
        subject="DNAString",        # TODO: support "RNAString" too
        pattern_length="integer",   # A single integer e.g. 36
        base1_code="integer",       # integer belonging to DNA_BASE_CODES
        base2_code="integer",
        base3_code="integer",
        base4_code="integer",
        base1_OCbuffer="SharedRaw",      # all buffers must be of length nchar(subject) - pattern_length + 1
        base2_OCbuffer="SharedRaw",
        base3_OCbuffer="SharedRaw",
	pre4buffer="SharedRaw",
        ## The "stats" slot is a named list with the following elements:
        ##   means: vector of 4 doubles
        ##   table1, table2, table3, table4: vectors of (pattern_length + 1) integers
        stats="list"
    )
)

setMethod("initialize", "BOC_SubjectString",
    function(.Object, subject, pattern_length, base_letters)
        .Defunct(msg="BOC_SubjectString objects are defunct")
)

### Typical use:
###   Biostrings:::plotBOC(chr1boc, "Human chr1")
plotBOC <- function(x, main)
    .Defunct(msg="BOC_SubjectString objects are defunct")

### Dispatch on 'subject' (see signature of generic).
### 'algorithm' is ignored.
setMethod("matchPattern", "BOC_SubjectString",
    function(pattern, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto")
    {
        .Defunct(msg="BOC_SubjectString objects are defunct")
    }
)

