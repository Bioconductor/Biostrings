### =========================================================================
### The matchWCP() generic & related functions
### -------------------------------------------------------------------------

setMethod("maxWeights", "WCP",
    function(x)
        .Defunct(msg="\"maxWeights\" method for WCP objects is defunct")
)

WCPscoreStartingAt <- function(wcp, subject, starting.at=1)
{
    .Defunct(msg="WCPscoreStartingAt() is defunct")
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "matchWCP" generic and methods.
###

setGeneric("matchWCP", signature="subject",
    function(wcp, subject, min.score="80%")
        standardGeneric("matchWCP")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchWCP", "character",
    function(wcp, subject, min.score="80%")
        matchWCP(wcp, XString(xsbasetype(wcp), subject), min.score)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchWCP", "XString",
    function(wcp, subject, min.score="80%")
        .Defunct(msg="matchWCP() is defunct")
)

### Dispatch on 'subject' (see signature of generic).
### WARNING: Unlike the other "matchWCP" methods, the XStringViews object
### returned by this method is not guaranteed to have its views ordered from
### left to right in general! One important particular case where this is
### guaranteed though is when 'isNormal(subject)' is TRUE (i.e. 'subject' is
### a normal XStringViews object).
### matchWCP does not support "out of limits"  matches.
setMethod("matchWCP", "XStringViews",
    function(wcp, subject, min.score="80%")
        .Defunct(msg="matchWCP() is defunct")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchWCP", "MaskedXString",
    function(wcp, subject, min.score="80%")
        .Defunct(msg="matchWCP() is defunct")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "countWCP" generic and methods.
###

setGeneric("countWCP", signature="subject",
    function(wcp, subject, min.score="80%")
        standardGeneric("countWCP")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countWCP", "character",
    function(wcp, subject, min.score="80%")
        countWCP(wcp, XString(xsbasetype(wcp), subject), min.score)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countWCP", "XString",
    function(wcp, subject, min.score="80%")
        .Defunct(msg="countWCP() is defunct")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countWCP", "XStringViews",
    function(wcp, subject, min.score="80%")
        .Defunct(msg="countWCP() is defunct")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countWCP", "MaskedXString",
    function(wcp, subject, min.score="80%")
        .Defunct(msg="countWCP() is defunct")
)

