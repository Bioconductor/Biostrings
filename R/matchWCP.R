### =========================================================================
### The matchWCP() generic & related functions
### -------------------------------------------------------------------------

setMethod("maxWeights", "WCP",
    function(x)
    {
        .Deprecated(msg="\"maxWeights\" method for WCP objects is deprecated")
        sapply(x@dictList, function(y) max(dataTable(y)[[1]]))
    })

WCPscoreStartingAt <- function(wcp, subject, starting.at=1)
{
    .Deprecated(msg="WCPscoreStartingAt() is deprecated")
    if (!is(wcp, "WCP"))
        stop("'wcp' must be a WCP object")
    ## checking 'subject'
    subject <- .normargSubject(subject)
    ## checking 'starting.at'
    if (!is.numeric(starting.at))
        stop("'starting.at' must be a vector of integers")
    if (!is.integer(starting.at))
        starting.at <- as.integer(starting.at)

    .Call2("WCP_score_starting_at", wcp, subject, starting.at, PACKAGE="Biostrings")
}

.XString.matchWCP <- function(wcp, subject, min.score, count.only=FALSE)
{
    if (!is(wcp, "WCP"))
        stop("'wcp' must be a WCP object")
    if (xsbasetype(wcp) != xsbasetype(subject))
        stop("'wcp' and 'subject' must have the same XString base type")
    min.score <- .normargMinScore(min.score, wcp)    
    C_ans <- .Call2("XString_match_WCP",
            wcp, subject, min.score, count.only,
            PACKAGE="Biostrings")
    if (count.only)
        return(C_ans)
    unsafe.newXStringViews(subject, start(C_ans), width(C_ans))
}

.XStringViews.matchWCP <- function(wcp, subject, min.score, count.only=FALSE)
{
    if (!is(wcp, "WCP"))
        stop("'wcp' must be a WCP object")
    if (xsbasetype(wcp) != xsbasetype(subject(subject)))
        stop("'wcp' and 'subject' must have the same XString base type")
    min.score <- .normargMinScore(min.score, wcp)
    C_ans <- .Call2("XStringViews_match_WCP",
            wcp,
            subject(subject), start(subject), width(subject),
            min.score, count.only,
            PACKAGE="Biostrings")
    if (count.only)
        return(C_ans)
    unsafe.newXStringViews(subject(subject), start(C_ans), width(C_ans))
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
    function(wcp, subject, min.score="80%") {
        .Deprecated(msg="matchWCP() is deprecated")
        .XString.matchWCP(wcp, subject, min.score)
    }
)

### Dispatch on 'subject' (see signature of generic).
### WARNING: Unlike the other "matchWCP" methods, the XStringViews object
### returned by this method is not guaranteed to have its views ordered from
### left to right in general! One important particular case where this is
### guaranteed though is when 'isNormal(subject)' is TRUE (i.e. 'subject' is
### a normal XStringViews object).
### matchWCP does not support "out of limits"  matches.
setMethod("matchWCP", "XStringViews",
    function(wcp, subject, min.score="80%") {
        .Deprecated(msg="matchWCP() is deprecated")
        .XStringViews.matchWCP(wcp, subject, min.score)
    }
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchWCP", "MaskedXString",
    function(wcp, subject, min.score="80%") {
        .Deprecated(msg="matchWCP() is deprecated")
        matchWCP(wcp, toXStringViewsOrXString(subject), min.score)
    }
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
    function(wcp, subject, min.score="80%") {
        .Deprecated(msg="countWCP() is deprecated")
        .XString.matchWCP(wcp, subject, min.score, count.only=TRUE)
    }
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countWCP", "XStringViews",
    function(wcp, subject, min.score="80%") {
        .Deprecated(msg="countWCP() is deprecated")
        .XStringViews.matchWCP(wcp, subject, min.score, count.only=TRUE)
    }
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countWCP", "MaskedXString",
    function(wcp, subject, min.score="80%") {
        .Deprecated(msg="countWCP() is deprecated")
        countWCP(wcp, toXStringViewsOrXString(subject), min.score)
    }
)
