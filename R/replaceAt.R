### =========================================================================
### replaceAt() & insertAt()
### -------------------------------------------------------------------------


### Only checks that the ranges in 'at' are within the bounds specified by
### 'min_start' and 'max_end'. Returns the nb of ranges in 'at'.
.checkargAt <- function(at, min_start, max_end)
{
    if (!is(at, "Ranges"))
        stop("'at' must be a Ranges object")
    at_len <- length(at)
    if (at_len == 0L)
        return(TRUE)
    at_min_start <- min(start(at))
    at_max_end <- max(end(at))
    if (at_min_start < min_start || at_max_end > max_end)
        stop("some ranges in 'at' are off-limits with respect to sequence 'x'")
    at_len
}

.normargAt <- function(at, x)
{
    if (!is(at, "RangesList"))
        stop("'at' must be a RangesList object")
    if (!is.null(names(at))) {
        names(at) <- NULL
        warning("'at' names were ignored")
    }
    at_len <- length(at)
    x_len <- length(x)
    if (at_len > x_len)
        stop("'at' cannot be longer than 'x'")
    if (at_len < x_len) {
        if (at_len == 0L)
            stop("'at' is a zero-length object but not 'x'")
        if (x_len %% at_len != 0L)
            warning("'length(x)' is not a multiple of 'length(at)'")
        at <- rep(at, length.out=x_len)
    }
    ## range() is fast only on a CompressedIRangesList.
    at_range <- range(at)
    idx <- which(elementLengths(at_range) != 0L)
    if (length(idx) != 0L) {
        tmp <- unlist(at_range)
        if (any(start(tmp) < 1L) || any(end(tmp) > width(x)[idx]))
            stop("some ranges in 'at' are off-limits with respect to ",
                 "their corresponding\n  sequence in 'x'")
    }
    at
}

.normargValue <- function(value, x_seqtype, NR)
{
    if (!is(value, "XStringSet")) {
        value_class <- paste0(x_seqtype, "StringSet")
        value <- as(value, value_class)
    } else if (seqtype(value) != x_seqtype) {
        seqtype(value) <- x_seqtype
    }
    value_len <- length(value)
    if (value_len > NR)
        stop("'value' cannot be longer than the number of replacements")
    if (value_len < NR) {
        if (value_len == 0L)
            stop("'length(value)' is zero but not the number ",
                 "of replacements")
        if (NR %% value_len != 0L)
            warning("the number of replacements is not ",
                    "a multiple of 'length(value)'")
        value <- rep(value, length.out=NR)
    }
    value
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### replaceAt()
###
### Performs replacements (a.k.a. substitutions) in an XString object 'x',
### or in each sequence of an XStringSet object 'x', at the position ranges
### specified in 'at'.
###

setGeneric("replaceAt", signature="x",
    function(x, at, value="") standardGeneric("replaceAt")
)

### 'at': a Ranges object containing the locations of the replacements.
###       The number of replacements (NR) is the length of 'at'.
### 'value': character vector or XStringSet object of length NR (or recycled
###       to length NR if necessary).
setMethod("replaceAt", "XString",
    function(x, at, value="")
    {
        x_len <- length(x)
        NR <- .checkargAt(at, 1L, x_len)  # nb of replacements

        ## Normalize 'value'.
        value <- .normargValue(value, seqtype(x), NR)
        if (NR == 0L)
            return(x)

        oo <- order(at)
        start2 <- start(at)[oo]
        end2 <- end(at)[oo]
        if (any(start2[-1L] <= end2[-NR]))
            stop("'at' must contain disjoint ranges (see '?isDisjoint')")
        value2 <- value[oo]
        ## The regions of 'x' that are preserved are the regions between the
        ## ranges in 'at'.
        preserved <- IRanges(c(1L, end2 + 1L), c(start2 - 1L, x_len))
        x_preserved <- as(Views(x, preserved), "XStringSet")
        last_value2 <- as("", class(value2))
        unlist(xscat(x_preserved, c(value2, last_value2)))
    }
)

### 'at': a RangesList object of the length of 'x' (recycled to the length
###       of 'x' if necessary) containing the locations of the replacements
###       for each sequence in 'x'.
###       The total number of replacements (NR) is the total number of ranges
###       in 'at' i.e. 'sum(elementLengths(at))' or 'length(unlist(at))'
###       after recycling of 'at'.
### 'value': character vector or XStringSet object of length NR (or recycled
###       to length NR if necessary).
### TODO: Support CharacterList or XStringSetList for 'value'. In that case
###       recycling should happen vertically and horizontally to bring it in
###       the same shape as 'at' (after recycling).
setMethod("replaceAt", "XStringSet",
    function(x, at, value="")
    {
        x_len <- length(x)

        ## Normalize 'at'.
        at <- .normargAt(at, x)
        unlisted_at <- unlist(at, use.names=FALSE)
        NR <- length(unlisted_at)  # nb of replacements

        ## Normalize 'value'.
        value <- .normargValue(value, seqtype(x), NR)
        if (NR == 0L)
            return(x)

        x_width <- width(x)
        at_eltlens <- elementLengths(at)

        unlisted_x <- unlist(x, use.names=FALSE)
        offsets <- cumsum(c(0L, x_width[-x_len]))
        offsets <- rep.int(offsets, at_eltlens)

        unlisted_ans <- replaceAt(unlisted_x,
                                  shift(unlisted_at, shift=offsets),
                                  value=value)

        ans_width <- x_width + sum(relist(width(value) - width(unlisted_at),
                                          PartitioningByEnd(at)))
        as(successiveViews(unlisted_ans, ans_width), "XStringSet") 
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### insertAt()
###
### Performs "multiple inserts" in an XString object 'x', or in each sequence
### of an XStringSet object 'x', *before* the positions passed in 'at'.
### TODO: Re-implement as a simple wrapper to insertAt().
###

setGeneric("insertAt", signature="x",
    function(x, at, value="") standardGeneric("insertAt")
)

### 'at': an integer vector containing the locations of the insertions. Each
###       insertion will happen right *before* the specified location, that is,
###       between the specified location and the location minus one.
###       The number of insertions NI is the length of 'at'.
### 'value': character vector or XStringSet object of length NI (or recycled
###       to NI if needed).
setMethod("insertAt", "XString",
    function(x, at, value="")
    {
        x_len <- length(x)

        ## Normalize 'at'.
        if (!is.numeric(at))
            stop("'at' must be a vector of integers")
        if (!is.integer(at))
            at <- as.integer(at)
        if (IRanges:::anyMissingOrOutside(at, 1L, x_len + 1L))
            stop("'at' contains invalid positions in 'x'")
        NI <- length(at)  # nb of insertions

        ## Normalize 'value'.
        value <- .normargValue(value, seqtype(x), NI)

        if (NI == 0L)
            return(x)

        oo <- IRanges:::orderInteger(at)
        at2 <- at[oo]
        value2 <- value[oo]
        last_width <- x_len - at2[NI] + 1L
        last_value2 <- as("", class(value2))
        chunks_width <- c(diff(c(1L, at2)), last_width)
        chunks <- as(successiveViews(x, chunks_width), "XStringSet")
        unlist(xscat(chunks, c(value2, last_value2)))
    }
)

### 'at': a list of integer vectors (or IntegerList object) of the length
###       of 'x' (recycled to the length of 'x' if shorter than 'x') containing
###       the locations of the insertions for each sequence in 'x'. Each
###       insertion will happen right *before* the specified location, that is,
###       between the specified location and the location minus one.
###       The total number of insertions NI is defined by 'length(unlist(at))'
###       after recycling of 'at'.
### 'value': character vector or XStringSet object of length NI (or recycled
###       to NI if needed).
### TODO: Like for replaceAt(), support CharacterList or XStringSetList for
###       'value'. In that case recycling should happen vertically and
###       horizontally to bring it in the same shape as 'at' (after recycling).
setMethod("insertAt", "XStringSet",
    function(x, at, value="")
    {
        x_len <- length(x)

        ## Normalize 'at'.
        if (!is(at, "IntegerList"))
            at <- as(at, "IntegerList")
        if (!is.null(names(at))) {
            names(at) <- NULL
            warning("'at' names were ignored")
        }
        at_len <- length(at)
        if (at_len > x_len)
            stop("'at' cannot have more list elements than 'x'")
        if (at_len < x_len) {
            if (at_len == 0L)
                stop("'at' is a zero-length list but not 'x'")
            if (x_len %% at_len != 0L)
                warning("'length(x)' is not a multiple of 'length(at)'")
            at <- rep(at, length.out=x_len)
        }
        unlisted_at <- unlist(at, use.names=FALSE)
        NI <- length(unlisted_at)  # nb of insertions

        ## Normalize 'value'.
        value <- .normargValue(value, seqtype(x), NI)

        if (NI == 0L)
            return(x)

        x_width <- width(x)
        at_eltlens <- elementLengths(at)
        max_at <- rep.int(x_width + 1L, at_eltlens)
        if (any(is.na(unlisted_at))
         || any(unlisted_at < 1L)
         || any(unlisted_at > max_at))
            stop("'at' contains invalid positions in 'x'")

        unlisted_x <- unlist(x, use.names=FALSE)
        offsets <- cumsum(c(0L, x_width[-x_len]))
        offsets <- rep.int(offsets, at_eltlens)
        unlisted_ans <- insertAt(unlisted_x, unlisted_at + offsets, value)
        ans_width <- x_width + sum(relist(width(value),
                                          PartitioningByEnd(at)))
        as(successiveViews(unlisted_ans, ans_width), "XStringSet") 
    }
)


if (FALSE) {  #      <<<--- begin testing block --->>>
### -------------------------------------------------------------------------
### Testing replaceAt()

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### On an XString object

x <- BString("abcdefghijklm")
at1 <- IRanges(c(1, 5, 12), c(3, 4, 12))
value <- c("X", "Y", "Z")
stopifnot(replaceAt(x, at1, value) == "XdYefghijkZm")
stopifnot(replaceAt(x, rev(at1), rev(value)) == "XdYefghijkZm")

at2 <- IRanges(c(14, 1, 1, 1, 1, 11), c(13, 0, 10, 0, 0, 10))
value <- 1:6
stopifnot(replaceAt(x, at2, value) == "24536klm1")
stopifnot(replaceAt(x, rev(at2), rev(value)) == "54236klm1")

### Deletions:
stopifnot(replaceAt(x, at1) == "defghijkm")
stopifnot(replaceAt(x, rev(at1)) == "defghijkm")
stopifnot(replaceAt(x, at2) == "klm")
stopifnot(replaceAt(x, rev(at2)) == "klm")

### Insertions:
at3 <- IRanges(c(6L, 10L, 2L, 5L), width=0L)
stopifnot(replaceAt(x, at3, value="-") == "a-bcd-e-fghi-jklm")

### No-ops:
stopifnot(replaceAt(x, at1, Views(x, at1)) == x)
stopifnot(replaceAt(x, at2, Views(x, at2)) == x)
stopifnot(replaceAt(x, at3, Views(x, at3)) == x)

### Testing performance with half-million small substitutions at random
### locations in Celegans chrI:
library(BSgenome.Celegans.UCSC.ce2)
genome <- BSgenome.Celegans.UCSC.ce2
chrI <- genome$chrI

### A very amateurish random disjoint ranges generator. Not for serious use!
randomDisjointRanges <- function(min_start, max_end, min_width, max_width, n)
{
    set.seed(33)
    offset <- min_start - 1L
    n1 <- n * (1.1 + max_width/(max_end-offset+1L))  # very approximative

    some_starts <- sample(max_end-offset+1L, n1, replace=TRUE) + offset
    some_widths <- sample(min_width:max_width, n1, replace=TRUE)
    some_ranges <- IRanges(some_starts, width=some_widths)
    some_ranges <- some_ranges[end(some_ranges) <= max_end]
    ans <- disjoin(some_ranges)
    if (min_width == 0L) {
        is_empty <- width(some_ranges) == 0L
        some_empty_ranges <- some_ranges[is_empty]
        ans <- sample(c(ans, some_empty_ranges))
    }
    if (length(ans) < n)
        stop("internal error, sorry")
    head(ans, n=n)
}

at4 <- randomDisjointRanges(1L, nchar(chrI), 0L, 20L, 500000L)
### Takes < 1s on my machine (64-bit Ubuntu).
system.time(current <- replaceAt(chrI, at4, Views(chrI, at4)))
stopifnot(current == chrI)

### Testing performance with half-million single-letter insertions at random
### locations in Celegans chrI:
at5 <- randomDisjointRanges(1L, nchar(chrI), 0L, 0L, 500000L)
### Takes < 1s on my machine (64-bit Ubuntu).
system.time(current <- replaceAt(chrI, at5, value="-"))
m <- matchPattern("-", current)
stopifnot(identical(sort(start(at5)), start(m) - seq_along(at5) + 1L))

system.time(current2 <- insertAt(chrI, start(at5), value="-"))
stopifnot(identical(current, current2))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### On an XStringSet object

x <- BStringSet(c("abcdefghijklmn", "ABCDE", "abcdef"))

at1 <- IRangesList(IRanges(1, width=0))
replaceAt(x, at1, value="++")  # prepend ++ to each sequence

at2 <- unname(split(IRanges(end=nchar(x), width=0), seq_along(x)))
replaceAt(x, at2, value="++")  # append ++ to each sequence

at3 <- IRangesList(IRanges(4, width=0))
replaceAt(x, at3, value="-+-")  # insert -+- in each sequence

at4 <- IRangesList(IRanges(c(6, 11, 13, 10, 2, 5), width=c(0, 2, 0, 0, 0, 0)),
                   IRanges(1:6, width=0),
                   IRanges(c(2, 4, 2), width=c(0, 3, 0)))
replaceAt(x, at4, value="-")
replaceAt(x, at4, value=paste0("[", seq_along(unlist(at4)), "]"))

### 'value' is recycled across sequences:
replaceAt(x, at4, value=paste0("[", 1:5, "]"))

### Testing performance with half-million single-letter insertions at random
### locations in Celegans upstream1000:
upstream1000 <- genome$upstream1000
set.seed(33)
at5 <- unname(split(IRanges(sample(1001L, 500000L, replace=TRUE), width=0),
                    sample(length(upstream1000), 500000L, replace=TRUE)))
### Takes < 1s on my machine.
system.time(upstream1000b <- replaceAt(upstream1000, at5, value="-"))


### -------------------------------------------------------------------------
### Testing insertAt()

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### On a DNAString object

x <- DNAString("ACCAGGATTYNAAG")
insertAt(x, 1L, value="++")  # prepend ++
insertAt(x, nchar(x)+1L, value="++")  # append ++
insertAt(x, 4L, value="-+-")  # insert -+-

at <- c(6L, 10L, 2L, 5L)
insertAt(x, at)  # no-op
insertAt(x, at, value="-+-")
insertAt(x, at, value=c("-", "--", "---", "----"))
value <- c("-", "--", "---", "----", "+", "++", "+++", "++++")
insertAt(x, c(at, at), value)

at <- c(5, 1, 6, 5)  # 2 insertions at position 5 
insertAt(x, at, value=c("+", "++", "-", "--"))

### Testing performance with half-million single-letter insertions at random
### locations in Celegans chrI:
library(BSgenome.Celegans.UCSC.ce2)
genome <- BSgenome.Celegans.UCSC.ce2
chrI <- genome$chrI
set.seed(33)
at <- sample(nchar(chrI)+1L, 500000L, replace=TRUE)
### Takes < 1s on my machine (64-bit Ubuntu). The naive method would take
### about 1h30 or more:
system.time(chrIb <- insertAt(chrI, at, value="-"))
matchPattern("---", chrIb)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### On a DNAStringSet object

x <- DNAStringSet(c("ACCAGGATTYNAAG", "CCAT", "ACCGAT"))
insertAt(x, 1L, value="++")  # prepend ++ to each sequence
insertAt(x, nchar(x)+1L, value="++")  # append ++ to each sequence
insertAt(x, 4L, value="-+-")  # insert -+- in each sequence

at <- list(c(6L, 10L, 2L, 5L) , 1:5 , 2L)
insertAt(x, at)  # no-op
insertAt(x, at, value="-")
### 'value' is recycled across sequences:
insertAt(x, at, value=c("-", "--", "---"))

### Testing performance with half-million single-letter insertions at random
### locations in Celegans upstream1000:
upstream1000 <- genome$upstream1000
set.seed(33)
at <- unname(split(sample(1001L, 500000L, replace=TRUE),
                   sample(length(upstream1000), 500000L, replace=TRUE)))
### Takes < 1s on my machine.
system.time(upstream1000b <- insertAt(upstream1000, at, value="-"))

### -------------------------------------------------------------------------
}  #                 <<<--- end testing block --->>>
