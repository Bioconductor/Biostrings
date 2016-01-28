### =========================================================================
### extractAt() & replaceAt()
### -------------------------------------------------------------------------


### Extracts multiple subsequences from XString object 'x', or from the
### sequences of XStringSet object 'x', at the ranges of positions specified
### thru 'at'.
setGeneric("extractAt", signature="x",
    function(x, at) standardGeneric("extractAt")
)

### Performs multiple subsequence replacements (a.k.a. substitutions) in
### XString object 'x', or in the sequences of XStringSet object 'x', at the
### ranges of positions specified thru 'at'.
setGeneric("replaceAt", signature="x",
    function(x, at, value="") standardGeneric("replaceAt")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions for checking that the ranges in a Ranges or RangesList
### object are within specified limits.
###

### Checks that the ranges in 'at' are within the limits specified by single
### integer values 'min_start' and 'max_end'.
.is_within_limits1 <- function(at, min_start, max_end)
{
    stopifnot(is(at, "Ranges"))
    stopifnot(isSingleInteger(min_start))
    stopifnot(isSingleInteger(max_end))

    at_len <- length(at)
    if (at_len == 0L)
        return(TRUE)
    at_min_start <- min(start(at))
    at_max_end <- max(end(at))
    at_min_start >= min_start && at_max_end <= max_end
}

### For all valid 'i', checks that the ranges in 'at[[i]]' are within the
### limits specified by 'limits[i]'.
.is_within_limits2 <- function(at, limits)
{
    stopifnot(is(at, "RangesList"))
    stopifnot(is(limits, "Ranges"))
    stopifnot(length(at) == length(limits))

    unlisted_at <- unlist(at, use.names=FALSE)
    tmp <- rep.int(limits, elementNROWS(at))
    min_starts <- start(tmp)
    max_ends <- end(tmp)
    all(start(unlisted_at) >= min_starts) && all(end(unlisted_at) <= max_ends)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions for "vertical" and "horizontal" recycling.
###
### TODO: This stuff is very generic. Move it to IRanges so it can be shared
### and re-used across the IRanges/GenomicRanges/XVector/Biostrings
### infrastructure.
###

.wrap_msg <- function(...)
    paste0(strwrap(paste0(...)), collapse="\n  ")

### Vertical recycling (of any vector-like object).
.V_recycle <- function(x, skeleton, x_what, skeleton_what)
{
    x_len <- length(x)
    skeleton_len <- length(skeleton)
    if (x_len == skeleton_len)
        return(x)
    if (x_len > skeleton_len && x_len != 1L)
        stop(.wrap_msg(
            "'", x_what, "' cannot be longer than ", skeleton_what
        ))
    if (x_len == 0L)
        stop(.wrap_msg(
            "'", x_what, "' is a zero-length object but ", skeleton_what,
            " is not zero"
        ))
    if (skeleton_len %% x_len != 0L)
        warning(.wrap_msg(
            skeleton_what, " is not a multiple of 'length(", x_what, ")'"
        ))
    rep(x, length.out=skeleton_len)
}

### Horizontal recycling (of a list-like object only).
.H_recycle <- function(x, skeleton, x_what, skeleton_what, more_blahblah=NA)
{
    ## TODO: Remove this when utils::relist() is fixed.
    ## See https://stat.ethz.ch/pipermail/r-devel/2013-June/066780.html for
    ## the original bug report.
    if (is.list(skeleton))
        stop(.wrap_msg(
            "because of a bug in utils::relist(), 'skeleton' cannot be ",
            "a list at the moment. Please use a List object instead ",
            "(e.g. by passing 'as(skeleton, \"List\")' instead of 'skeleton')."
        ))

    stopifnot(is.list(x) || is(x, "List"))
    stopifnot(is.list(skeleton) || is(skeleton, "List"))
    x_len <- length(x)
    skeleton_len <- length(skeleton)
    stopifnot(x_len == skeleton_len)

    x_what2 <- paste0("some list elements in '", x_what, "'")
    if (!is.na(more_blahblah))
        x_what2 <- paste0(x_what2, " (", more_blahblah, ")")

    x_eltNROWS <- unname(elementNROWS(x))
    skeleton_eltNROWS <- unname(elementNROWS(skeleton))
    idx <- which(x_eltNROWS != skeleton_eltNROWS)
    if (length(idx) == 0L)
        return(x)

    longer_idx <- which(x_eltNROWS > skeleton_eltNROWS)
    shorter_idx <- which(x_eltNROWS < skeleton_eltNROWS)
    if (length(longer_idx) == 0L && length(shorter_idx) == 0L)
        return(x)
    if (length(longer_idx) != 0L) {
        if (max(x_eltNROWS[longer_idx]) >= 2L)
            stop(.wrap_msg(
                x_what2, " are longer than their corresponding ",
                "list element in '", skeleton_what, "'"
            ))
    }
    if (length(shorter_idx) != 0L) {
        tmp <- x_eltNROWS[shorter_idx]
        if (min(tmp) == 0L)
            stop(.wrap_msg(
                x_what2, " are of length 0, but their corresponding ",
                "list element in '", skeleton_what, "' is not"
            ))
        if (max(tmp) >= 2L)
            stop(.wrap_msg(
                x_what2, " are shorter than their corresponding ",
                "list element in '", skeleton_what, "', but have ",
                "a length >= 2. \"Horizontal\" recycling only supports ",
                "list elements of length 1 at the moment."
            ))
    }

    ## From here 'x[idx]' is guaranteed to contain list elements of length 1.

    ## We use an "unlist => stretch => relist" algo to perform the horizontal
    ## recycling. Because of this, the returned value is not necessary of the
    ## same class as 'x' (e.g. can be an IntegerList if 'x' is an ordinary
    ## list of integers and 'skeleton' a List object).
    unlisted_x <- unlist(x, use.names=FALSE)
    times <- rep.int(1L, length(unlisted_x))
    idx2 <- cumsum(x_eltNROWS)[idx]
    times[idx2] <- skeleton_eltNROWS[idx]
    unlisted_ans <- rep.int(unlisted_x, times)
    ans <- relist(unlisted_ans, skeleton)
    names(ans) <- names(x)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions for normalizing the 'at' argument.
###

### The integers in 'at' are interpreted as the start positions of zero-width
### ranges.
.make_IRanges_from_numeric <- function(at) IRanges(at, width=0L)

.make_CompressedIRangesList_from_IntegerList <- function(at)
{
    relist(.make_IRanges_from_numeric(unlist(at, use.names=FALSE)), at)
}

.make_IRanges_from_at <- function(at)
{
    if (is.numeric(at))
        at <- .make_IRanges_from_numeric(at)
    if (!is(at, "Ranges"))
        stop(.wrap_msg(
            "'at' must be a Ranges object (or a numeric vector containing ",
            "the start positions of zero-width ranges)"
        ))
    as(at, "IRanges", strict=FALSE)
}

.make_CompressedIRangesList_from_at <- function(at)
{
    if (is.numeric(at))
        at <- .make_IRanges_from_numeric(at)
    if (is(at, "Ranges"))
        at <- IRangesList(at)
    if (is.list(at))
        at <- IntegerList(at)
    if (is(at, "IntegerList"))
        at <- .make_CompressedIRangesList_from_IntegerList(at)
    if (!is(at, "RangesList"))
        stop(.wrap_msg(
            "'at' must be a RangesList object (or an IntegerList object ",
            "or a list of numeric vectors, containing the start positions ",
            "of zero-width ranges). ",
            "Also it can be a Ranges object (or a numeric vector containing ",
            "the start positions of zero-width ranges) and in that case ",
            "is interpreted as a RangesList object of length 1."
        ))
    as(at, "CompressedIRangesList", strict=FALSE)
}

.normarg_at1 <- function(at, x)
{
    at <- .make_IRanges_from_at(at)
    if (!.is_within_limits1(at, 1L, length(x)))
        stop("some ranges in 'at' are off-limits with respect to sequence 'x'")
    at
}

### Returns a RangesList object of the same length as 'x'.
.normarg_at2 <- function(at, x)
{
    at <- .make_CompressedIRangesList_from_at(at)
    if (!is.null(names(at))) {
        names(at) <- NULL
        warning("'at' names were ignored")
    }
    at <- .V_recycle(at, x, "at", "'length(x)'")
    if (!.is_within_limits2(at, IRanges(1L, width(x))))
        stop(.wrap_msg(
            "some ranges in 'at' are off-limits with respect to ",
            "their corresponding sequence in 'x'"
        ))
    at
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions for normalizing the 'value' argument.
###

.make_XStringSet_from_value <- function(value, x_seqtype)
{
    if (!is(value, "XStringSet")) {
        value_class <- paste0(x_seqtype, "StringSet")
        value <- try(as(value, value_class), silent=TRUE)
        if (is(value, "try-error"))
            stop("failed to coerce 'value' to a ", value_class, " object")
    } else if (seqtype(value) != x_seqtype) {
        seqtype(value) <- x_seqtype
    }
    value
}

.make_XStringSetList_from_value <- function(value, x_seqtype)
{
    if (is.character(value)) {
        value_class <- paste0(x_seqtype, "StringSet")
        value <- as(value, value_class)
    }
    if (is(value, "XStringSet"))
        value <- relist(value, list(seq_along(value)))
    if (is.list(value))
        value <- CharacterList(value)
    if (is(value, "CharacterList")) {
        unlisted_value <- unlist(value, use.names=FALSE)
        unlisted_value <- .make_XStringSet_from_value(unlisted_value, x_seqtype)
        value <- relist(unlisted_value, value)
    }
    if (!is(value, "XStringSetList"))
        stop("invalid type of 'value'")
    if (seqtype(value) != x_seqtype)
        seqtype(value) <- x_seqtype
    value
}

.normarg_value1 <- function(value, at, x_seqtype)
{
    value <- .make_XStringSet_from_value(value, x_seqtype)
    .V_recycle(value, at, "value", "the number of replacements")
}

### 'at' is assumed to be normalized so it has the length of 'x'.
.normarg_value2 <- function(value, at, x_seqtype)
{
    value <- .make_XStringSetList_from_value(value, x_seqtype)
    value <- .V_recycle(value, at, "value", "'length(x)'")
    .H_recycle(value, at, "value", "at",
        "after recycling of 'at' and 'value' to the length of 'x'")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extractAt()
###

setMethod("extractAt", "XString",
    function(x, at)
    {
        at <- .make_IRanges_from_at(at)
        ## extractList() will check that all the ranges in 'at' are within
        ## the limits of sequence 'x'.
        extractList(x, at)
    }
)

setMethod("extractAt", "XStringSet",
    function(x, at)
    {
        at <- .normarg_at2(at, x)
        at_eltNROWS <- elementNROWS(at)
        x2 <- rep.int(unname(x), at_eltNROWS)
        unlisted_at <- unlist(at, use.names=FALSE)
        unlisted_ans <- subseq(x2, start=start(unlisted_at),
                                   width=width(unlisted_at))
        names(unlisted_ans) <- names(unlisted_at)
        ans <- relist(unlisted_ans, at)
        names(ans) <- names(x)
        ans
    }
)


if (FALSE) {  #     <<<--- begin testing extractAt() --->>>

library(Biostrings)

### From an XStringSet object
### =========================

### Only compare the classes, the shapes (i.e. lengths + elementNROWS +
### names), the inner names, and the sequence contents. Doesn't look at
### the metadata or the metadata columns (outer or inner).
identical_XStringSetList <- function(target, current)
{
    ok1 <- identical(class(target), class(current))
    ok2 <- identical(elementNROWS(target), elementNROWS(current))
    unlisted_target <- unlist(target, use.names=FALSE)
    unlisted_current <- unlist(current, use.names=FALSE)
    ok3 <- identical(names(unlisted_target), names(unlisted_current))
    ok4 <- all(unlisted_target == unlisted_current)
    ok1 && ok2 && ok3 && ok4
}

x <- BStringSet(c(seq1="ABCD", seq2="abcdefghijk", seq3="XYZ"))
at <- IRangesList(IRanges(c(2, 1), c(3, 0)),
                  IRanges(c(7, 2, 12, 7), c(6, 5, 11, 8)),
                  IRanges(2, 2))
### Set inner names on 'at'.
unlisted_at <- unlist(at)
names(unlisted_at) <- paste0("rg", sprintf("%02d", seq_along(unlisted_at)))
at <- relist(unlisted_at, at)

res1a <- extractAt(x, at)
res1b <- BStringSetList(mapply(extractAt, x, at))
stopifnot(identical_XStringSetList(res1a, res1b))

res2a <- extractAt(x, at[3])
res2b <- BStringSetList(mapply(extractAt, x, at[3]))
stopifnot(identical_XStringSetList(res2a, res2b))
res2c <- extractAt(x, at[[3]])
stopifnot(identical_XStringSetList(res2a, res2c))

### Testing performance with 2 millions small extractions at random
### locations in Fly upstream sequences:
dm3_upstream_filepath <- system.file("extdata",
                                     "dm3_upstream2000.fa.gz",
                                     package="Biostrings")
dm3_upstream <- readDNAStringSet(dm3_upstream_filepath)
dm3_upstream <- dm3_upstream[width(dm3_upstream) == 2000L]

set.seed(33)
NE <- 2000000L  # total number of extractions
sample_size <- NE * 1.1
some_ranges <- IRanges(sample(2001L, sample_size, replace=TRUE),
                       width=sample(0:75, sample_size, replace=TRUE))
some_ranges <- head(some_ranges[end(some_ranges) <= 2000L], n=NE)
split_factor <- factor(sample(length(dm3_upstream), NE, replace=TRUE),
                       levels=seq_along(dm3_upstream))
at <- unname(split(some_ranges, split_factor))

### Takes < 1s on my machine.
system.time(res3a <- extractAt(dm3_upstream, at)) 

### Equivalent but about 250x slower than the above on my machine.
system.time(res3b <- DNAStringSetList(mapply(extractAt, dm3_upstream, at)))
stopifnot(identical_XStringSetList(res3a, res3b))

}  #                <<<--- end testing extractAt() --->>>


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### replaceAt()
###

setMethod("replaceAt", "XString",
    function(x, at, value="")
    {
        if (length(at) == 0L && length(value) == 0L)
            return(x)
        at <- .normarg_at1(at, x)
        value <- .normarg_value1(value, at, seqtype(x))
        NR <- length(at)  # same as length(value) -- nb of replacements
        if (NR == 0L)
            return(x)
        .Call2("XString_replaceAt", x, at, value,
               PACKAGE="Biostrings")
    }
)

setMethod("replaceAt", "XStringSet",
    function(x, at, value="")
    {
        if (length(at) == 0L && length(value) == 0L)
            return(x)
        at <- .normarg_at2(at, x)
        value <- .normarg_value2(value, at, seqtype(x))
        ans <- .Call2("XStringSet_replaceAt", x, at, value,
                      PACKAGE="Biostrings")
        names(ans) <- names(x)
        mcols(ans) <- mcols(x)
        ans
    }
)


if (FALSE) {  #     <<<--- begin testing replaceAt() --->>>

library(Biostrings)

### On an XString object
### ====================

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

system.time(current2 <- replaceAt(chrI, start(at5), value="-"))
stopifnot(identical(current, current2))

matchPattern("---", current)

### On an XStringSet object
### =======================

### Only compare the classes, lengths, names, and sequence contents.
### Doesn't look at the metadata or the metadata columns.
identical_XStringSet <- function(target, current)
{
    ok1 <- identical(class(target), class(current))
    ok2 <- identical(names(target), names(current))
    ok3 <- all(target == current)
    ok1 && ok2 && ok3
}

x <- BStringSet(c(seq1="ABCD", seq2="abcdefghijk", seq3="XYZ"))
at <- IRangesList(IRanges(c(2, 1), c(3, 0)),
                  IRanges(c(7, 2, 12, 7), c(6, 5, 11, 8)),
                  IRanges(2, 2))
### Set inner names on 'at'.
unlisted_at <- unlist(at)
names(unlisted_at) <- paste0("rg", sprintf("%02d", seq_along(unlisted_at)))
at <- relist(unlisted_at, at)

current <- replaceAt(x, at, value=extractAt(x, at))  # no-op
stopifnot(identical_XStringSet(x, current))

res1a <- replaceAt(x, at)  # deletions
res1b <- mendoapply(replaceAt, x, at)
stopifnot(identical_XStringSet(res1a, res1b))

### Testing performance with half-million single-letter insertions at random
### locations in Fly upstream sequences:
set.seed(33)
split_factor <- factor(sample(length(dm3_upstream), 500000L, replace=TRUE),
                       levels=seq_along(dm3_upstream))
at5 <- unname(split(sample(2001L, 500000L, replace=TRUE),
                    split_factor))
### Takes < 1s on my machine.
system.time(res5a <- replaceAt(dm3_upstream, at5, value="-"))
### Equivalent but about 1400x slower than the above on my machine.
system.time(res5b <- mendoapply(replaceAt,
                                dm3_upstream, as(at5, "List"), as("-", "List")))
stopifnot(identical_XStringSet(res5a, res5b))

}  #                <<<--- end testing replaceAt() --->>>

