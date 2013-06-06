### =========================================================================
### extractAt() & replaceAt()
### -------------------------------------------------------------------------


### Extracts multiple subsequences from XString object 'x', or from the
### sequences of XStringSet object 'x', at the ranges of positions specified
### in 'at'.
setGeneric("extractAt", signature="x",
    function(x, at) standardGeneric("extractAt")
)

### Performs multiple subsequence replacements (a.k.a. substitutions) in
### XString object 'x', or in the sequences of XStringSet object 'x', at the
### ranges of positions specified in 'at'.
setGeneric("replaceAt", signature="x",
    function(x, at, value="") standardGeneric("replaceAt")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions for normalizing the 'at' argument.
###

### The integers in 'at' are interpreted as the start positions of zero-width
### ranges.
.make_Ranges_from_numeric <- function(at) IRanges(at, width=0L)

.make_RangesList_from_IntegerList <- function(at)
{
    relist(.make_Ranges_from_numeric(unlist(at, use.names=FALSE)), at)
}

.make_Ranges_from_at <- function(at)
{
    if (is.numeric(at))
        at <- .make_Ranges_from_numeric(at)
    if (!is(at, "Ranges"))
        stop("'at' must be a Ranges object (or a numeric vector containing ",
             "the start\n  positions of zero-width ranges)")
    at
}

.make_RangesList_from_at <- function(at)
{
    if (is.numeric(at))
        at <- .make_Ranges_from_numeric(at)
    if (is(at, "Ranges"))
        at <- IRangesList(at)
    if (is.list(at))
        at <- IntegerList(at)
    if (is(at, "IntegerList"))
        at <- .make_RangesList_from_IntegerList(at)
    if (!is(at, "RangesList"))
        stop("'at' must be a RangesList object (or an IntegerList object ",
             "or a list of\n  numeric vectors, containing the start ",
             "positions of zero-width ranges).\n  ",
             "Also it can be a Ranges object (or a numeric vector containing ",
             "the start\n  positions of zero-width ranges) and in that case ",
             "is interpreted as a\n  RangesList object of length 1.")
    at
}

### Checks that the ranges in 'at' are within the bounds specified by
### 'min_start' and 'max_end'.
.check_at_bounds <- function(at, min_start, max_end)
{
    at_len <- length(at)
    if (at_len == 0L)
        return()
    at_min_start <- min(start(at))
    at_max_end <- max(end(at))
    if (at_min_start < min_start || at_max_end > max_end)
        stop("some ranges in 'at' are off-limits with respect to sequence 'x'")
}

.normarg_at1 <- function(at, x)
{
    at <- .make_Ranges_from_at(at)
    .check_at_bounds(at, 1L, length(x))
    at
}

### Returns a RangesList object of the same length as 'x'.
.normarg_at2 <- function(at, x)
{
    at <- .make_RangesList_from_at(at)
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
            stop("'at' is a zero-length object but 'x' is not")
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions for normalizing the 'value' argument.
###
### 2 helper functions to bring 'value' in the same "shape" as 'at'.
###   - For .normarg_value1(): 'at' must be a Ranges object (not checked).
###     The returned 'value' is an XStringSet object obtained by coercing the
###     input 'value' to XStringSet (if necessary) and recycling it to the
###     length of 'at' (if necessary).
###   - For .normarg_value2(): 'at' must be a RangesList object (not checked).
###     The input 'value' must be an XStringSetList object, a CharacterList
###     object, or a list of character vectors.
###     The returned 'value' is an XStringSetList object obtained by coercing
###     the input 'value' to XStringSetList (if necessary) and recycling it
###     "vertically" and "horizontally" to bring it in the same "shape" as
###     'at' (i.e. same length and same elementLengths).

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
    value_len <- length(value)
    at_len <- length(at)
    if (value_len > at_len)
        stop("'value' cannot be longer than the number of replacements")
    if (value_len < at_len) {
        if (value_len == 0L)
            stop("'length(value)' is zero but the number ",
                 "of replacements is not")
        if (at_len %% value_len != 0L)
            warning("the number of replacements is not ",
                    "a multiple of 'length(value)'")
        value <- rep(value, length.out=at_len)
    }
    value
}

### 'at' is assumed to be normalized so it has the length of 'x'.
.normarg_value2 <- function(value, at, x_seqtype)
{
    value <- .make_XStringSetList_from_value(value, x_seqtype)
    ## Vertical recycling.
    value_len <- length(value)
    at_len <- length(at)  # same as length(x)
    if (value_len > at_len)
        stop("'value' cannot be longer than 'x'")
    if (value_len < at_len) {
        if (value_len == 0L)
            stop("'value' is a zero-length object but 'x' is not")
        if (at_len %% value_len != 0L)
            warning("'length(x)' is not a multiple of 'length(value)'")
        value <- rep(value, length.out=at_len)
    }
    if (at_len != 0L) {
        ## Horizontal recycling.
        value_eltlens <- elementLengths(value)
        at_eltlens <- elementLengths(at)
        is_longer <- value_eltlens > at_eltlens
        if (any(is_longer))
            stop("some list elements in 'value' (after recycling of ",
                 "'value' to the length of\n  'x') are longer than ",
                 "their corresponding list element in 'at'")
        is_shorter <- value_eltlens < at_eltlens
        if (any(is_shorter)) {
            is_empty <- value_eltlens == 0L
            if (any(is_shorter & is_empty))
                stop("some list elements in 'value' (after recycling of ",
                     "'value' to the length of\n  'x') are of length 0, ",
                     "but their corresponding list element in 'at' is not")
            is_not_singleton <- value_eltlens != 1L
            if (any(is_shorter & is_not_singleton))
                stop("some list elements in 'value' (after \"vertical\" ",
                     "recycling of 'value' i.e.\n  recycling to the length ",
                     "of 'x') are shorter than their corresponding list\n  ",
                     "element in 'at', but have a length != 1. ",
                     "Only list elements of length 1\n  can be ",
                     "\"horizontally\" recycled at the moment.")
            idx <- which(is_shorter)
            value[idx] <- relist(rep.int(unlist(value[idx], use.names=FALSE),
                                         at_eltlens[idx]),
                                 at[idx])
        }
    }
    value
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other helper functions.
###

.unlist_and_shift_at <- function(at, x)
{
    unlisted_at <- unlist(at, use.names=FALSE)
    x_len <- length(x)
    x_width <- width(x)
    at_eltlens <- elementLengths(at)
    offsets <- cumsum(c(0L, x_width[-x_len]))
    offsets <- rep.int(offsets, at_eltlens)
    shift(unlisted_at, shift=offsets)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extractAt()
###

setMethod("extractAt", "XString",
    function(x, at)
    {
        at <- .normarg_at1(at, x)
        as(Views(x, at), "XStringSet")
    }
)

setMethod("extractAt", "XStringSet",
    function(x, at)
    {
        at <- .normarg_at2(at, x)
        unlisted_x <- unlist(x, use.names=FALSE)
        unlisted_at <- .unlist_and_shift_at(at, x)
        unlisted_ans <- extractAt(unlisted_x, unlisted_at)
        names(at) <- names(x)
        relist(unlisted_ans, at)
    }
)


if (FALSE) {  #     <<<--- begin testing extractAt() --->>>

### From an XString object
### ======================

x <- BString("abcdefghijklm")
at1 <- IRanges(5:1, width=3)
extractAt(x, at1)
names(at1) <- LETTERS[22:26]
extractAt(x, at1)

at2 <- IRanges(c(1, 5, 12), c(3, 4, 12), names=c("X", "Y", "Z"))
extractAt(x, at2)
extractAt(x, rev(at2))

### From an XStringSet object
### =========================

### Only compare the classes, the shapes (i.e. lengths + elementLengths +
### names), the inner names, and the sequence contents. Doesn't look at
### the metadata or the metadata columns (outer or inner).
identical_XStringSetList <- function(target, current)
{
    ok1 <- identical(class(target), class(current))
    ok2 <- identical(elementLengths(target), elementLengths(current))
    unlisted_target <- unlist(target, use.names=FALSE)
    unlisted_current <- unlist(current, use.names=FALSE)
    ok3 <- identical(names(unlisted_target), names(unlisted_current))
    ok4 <- all(unlisted_target == unlisted_current)
    ok1 && ok2 && ok3 && ok4
}

x <- BStringSet(c(seq1="AAAA", seq2="abcdefghijk", seq3="XYZ"))
at <- IRangesList(IRanges(4, 4),
                  IRanges(c(2, 6:3), 5),
                  IRanges(c(2, 1), c(3, 1)))
### Set the inner names.
unlisted_at <- unlist(at)
names(unlisted_at) <- paste0("rg", sprintf("%02d", seq_along(unlisted_at)))
at <- relist(unlisted_at, at)

res1a <- extractAt(x, at)
res1b <- as(mapply(extractAt, x, at), "List")
stopifnot(identical_XStringSetList(res1a, res1b))

res2a <- extractAt(x, at[3])
res2b <- as(mapply(extractAt, x, at[3]), "List")
stopifnot(identical_XStringSetList(res2a, res2b))
res2c <- extractAt(x, at[[3]])
stopifnot(identical_XStringSetList(res2a, res2c))

### Testing performance with half-million small extractions at random
### locations in Celegans upstream1000:
library(BSgenome.Celegans.UCSC.ce2)
genome <- BSgenome.Celegans.UCSC.ce2
upstream1000 <- genome$upstream1000
set.seed(33)
NE <- 500000L  # total number of extractions
sample_size <- NE * 1.1
some_ranges <- IRanges(sample(1001L, sample_size, replace=TRUE),
                       width=sample(0:75, sample_size, replace=TRUE))
some_ranges <- head(some_ranges[end(some_ranges) <= 1000L], n=NE)
split_factor <- factor(sample(length(upstream1000), NE, replace=TRUE),
                       levels=seq_along(upstream1000))
at <- unname(split(some_ranges, split_factor))
### Takes < 1s on my machine.
system.time(res3a <- extractAt(upstream1000, at))
### This is about 20x slower than the above on my machine.
system.time(res3b <- as(mapply(extractAt, upstream1000, at), "List"))
stopifnot(identical_XStringSetList(res3a, res3b))

}  #                <<<--- end testing extractAt() --->>>


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### replaceAt()
###

setMethod("replaceAt", "XString",
    function(x, at, value="")
    {
        at <- .normarg_at1(at, x)
        value <- .normarg_value1(value, at, seqtype(x))
        NR <- length(at)  # same as length(value) -- nb of replacements
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
        preserved <- IRanges(c(1L, end2 + 1L), c(start2 - 1L, length(x)))
        x_preserved <- as(Views(x, preserved), "XStringSet")
        last_value2 <- as("", class(value2))
        unlist(xscat(x_preserved, c(value2, last_value2)))
    }
)

setMethod("replaceAt", "XStringSet",
    function(x, at, value="")
    {
        at <- .normarg_at2(at, x)
        value <- .normarg_value2(value, at, seqtype(x))
        unlisted_x <- unlist(x, use.names=FALSE)
        unlisted_at <- .unlist_and_shift_at(at, x)
        unlisted_value <- unlist(value, use.names=FALSE)
        unlisted_ans <- replaceAt(unlisted_x, unlisted_at,
                                  value=unlisted_value)
        delta_width <- width(unlisted_value) - width(unlisted_at)
        ans_width <- width(x) + sum(relist(delta_width,
                                           PartitioningByEnd(at)))
        ans <- as(successiveViews(unlisted_ans, ans_width), "XStringSet")
        names(ans) <- names(x)
        ans
    }
)


if (FALSE) {  #     <<<--- begin testing replaceAt() --->>>

### On an XString object
### ====================

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
at4 <- c(5, 1, 6, 5)  # 2 insertions before position 5 
replaceAt(x, at4, value=c("+", "-", "*", "/"))

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
### Takes < 2s on my machine (64-bit Ubuntu).
system.time(current <- replaceAt(chrI, at4, Views(chrI, at4)))
stopifnot(current == chrI)

### Testing performance with half-million single-letter insertions at random
### locations in Celegans chrI:
at5 <- randomDisjointRanges(1L, nchar(chrI), 0L, 0L, 500000L)
### Takes < 2s on my machine (64-bit Ubuntu).
system.time(current <- replaceAt(chrI, at5, value="-"))
m <- matchPattern("-", current)
stopifnot(identical(sort(start(at5)), start(m) - seq_along(at5) + 1L))

system.time(current2 <- replaceAt(chrI, start(at5), value="-"))
stopifnot(identical(current, current2))

matchPattern("---", current)

### On an XStringSet object
### =======================

### Only compare the classes, the lengths, the names, and the sequence
### contents. Doesn't look at the metadata or the metadata columns.
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
### Set the inner names.
unlisted_at <- unlist(at)
names(unlisted_at) <- paste0("rg", sprintf("%02d", seq_along(unlisted_at)))
at <- relist(unlisted_at, at)

current <- replaceAt(x, at, value=extractAt(x, at))  # no-op
stopifnot(identical_XStringSet(x, current))

res1a <- replaceAt(x, at)  # deletions
res1b <- mendoapply(replaceAt, x, at)
stopifnot(identical_XStringSet(res1a, res1b))

x <- BStringSet(c("abcdefghijklmn", "ABCDE", "abcdef"))

replaceAt(x, 1, value="++")  # prepend ++ to sequences
replaceAt(x, as(nchar(x) + 2L, "List"), value="++")  # append ++ to sequences
replaceAt(x, 4, value="-+-")  # insert -+- in each sequence

at4 <- IRangesList(IRanges(c(6, 11, 13, 10, 2, 5), width=c(0, 2, 0, 0, 0, 0)),
                   IRanges(1:6, width=0),
                   IRanges(c(2, 4, 2), width=c(0, 3, 0)))
replaceAt(x, at4, value="-")
value4 <- relist(paste0("[", seq_along(unlist(at4)), "]"), at4)
replaceAt(x, at4, value=value4)
replaceAt(x, at4, value=as(c("+", "-", "*"), "List"))

### Testing performance with half-million single-letter insertions at random
### locations in Celegans upstream1000:
upstream1000 <- genome$upstream1000
set.seed(33)
split_factor <- factor(sample(length(upstream1000), 500000L, replace=TRUE),
                       levels=seq_along(upstream1000))
at5 <- unname(split(sample(1001L, 500000L, replace=TRUE),
                    split_factor))
### Takes < 2s on my machine.
system.time(res5a <- replaceAt(upstream1000, at5, value="-"))
### This is about 40x slower than the above on my machine.
system.time(res5b <- mendoapply(replaceAt,
                                upstream1000, as(at5, "List"), as("-", "List")))
stopifnot(identical_XStringSet(res5a, res5b))

}  #                <<<--- end testing replaceAt() --->>>

