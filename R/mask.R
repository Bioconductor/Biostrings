### =========================================================================
### The mask() generic and methods
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###
###
### "Normalized" BStringViews objects
### ---------------------------------
###
### Definition: BStringViews object x is "normalized" means
###  (1) no "out of limits" views,
###  (2) views are sorted from left to right (start(x) is ascending),
###  (3) views don't overlap and can't be adjacents i.e. there is at least
###      1 letter between 2 any given views.
### If length(x) >= 2, then the 3 above conditions are equivalent to:
###   1 <= start(x)[i] <= end(x)[i] <
###        start(x)[i+1] <= end(x)[i+1] <= nchar(subject(x))
### for every 1 <= i < length(x).
### If length(x) == 0, then x is normalized.
### If length(x) == 1, then x is normalized <=> view is not "out of limits".
###
### "Normalizing" a BStringViews object:
### For any given BStringViews object x, let's call S(x) the subset of
### integers defined by:
###   (union of all [start(x)[i],end(x)[i]]) inter [1,nchar(subject(x))]
### We can see that there is a unique "normalized" BStringViews object x0
### (with same subject as x) such that S(x0) = S(x).
### So we can define the normalize() function: x -> x0 = normalize(x).
### An interesting property of this function is that, for any x,
### x0 = normalize(x) is the shortest BStringViews object (with same
### subject than x) such that S(x0) = S(x).
### It can also been shown that length(x0) <= (nchar(subject(x)) + 1) / 2.
###
### Some basic operations on "normalized" BStringViews objects with same
### subject (the results of these operations are "normalized" too):
###   x2 <- !x: (unary operation) S(x2) =  [1:nchar(subject(x))] - S(x)
###   x3 <- x | y: s(x3) = (S(x) union S(y))
###   x4 <- x & y: s(x4) = (S(x) inter S(y))
### So we end up having the equivalent of the fundamental set operations
### (complementary, union, intersection).
### No need to define an equivalent for the set difference (S(x) - S(y)),
### or the set symetric difference ((S(x) -S(y)) union (S(y) - S(x)))
### since they can be achieved by using the fundamental operations.
###
### Notes:
### - !x is actually obtained with mask(x) (see below)
### - The results of the | and & operations are undefined if one of
###   the operand is not "normalized" or if the 2 operands don't have the
###   same subject! (no need to do any check of any sort, or to try to fail
###   gracefully, this will just result in slowing down the operators).
###   The result has always the same subject as the operands.
###
###
### The "mask" generic and methods
### ------------------------------
###
### For any given BStringViews object x, y <- mask(x) is the shortest (in
### term of number of views) BStringViews object such that:
###   (a) subject(y) == subject(x),
###   (b) y views cover the portions of the subject that are not covered by x
###       views,
###   (c) y views are sorted from left to right,
###   (d) y views are not "out of limits".
### Relationship with the concept of "normalized" BStringViews objects (see
### above):
###   - mask(x) is !x
###   - mask(x) is normalized.
###   - mask(mask(x)) is normalize(x).
###

setMethod("mask", "BStringViews",
    function(x, start, end, ...)
    {
        if (missing(start))
            start <- 1L
        if (missing(end))
            end <- nchar(subject(x))
        callNextMethod(x, start, end, ...)
    }
)

### We wan't to be able to do this:
###   > mask(BString("AbcbcbDE"), 2, 6)
### or this:
###   > mask(BString("AbcbcbDE"), "bcb")
setMethod("mask", "XString",
    function(x, start, end, pattern)
    {
        if (missing(pattern)) {
            if (missing(start))
                start <- NA
            if (isNumericOrNAs(start)) {
                if (missing(end))
                    end <- NA
                return(mask(views(x, start, end)))
            }
            if (!missing(end))
                stop("invalid 'start' argument")
            pattern <- start
        } else {
            if (!missing(start) || !missing(end))
                stop("can't give 'start' (or 'end') when 'pattern' is given")
        }
        mask(matchPattern(pattern, x))
    }
)

setMethod("mask", "character",
    function(x, start, end, pattern)
    {
        if (length(x) != 1 || is.na(x))
            stop("can't mask a character vector that is not a single string")
        if (missing(pattern))
            mask(BString(x), start, end)
        else
            mask(BString(x), start, end, pattern)
    }
)

