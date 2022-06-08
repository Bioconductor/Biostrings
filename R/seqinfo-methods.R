### =========================================================================
### seqinfo() methods
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqinfo() getter and setter for DNAStringSet objects
###

.get_DNAStringSet_seqinfo <- function(x)
{
    ans_seqnames <- names(x)
    if (is.null(ans_seqnames))
        ans_seqnames <- as.character(seq_along(x))

    ans_seqlengths <- width(x)

    ans_is_circular <- metadata(x)$is_circular
    if (is.null(ans_is_circular))
        ans_is_circular <- NA

    ans_genome <- metadata(x)$genome
    if (is.null(ans_genome))
        ans_genome <- NA

    Seqinfo(ans_seqnames, ans_seqlengths, ans_is_circular, ans_genome)
}

setMethod("seqinfo", "DNAStringSet", .get_DNAStringSet_seqinfo)

### We implement a restricted seqinfo() setter for DNAStringSet object 'x'
### that supports altering **only** the seqlevels and/or circularity flags
### and/or genome of 'seqinfo(x)'.
### It does NOT allow subsetting 'seqinfo(x)' (by dropping/reordering some
### of its seqlevels), or altering its seqlengths!
### In other words, except for their seqnames() or isCircular() or genome(),
### Seqinfo objects 'new_seqinfo' and 'old_seqinfo' must be identical. This
### is all we need to make the seqlevelsStyle() setter work on a DNAStringSet
### object.
.check_new2old_and_new_seqinfo <-
    function(new2old, new_seqinfo, old_seqinfo, context="")
{
    if (length(new_seqinfo) != length(old_seqinfo))
        stop(wmsg("the supplied 'seqinfo' must have the same ",
                  "length as the current 'seqinfo'", context))
    if (!(is.null(new2old) || identical(new2old, seq_along(new_seqinfo))))
        stop(wmsg("'new2old' can only be set to NULL or ",
                  "'seq_along(seqinfo(x))'", context))
    seqnames(old_seqinfo) <- seqnames(new_seqinfo)
    isCircular(old_seqinfo) <- isCircular(new_seqinfo)
    genome(old_seqinfo) <- genome(new_seqinfo)
    if (!identical(new_seqinfo, old_seqinfo))
        stop(wmsg("the seqlengths() of the supplied 'seqinfo' must be ",
                  "identical to those of the current 'seqinfo'", context))
}

.set_DNAStringSet_seqinfo <-
    function(x, new2old=NULL,
                pruning.mode=c("error", "coarse", "fine", "tidy"),
                value)
{
    if (!is(value, "Seqinfo"))
        stop("the supplied 'seqinfo' must be a Seqinfo object")

    context <- paste0(" when replacing the 'seqinfo' of a ",
                      classNameForDisplay(x), " object")
    pruning.mode <- match.arg(pruning.mode)
    if (pruning.mode != "error")
        stop(wmsg("'pruning.mode' is not supported", context))
    .check_new2old_and_new_seqinfo(new2old, value, seqinfo(x), context)

    new_names <- seqnames(value)
    if (identical(new_names, as.character(seq_along(x))))
        new_names <- NULL
    names(x) <- new_names

    new_is_circular <- isCircular(value)
    if (all(is.na(new_is_circular)))
        new_is_circular <- NULL
    metadata(x)$is_circular <- new_is_circular

    new_genome <- genome(value)
    if (all(is.na(new_genome)))
        new_genome <- NULL
    metadata(x)$genome <- new_genome

    x
}

setReplaceMethod("seqinfo", "DNAStringSet", .set_DNAStringSet_seqinfo)

