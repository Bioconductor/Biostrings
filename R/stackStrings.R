### =========================================================================
### stackStrings()
### -------------------------------------------------------------------------


### Return an XString object of length 1.
.normarg_padding.letter <- function(padding.letter, seqtype)
{
    class <- paste0(seqtype, "String")
    if (isSingleString(padding.letter) && nchar(padding.letter) == 1L)
        return(as(padding.letter, class))
    if (is(padding.letter, "XString") && length(padding.letter) == 1L)
        return(as(padding.letter, class))
    if (is(padding.letter, "XStringSet") &&
        length(padding.letter) == 1L && width(padding.letter) == 1L)
        return(as(padding.letter[[1L]], class))
    stop("'padding.letter' must be a single letter")
}

### 'filler_width' must be an integer vector and 'letter' an XString object
### of length 1.
.make_sequence_fillers_from_widths <- function(filler_width, letter)
{
    if (length(filler_width) == 0L) {
        max_width <- 0L
        at <- IRanges()
    } else {
        max_width <- max(filler_width)
        at <- IRanges(1L, filler_width)
    }
    biggest_filler <- rep.int(letter, max_width)
    extractAt(biggest_filler, at)
}

stackStrings <- function(x, from, to, padding.letter, shift=0L)
{
    if (!is(x, "XStringSet"))
        stop("'x' must be an XStringSet object")
    if (!isSingleNumber(from))
        stop("'from' must be a single integer")
    if (!is.integer(from))
        from <- as.integer(from)
    if (!isSingleNumber(to))
        stop("'to' must be a single integer")
    if (!is.integer(to))
        to <- as.integer(to)
    width0 <- to - from + 1L
    if (width0 < 0L)
        stop("'to' must be >= 'from - 1L'")
    padding.letter <- .normarg_padding.letter(padding.letter, seqtype(x))
    if (!is.numeric(shift))
        stop("'shift' must be a vector of integers")
    if (!is.integer(shift)) 
        shift <- as.integer(shift)
    ## .V_recycle() is currently defined in replaceAt.R but it needs to move
    ## to a place more appropriate for sharing (see TODO note in replaceAt.R).
    shift <- .V_recycle(shift, x, "shift", "'length(x)'")

    left_margin <- shift + 1L - from
    right_margin <- to - (shift + width(x))

    left_pad <- pmin(pmax(left_margin, 0L), width0)
    right_pad <- pmin(pmax(right_margin, 0L), width0)
    left_trim <- pmin(pmax(-left_margin, 0L), width(x))
    right_trim <- pmin(pmax(-right_margin, 0L), width(x))

    left <- .make_sequence_fillers_from_widths(left_pad, padding.letter)
    right <- .make_sequence_fillers_from_widths(right_pad, padding.letter)
    middle <- narrow(x, start=1L+left_trim, end=-(1L+right_trim),
                        use.names=FALSE)
    ans <- xscat(left, middle, right)
    names(ans) <- names(x)
    mcols(ans) <- mcols(x)
    ans
}

