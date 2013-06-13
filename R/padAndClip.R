### =========================================================================
### padAndClip() and stackStrings()
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

### 'filler_width' must be an integer vector, and 'letter' an XString object
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

padAndClip <- function(x, views, padding.letter)
{
    if (!is(x, "XStringSet"))
        stop("'x' must be an XStringSet object")
    x_width <- width(x)

    if (!is(views, "Ranges"))
        stop("'views' must be a Ranges object")
    ## .V_recycle() is currently defined in replaceAt.R but it needs to move
    ## to a place more appropriate for sharing (see TODO note in replaceAt.R).
    views_start <- .V_recycle(start(views), x, "views", "'length(x)'")
    ## We don't want to generate the same warning twice.    
    views_width <- suppressWarnings(.V_recycle(width(views), x,
                                               "views", "'length(x)'"))

    padding.letter <- .normarg_padding.letter(padding.letter, seqtype(x))

    ## Left and right margins.
    Lmargin <- 1L - views_start
    Rmargin <- views_width - x_width - Lmargin

    ## Clip. Eventhough padAndClip() conceptually pads first (with an
    ## infinite number of padding letters on the left and right) and then
    ## clips, in practice we clip first and then pad.
    Lclipping <- pmin(pmax(-Lmargin, 0L), x_width)
    Rclipping <- pmin(pmax(-Rmargin, 0L), x_width)
    clipped_x <- narrow(x, start=1L+Lclipping, end=x_width-Rclipping,
                           use.names=FALSE)

    ## Left and right padding.
    Lpad_width <- pmin(pmax(Lmargin, 0L), views_width)
    Lpad <- .make_sequence_fillers_from_widths(Lpad_width, padding.letter)
    Rpad_width <- pmin(pmax(Rmargin, 0L), views_width)
    Rpad <- .make_sequence_fillers_from_widths(Rpad_width, padding.letter)
    ans <- xscat(Lpad, clipped_x, Rpad)

    names(ans) <- names(x)
    mcols(ans) <- mcols(x)
    ans
}

### Convenience wrapper to padAndClip(). Returned object is always rectangular
### (i.e. constant-width).
stackStrings <- function(x, from, to, padding.letter, shift=0L)
{
    if (!is(x, "XStringSet"))
        stop("'x' must be an XStringSet object")

    ## Normalize 'from'.
    if (!isSingleNumber(from))
        stop("'from' must be a single integer")
    if (!is.integer(from))
        from <- as.integer(from)

    ## Normalize 'to'.
    if (!isSingleNumber(to))
        stop("'to' must be a single integer")
    if (!is.integer(to))
        to <- as.integer(to)
    if (to < from - 1L)
        stop("'to' must be >= 'from - 1L'")

    ## Normalize 'shift'.
    if (!is.numeric(shift))
        stop("'shift' must be a vector of integers")
    if (!is.integer(shift)) 
        shift <- as.integer(shift)
    ## .V_recycle() is currently defined in replaceAt.R but it needs to move
    ## to a place more appropriate for sharing (see TODO note in replaceAt.R).
    shift <- .V_recycle(shift, x, "shift", "'length(x)'")

    views <- IRanges(from - shift, to - shift)
    padAndClip(x, views, padding.letter)
}

