### =========================================================================
### The maskMotif() generic & methods
### -------------------------------------------------------------------------


setGeneric("maskMotif", signature=c("x", "motif"),
    function(x, motif, min.block.width=1, ...) standardGeneric("maskMotif")
)

setMethod("maskMotif", signature(x="MaskedXString", motif="XString"),
    function(x, motif, min.block.width=1, ...)
    {
        if (!isSingleNumber(min.block.width))
            stop("'min.block.width' must be a single integer")
        if (!is.integer(min.block.width))
            min.block.width <- as.integer(min.block.width)
        nir1 <- as(matchPattern(motif, x, ...), "NormalIRanges")
        desc1 <- paste(as.character(motif), "-blocks", sep="")
        if (min.block.width > length(motif)) {
            nir1 <- nir1[width(nir1) >= min.block.width]
            desc1 <- paste(desc1, " [width>=", min.block.width, "]", sep="")
        }
        mask1 <- new2("MaskCollection", nir_list=list(nir1),
                                        width=length(x),
                                        active=TRUE,
                                        desc=desc1,
                                        check=FALSE)
        masks(x) <- append(masks(x), mask1)
        x
    }
)

setMethod("maskMotif", signature(x="MaskedXString", motif="character"),
    function(x, motif, min.block.width=1, ...)
    {
        motif <- XString(seqtype(x), motif)
        maskMotif(x, motif, min.block.width=min.block.width, ...)
    }
)

setMethod("maskMotif", signature(x="XString", motif="ANY"),
    function(x, motif, min.block.width=1, ...)
    {
        x <- as(x, paste("Masked", xsbaseclass(x), sep=""))
        maskMotif(x, motif, min.block.width=min.block.width, ...)
    }
)

### Used in Robert's book!
### mask() was introduced at the time where only XStringViews objects were
### available in Biostrings so it returns one instead of a MaskedXString
### object.
### Deprecate in Biostrings 2.9!
mask <- function(x, start=NA, end=NA, pattern)
{
    if (!is(x, "XString"))
        x <- XString(NULL, x)
    if (missing(pattern)) {
        if (isNumericOrNAs(start)) {
            if (length(start) == 1L && is.na(start))
                start <- 1L
            if (length(end) == 1L && is.na(end))
                end <- length(x)
            return(gaps(Views(x, start=start, end=end)))
        }
        if (!missing(end))
            stop("invalid 'start' argument")
        pattern <- start
    } else {
        if (!missing(start) || !missing(end))
            stop("can't give 'start' (or 'end') when 'pattern' is given")
    }
    as(maskMotif(x, pattern), "XStringViews")
}

