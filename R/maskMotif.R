### =========================================================================
### The maskMotif() generic & methods
### -------------------------------------------------------------------------


setGeneric("maskMotif", signature=c("x", "motif"),
    function(x, motif, min.block.width=1) standardGeneric("maskMotif")
)

setMethod("maskMotif", signature(x="MaskedXString", motif="XString"),
    function(x, motif, min.block.width=1)
    {
        if (!isSingleNumber(min.block.width))
            stop("'min.block.width' must be a single integer")
        if (!is.integer(min.block.width))
            min.block.width <- as.integer(min.block.width)
        nir1 <- toNormalIRanges(matchPattern(motif, unmasked(x)))
        name1 <- paste(as.character(motif), "-blocks", sep="")
        if (min.block.width > length(motif)) {
            nir1 <- nir1[width(nir1) >= min.block.width]
            name1 <- paste(name1, " [width>=", min.block.width, "]", sep="")
        }
        mask1 <- new("MaskCollection", nir_list=list(nir1),
                                       width=length(x),
                                       active=TRUE,
                                       NAMES=name1)
        masks(x) <- append(masks(x), mask1)
        x
    }
)

setMethod("maskMotif", signature(x="MaskedXString", motif="character"),
    function(x, motif, min.block.width=1)
    {
        motif <- XString(baseXStringSubtype(x), motif)
        maskMotif(x, motif, min.block.width=min.block.width)
    }
)

setMethod("maskMotif", signature(x="XString", motif="ANY"),
    function(x, motif, min.block.width=1)
    {
        x <- as(x, paste("Masked", baseXStringSubtype(x), sep=""))
        maskMotif(x, motif, min.block.width=min.block.width)
    }
)

