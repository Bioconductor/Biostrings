### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "replaceLetterAt" generic function and methods.
###

setGeneric("replaceLetterAt", signature="x",
    function(x, at, letter, if.not.extending="replace", verbose=FALSE)
        standardGeneric("replaceLetterAt")
)

setMethod("replaceLetterAt", "DNAString",
    function(x, at, letter, if.not.extending="replace", verbose=FALSE)
    {
        if (is(at, "Rle"))
            at <- as.vector(at)
        if (is.logical(at)) {
            if (length(at) != length(x))
                stop("when 'at' is a logical sequence, it must have the ",
                     "same length as 'x'")
            at <- which(at)
        } else {
            if (!is.numeric(at))
                stop("'at' must be a vector of integers")
            if (!is.integer(at))
                at <- as.integer(at)
        }
        if (is(letter, "DNAString"))
            letter <- as.character(letter)
        else if (!is.character(letter))
            stop("'letter' must be a DNAString object or a character vector")
        lkup <- get_seqtype_conversion_lookup("B", seqtype(x))
        if (!isSingleString(if.not.extending))
            stop("'if.not.extending' must be a single string")
        if.not.extending <- match.arg(if.not.extending, c("replace", "skip", "merge", "error"))
        if (!isTRUEorFALSE(verbose))
            stop("'verbose' must be TRUE or FALSE")
        .Call2("XString_replace_letter_at",
              x, at, letter, lkup, if.not.extending, verbose,
              PACKAGE="Biostrings")
    }
)

## Current restrictions: 'x' and 'at' must be rectangular i.e. 'x' must have
## a constant width and 'at' must be a logical matrix.
## TODO: Get rid of these restrictions.
setMethod("replaceLetterAt", "DNAStringSet",
    function(x, at, letter, if.not.extending="replace", verbose=FALSE)
    {
        if (length(x) == 0L)
            stop("'x' has no element")
        x_width <- width(x)
        if (!isConstant(x_width))
            stop("'x' must be rectangular (i.e. have a constant width)")
        if (!is.logical(at) || !is.matrix(at))
            stop("'at' must be a matrix of logicals")
        if (nrow(at) != length(x) || ncol(at) != x_width[1])
            stop("'x' and 'at' must have the same dimensions")
        if (is(letter, "DNAStringSet"))
            letter <- as.character(letter)
        else if (!is.character(letter))
            stop("'letter' must be a DNAStringSet object or a character vector")
        if (length(letter) != length(x))
            stop("'x' and 'letter' must have the same length")
        if (!all(width(letter) == rowSums(at)))
            stop("width(letter) and rowSums(at) must be the same")
        unlisted_x <- unlist(x, use.names=FALSE)
        unlisted_ans <- replaceLetterAt(unlisted_x, as.vector(t(at)), letter,
                                        if.not.extending=if.not.extending,
                                        verbose=verbose)
        relist(unlisted_ans, x)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The ".inplaceReplaceLetterAt" function.
###
### The user should NEVER use this function!
### This function is used by the BSgenome package for injecting SNPs into the
### sequences of a BSgenome object at sequence-load time.
###

.inplaceReplaceLetterAt <- function(x, at, letter)
{
    lkup <- get_seqtype_conversion_lookup("B", seqtype(x))
    .Call2("XString_inplace_replace_letter_at",
          x, at, letter, lkup,
          PACKAGE="Biostrings")
}

