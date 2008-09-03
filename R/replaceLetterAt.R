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
        if (!is.numeric(at))
            stop("'at' must be a vector of integers")
        if (!is.integer(at))
            at <- as.integer(at)
        if (!is.character(letter))
            stop("'letter' must be a character vector")
        lkup <- getXStringSubtypeConversionLookup("BString", class(x))
        if (!isSingleString(if.not.extending))
            stop("'if.not.extending' must be a single string")
        if.not.extending <- match.arg(if.not.extending, c("replace", "skip", "merge", "error"))
        if (!isTRUEorFALSE(verbose))
            stop("'verbose' must be 'TRUE' or 'FALSE'")
        .Call("XString_replace_letter_at",
              x, at, letter, lkup, if.not.extending, verbose,
              PACKAGE="Biostrings")
    }
)

replaceLetterAtLoc <- function(x, loc, letter, if.not.extending="replace", verbose=FALSE)
{
    .Deprecated("replaceLetterAt")
    replaceLetterAt(x, loc, letter, if.not.extending=if.not.extending, verbose=verbose)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The ".inplaceReplaceLetterAt" function.
###
### The user should NEVER use this function!
### This function is used by the BSgenome package for injecting SNPs into the
### sequences of a BSgenome object at sequence-load time.
###

.inplaceReplaceLetterAt <- function(x, at, letter)
{
    lkup <- getXStringSubtypeConversionLookup("BString", class(x))
    .Call("XString_inplace_replace_letter_at",
          x, at, letter, lkup,
          PACKAGE="Biostrings")
}

