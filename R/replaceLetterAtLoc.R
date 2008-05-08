### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "replaceLetterAtLoc" generic function and methods.
###

setGeneric("replaceLetterAtLoc", signature="x",
    function(x, loc, letter, if.not.extending="replace", verbose=FALSE)
        standardGeneric("replaceLetterAtLoc")
)

setMethod("replaceLetterAtLoc", "DNAString",
    function(x, loc, letter, if.not.extending="replace", verbose=FALSE)
    {
        if (!is.numeric(loc))
            stop("'loc' must be a vector of integers")
        if (!is.integer(loc))
            loc <- as.integer(loc)
        if (!is.character(letter))
            stop("'letter' must be a character vector")
        lkup <- getXStringSubtypeConversionLookup("BString", class(x))
        if (!isSingleString(if.not.extending))
            stop("'if.not.extending' must be a single string")
        if.not.extending <- match.arg(if.not.extending, c("replace", "skip", "merge", "error"))
        if (!isTRUEorFALSE(verbose))
            stop("'verbose' must be 'TRUE' or 'FALSE'")
        .Call("XString_replace_locs_bySTRSXP",
              x, loc, letter, lkup, if.not.extending, verbose,
              PACKAGE="Biostrings")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The ".inplaceReplaceLetterAtLoc" function.
###
### The user should NEVER use this function!
### This function is used by the BSgenome package for injecting SNPs into the
### sequences of a BSgenome object at sequence-load time.
###

.inplaceReplaceLetterAtLoc <- function(x, loc, letter)
{
    lkup <- getXStringSubtypeConversionLookup("BString", class(x))
    .Call("XString_inplace_replace_locs_bySTRSXP",
          x, loc, letter, lkup,
          PACKAGE="Biostrings")
}

