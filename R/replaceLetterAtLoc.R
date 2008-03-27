### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "replaceLetterAtLoc" generic function and methods.
###

setGeneric("replaceLetterAtLoc", signature="x",
    function(x, loc, letter, mode="replace", on.incompatible="replace", verbose=FALSE)
        standardGeneric("replaceLetterAtLoc")
)

### loc: integer vect with no NAs (locations can be repeated, in this case last replacement prevails)
### letter: character vector with no NAs
setMethod("replaceLetterAtLoc", "DNAString",
    function(x, loc, letter, mode="replace", on.incompatible="replace", verbose=FALSE)
    {
        if (!is.numeric(loc))
            stop("'loc' must be a vector of integers")
        if (!is.integer(loc))
            loc <- as.integer(loc)
        if (!is.character(letter))
            stop("'letter' must be a character vector")
        if (!isSingleString(mode))
            stop("'mode' must be a single string")
        lkup <- getXStringSubtypeConversionLookup("BString", class(x))
        mode <- match.arg(mode, c("replace", "merge"))
        if (!isSingleString(on.incompatible))
            stop("'on.incompatible' must be a single string")
        on.incompatible <- match.arg(on.incompatible, c("replace", "skip", "error"))
        if (!isTRUEorFALSE(verbose))
            stop("'verbose' must be TRUE or FALSE")
        .Call("XString_replace_locs_bySTRSXP",
              x, loc, letter, lkup, mode, on.incompatible, verbose,
              PACKAGE="Biostrings")
    }
)

