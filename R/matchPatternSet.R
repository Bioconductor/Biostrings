### =========================================================================
### The matchPatternSet() generic & related functions
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "PatternSet" class.
###
### Very general container for the dictionary problem.
###

setClass("PatternSet", representation("VIRTUAL"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "UniLenDNAPatternSet" class.
###
### A container for storing a preprocessed uniform-length dictionary (or set)
### of DNA patterns.
###
### Slot description:
###
###   nchar: a single integer N (e.g. N=25)
###
###   length: the length L of the original dictionary (e.g. L=10^6)
###
###   AC_tree: an external integer vector (XInteger object) for the storage
###       of the Aho-Corasick 4-ary tree built from the original dictionary.
### 
###   base_codes: 4 integers, normally the DNA base internal codes (see
###       DNA_BASE_CODES), attached to the 4 internal child slots of any node
###       in the AC_tree slot (see typedef ACNode in the C code for more
###       info). 
###
###   dups: an unnamed (and eventually empty) list of integer vectors for
###       efficient (compact) representation of the duplicated words (or
###       patterns) found in the original dictionary.
###       Example:
###           list(c(1, 8, 9), c(5, 7), c(6, 13))
###       All the integer values in this list are non NAs, >= 1, <= length
###       and unique across the entire list. Each element of the list is a
###       vector of length >= 2. In addition, the values within a vector are
###       sorted (in ascending order) and the list elements are sorted by
###       ascending first value.
###

setClass("UniLenDNAPatternSet",
    contains="PatternSet",
    representation(
        nchar="integer",
        length="integer",
        AC_tree="XInteger",
        base_codes="integer",
        dups="list" 
    )
)

### Supported 'dict':
###   - character vector
###   - list of character strings or BString or DNAString or RNAString objects
###   - BStringViews object
### Typical use:
###   library(hgu95av2probe)
###   dict <- hgu95av2probe$sequence # the original dictionary
###   ps <- new("UniLenDNAPatternSet", dict)
###

.UniLenDNAPatternSet.init_with_character <- function(.Object, dict)
{
    if (any(is.na(dict)))
        stop("'dict' contains NAs")
    pattern_length <- unique(nchar(dict, type="bytes"))
    if (length(pattern_length) != 1)
        stop("all strings in 'dict' must have the same length")
    if (pattern_length == 0)
        stop("strings in 'dict' are empty")
    .Object@nchar <- pattern_length
    .Object@length <- length(dict)
    ## init <- .Call("match_AC_init_with_stringvect", PACKAGE="Biostrings")
    .Object
}

setMethod("initialize", "UniLenDNAPatternSet",
    function(.Object, dict)
    {
        if (is.character(dict)) {
            if (length(dict) == 0)
                stop("'dict' is an empty character vector")
            .Object <- .UniLenDNAPatternSet.init_with_character(.Object, dict)
	} else if (is.list(dict)) {
            if (length(dict) == 0)
                stop("'dict' is an empty list")
            pattern_class <- unique(sapply(dict, class))
            if (length(pattern_class) != 1)
                stop("all elements in 'dict' must belong to the same class")
            if (pattern_class == "character") {
                if (!all(sapply(dict, length) == 1))
                    stop("all character vectors in 'dict' must be of length 1")
                dict <- unlist(dict, recursive=FALSE, use.names=FALSE)
                .Object <- .UniLenDNAPatternSet.init_with_character(.Object, dict)
            } else if (extends(pattern_class, "BString")) {
                pattern_length <- unique(sapply(dict, nchar))
                if (length(pattern_length) != 1)
                    stop("all ", pattern_class, " objects in 'dict' must have the same length")
                .Object@nchar <- pattern_length
                .Object@length <- length(dict)
                ## init <- .Call("match_AC_init_with_BStringList", PACKAGE="Biostrings")
            } else
                stop("invalid 'dict' (type '?UniLenDNAPatternSet' for more information)")
        } else if (is(dict, "BStringViews")) {
            if (length(dict) == 0)
                stop("'dict' has no views")
            pattern_length <- width(dict)
            if (any(nchar(dict) != pattern_length))
                stop("'dict' has out of limits views")
            pattern_length <- unique(pattern_length)
            if (length(pattern_length) != 1)
                stop("all views in 'dict' must have the same width")
            .Object@nchar <- pattern_length
            .Object@length <- length(dict)
            dict <- as.list(dict)
            ## init <- .Call("match_AC_init_with_BStringList", PACKAGE="Biostrings")
        } else
            stop("invalid 'dict' (type '?UniLenDNAPatternSet' for more information)")
        .Object
    }
)

setMethod("nchar", "UniLenDNAPatternSet",
    function(x, type = "chars", allowNA = FALSE) x@nchar)

setMethod("length", "UniLenDNAPatternSet",
    function(x) x@length)

