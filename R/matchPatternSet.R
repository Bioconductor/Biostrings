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
### The "ULdna_PatternSet" class.
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
###   AC_base_codes: 4 integers, normally the DNA base internal codes (see
###       DNA_BASE_CODES), attached to the 4 internal child slots of any node
###       in the AC_tree slot (see typedef ACNode in the C code for more
###       info). 
###
###   dups: an unnamed (and eventually empty) list of integer vectors
###       containing the indices of the duplicated words (patterns) found in
###       the original dictionary. For example:
###           list(c(6, 8), c(2, 9, 13), c(4, 12))
###       This allows efficient (compact) representation considering that the
###       number of duplicated is expected to be small.
###       All the integer values in this list are non NAs, >= 1, <= length
###       and unique across the entire list. Each element of the list is a
###       vector of length >= 2. In addition, the values within a vector are
###       sorted (in ascending order) <NO MORE TRUE> and the list elements are
###       sorted by ascending first value </NO MORE TRUE>.
###

setClass("ULdna_PatternSet",
    contains="PatternSet",
    representation(
        nchar="integer",
        length="integer",
        AC_tree="XInteger",
        AC_base_codes="integer",
        dups="list" 
    )
)

debug_ACuldna <- function()
{
    invisible(.Call("match_ACuldna_debug", PACKAGE="Biostrings"))
}

### 'dict' must be a string vector (aka character vector) with at least 1
### element.
.ULdna_PatternSet.init_with_StrVect <- function(.Object, dict)
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
    init <- .Call("ACuldna_init_with_StrVect", dict, PACKAGE="Biostrings")
    .Object@AC_tree <- XInteger(1)
    .Object@AC_tree@xp <- init$AC_tree_xp
    .Object@AC_base_codes <- init$AC_base_codes
    .Object@dups <- init$dups
    .Object
}

### 'dict' must be a list of BString objects of the same class (i.e. all
### BString instances or all DNAString instances or etc...) with at least 1
### element.
.ULdna_PatternSet.init_with_BStringList <- function(.Object, dict)
{
    pattern_length <- unique(sapply(dict, nchar))
    if (length(pattern_length) != 1)
        stop("all ", class(dict[[1]]), " objects in 'dict' must have the same length")
    .Object@nchar <- pattern_length
    .Object@length <- length(dict)
    dict0 <- lapply(dict, function(pattern) list(pattern@data@xp, pattern@offset, pattern@length))
    init <- .Call("ACuldna_init_with_BStringList", dict0, PACKAGE="Biostrings")
    .Object@AC_tree <- XInteger(1)
    .Object@AC_tree@xp <- init$AC_tree_xp
    .Object@AC_base_codes <- init$AC_base_codes
    .Object@dups <- init$dups
    .Object
}

### Supported 'dict':
###   - character vector
###   - list of character strings or BString or DNAString or RNAString objects
###   - BStringViews object
### Typical use:
###   library(hgu95av2probe)
###   dict <- hgu95av2probe$sequence # the original dictionary
###   patset <- new("ULdna_PatternSet", dict)
###
setMethod("initialize", "ULdna_PatternSet",
    function(.Object, dict)
    {
        if (is.character(dict)) {
            if (length(dict) == 0)
                stop("'dict' is an empty character vector")
            return(.ULdna_PatternSet.init_with_StrVect(.Object, dict))
	}
        if (is.list(dict)) {
            if (length(dict) == 0)
                stop("'dict' is an empty list")
            pattern_class <- unique(sapply(dict, class))
            if (length(pattern_class) != 1)
                stop("all elements in 'dict' must belong to the same class")
            if (pattern_class == "character") {
                if (!all(sapply(dict, length) == 1))
                    stop("all character vectors in 'dict' must be of length 1")
                dict0 <- unlist(dict, recursive=FALSE, use.names=FALSE)
                return(.ULdna_PatternSet.init_with_StrVect(.Object, dict0))
            }
            if (extends(pattern_class, "BString"))
                return(.ULdna_PatternSet.init_with_BStringList(.Object, dict))
        }
        if (is(dict, "BStringViews")) {
            if (length(dict) == 0)
                stop("'dict' has no views")
            pattern_length <- width(dict)
            if (any(nchar(dict) != pattern_length))
                stop("'dict' has out of limits views")
            pattern_length <- unique(pattern_length)
            if (length(pattern_length) != 1)
                stop("all views in 'dict' must have the same width")
            return(.ULdna_PatternSet.init_with_BStringList(.Object, as.list(dict)))
        }
        stop("invalid 'dict' (type '?ULdna_PatternSet' for more information)")
    }
)

setMethod("nchar", "ULdna_PatternSet",
    function(x, type = "chars", allowNA = FALSE) x@nchar)

setMethod("length", "ULdna_PatternSet",
    function(x) x@length)

