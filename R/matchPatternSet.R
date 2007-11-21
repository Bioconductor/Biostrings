### =========================================================================
### The matchPatternSet() generic & related functions
### -------------------------------------------------------------------------

###
setClass("PatternSet", representation("VIRTUAL"))

### Slot description
###   dups: an unnamed (and eventually empty) list of integer vectors where
###         all the integer values are non NAs, >= 1, <= npattern and unique
###         across all vectors. Each vector is of length >= 2.
###         In addition, each vector is sorted (by ascending value) and the
###         list elements are sorted by ascending first value.
###         Example:
###           list(c(1L, 8L, 9L), c(5L, 7L), c(6L, 13L))
setClass("UniLengthDNAPatternSet",
    contains="PatternSet",
    representation(
        pattern_length="integer",   # A single integer L (e.g. L=36)
        base1_code="integer",       # integer belonging to DNA_BASE_CODES
        base2_code="integer",
        base3_code="integer",
        base4_code="integer",
        npattern="integer",         # A single integer N (e.g. N=10^7)
        AC_tree="XInteger",         # Aho-Corasick 4-ary tree
        dups="list" 
    )
)

### Supported 'src':
###   - character vector
###   - list of character strings or BString or DNAString or RNAString objects
###   - BStringViews object
### Typical use:
###   library(hgu95av2probe)
###   src <- hgu95av2probe$sequence
###   dict <- new("UniLengthDNAPatternSet", src, names(Biostrings:::DNA_BASE_CODES))
###
setMethod("initialize", "UniLengthDNAPatternSet",
    function(.Object, src, base_letters)
    {
        if (is.character(src)) {
            if (length(src) == 0)
                stop("'src' is an empty character vector")
            if (any(is.na(src)))
                stop("'src' contains NAs")
            pattern_length <- unique(nchar(src, type="bytes"))
            if (length(pattern_length) != 1)
                stop("all strings in 'src' must have the same length")
            if (pattern_length == 0)
                stop("strings in 'src' are empty")
	} else if (is.list(src)) {
        } else if (is(src, "BStringViews")) {
        } else
            stop("invalid 'src' (type '?UniLengthDNAPatternSet' for more information)")
        .Object
    }
)

