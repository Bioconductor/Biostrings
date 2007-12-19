### =========================================================================
### The matchPDict() generic & related functions
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "PDict" class.
###
### Very general container for any kind of pattern dictionary.
###

setClass("PDict", representation("VIRTUAL"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "ULdna_PDict" class.
###
### A container for storing a preprocessed uniform-length dictionary (or set)
### of DNA patterns.
###
### Slot description:
###
###   width: a single integer W (e.g. W=25)
###
###   length: the length L of the original dictionary (e.g. L=10^6)
###
###   ACtree: an external integer vector (XInteger object) for the storage
###       of the Aho-Corasick 4-ary tree built from the original dictionary.
### 
###   AC_base_codes: 4 integers, normally the DNA base internal codes (see
###       DNA_BASE_CODES), attached to the 4 internal child slots of any node
###       in the ACtree slot (see typedef ACNode in the C code for more
###       info). 
###
###   dups: an integer vector of length L containing the duplicate info.
###

setClass("ULdna_PDict",
    contains="PDict",
    representation(
        width="integer",
        length="integer",
        ACtree="XInteger",
        AC_base_codes="integer",
        dups="integer"
    )
)

debug_ULdna <- function()
{
    invisible(.Call("match_ULdna_debug", PACKAGE="Biostrings"))
}

### 'dict' must be a string vector (aka character vector) with at least 1
### element.
.ULdna_PDict.init_with_StrVect <- function(.Object, dict)
{
    if (any(is.na(dict)))
        stop("'dict' contains NAs")
    pattern_length <- unique(nchar(dict, type="bytes"))
    if (length(pattern_length) != 1)
        stop("all strings in 'dict' must have the same length")
    if (pattern_length == 0)
        stop("strings in 'dict' are empty")
    .Object@width <- pattern_length
    .Object@length <- length(dict)
    init <- .Call("ULdna_init_with_StrVect", dict, PACKAGE="Biostrings")
    .Object@ACtree <- XInteger(1)
    .Object@ACtree@xp <- init$ACtree_xp
    .Object@AC_base_codes <- init$AC_base_codes
    .Object@dups <- init$dups
    .Object
}

### 'dict' must be a list of BString objects of the same class (i.e. all
### BString instances or all DNAString instances or etc...) with at least 1
### element.
.ULdna_PDict.init_with_BStringList <- function(.Object, dict)
{
    pattern_length <- unique(sapply(dict, nchar))
    if (length(pattern_length) != 1)
        stop("all ", class(dict[[1]]), " objects in 'dict' must have the same length")
    .Object@width <- pattern_length
    .Object@length <- length(dict)
    dict0 <- lapply(dict, function(pattern) list(pattern@data@xp, pattern@offset, pattern@length))
    init <- .Call("ULdna_init_with_BStringList", dict0, PACKAGE="Biostrings")
    .Object@ACtree <- XInteger(1)
    .Object@ACtree@xp <- init$ACtree_xp
    .Object@AC_base_codes <- init$AC_base_codes
    .Object@dups <- init$dups
    .Object
}

### 'dict' must be a BStringViews object.
.ULdna_PDict.init_with_BStringViews <- function(.Object, dict)
{
    if (length(dict) == 0)
        stop("'dict' has no views")
    pattern_length <- width(dict)
    if (any(nchar(dict) != pattern_length))
        stop("'dict' has out of limits views")
    pattern_length <- unique(pattern_length)
    if (length(pattern_length) != 1)
        stop("all views in 'dict' must have the same width")
    .Object@width <- pattern_length
    .Object@length <- length(dict)
    init <- .Call("ULdna_init_with_views",
                  subject(dict)@data@xp, subject(dict)@offset,
                  start(dict), end(dict),
                  PACKAGE="Biostrings")
    .Object@ACtree <- XInteger(1)
    .Object@ACtree@xp <- init$ACtree_xp
    .Object@AC_base_codes <- init$AC_base_codes
    .Object@dups <- init$dups
    .Object
}

### The input dictionary 'dict' must be:
###   - of length >= 1
###   - uniform-length (i.e. all words have the same length)
### The supported types for the input dictionary are:
###   - character vector
###   - list of character strings or BString or DNAString or RNAString objects
###   - BStringViews object
### Typical use:
###   > library(Biostrings)
###   > dict <- c("abb", "aca", "bab", "caa", "abd", "aca")
###   > pdict <- new("ULdna_PDict", dict)
###
setMethod("initialize", "ULdna_PDict",
    function(.Object, dict)
    {
        on.exit(.Call("ULdna_free_ACnodebuf", PACKAGE="Biostrings"))
        if (is.character(dict)) {
            if (length(dict) == 0)
                stop("'dict' is an empty character vector")
            return(.ULdna_PDict.init_with_StrVect(.Object, dict))
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
                return(.ULdna_PDict.init_with_StrVect(.Object, dict0))
            }
            if (extends(pattern_class, "BString"))
                return(.ULdna_PDict.init_with_BStringList(.Object, dict))
        }
        if (is(dict, "BStringViews"))
            return(.ULdna_PDict.init_with_BStringViews(.Object, dict))
        stop("invalid 'dict' (type '?ULdna_PDict' for more information)")
    }
)

setMethod("width", "ULdna_PDict", function(x) x@width)

setMethod("length", "ULdna_PDict", function(x) x@length)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "PDictMatches" class.
###
### A container for storing the matches returned by matchPDict().
###
### Slot description:
###
###   subject: the searched string.
###
###   matches: a list of integer vectors.

setClass("PDictMatches",
    representation(
        subject="BString",
        matches="list"
    )
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "matchPDict" generic and methods.
###

setGeneric(
    "matchPDict", signature="subject",
    function(pdict, subject, algorithm="auto", mismatch=0, fixed=TRUE)
        standardGeneric("matchPDict")
)

### Return a list of integer vectors.
###
### A real use-case:   
###   > library(hgu95av2probe)
###   > dict <- hgu95av2probe$sequence # the original dictionary
###   > pdict <- new("ULdna_PDict", dict)
###   > library(BSgenome.Hsapiens.UCSC.hg18)
###   > chr1 <- BString(Hsapiens$chr1) # because dict is a character vector
###   > system.time(pid2matchends <- Biostrings:::.match.ULdna_PDict.exact(pdict, chr1))
###      user  system elapsed 
###    50.663   0.000  50.763
###   > nmatches <- sapply(pid2matchends, length)
###   > table(nmatches)
###   > id0 <- which(nmatches == max(nmatches))
###   > p0 <- BString(dict[id0])
###   > p0
###     25-letter "BString" instance
###   Value: CTGTAATCCCAGCACTTTGGGAGGC
###   > subBString(chr1, pid2matchends[[id0]][1]-24, pid2matchends[[id0]][1]) == p0
###   [1] TRUE
### For a more extensive validation:
###   > pidOK <- sapply(seq_len(length(pid2matchends)), function(pid) identical(pid2matchends[[pid]], end(matchPattern(BString(dict[pid]), chr1))))
###   > all(pidOK)
### but be aware that THIS WILL TAKE THE WHOLE DAY!!! (20-24 hours)
###
### With a big random dictionary, on george1:
###   Trying to simulate Solexa data:
###   > library(Biostrings)
###   > dict_length <- 10^6
###   > s <- BString(paste(sample(c("A", "C", "G", "T"), 36*dict_length, replace=TRUE), collapse="")) # takes < 5 seconds
###   > views_start <- (0:(dict_length-1)) * 36 + 1
###   > dict <- views(s, views_start, views_start + 35)
###   Building the Aho-Corasick 4-ary tree from the input dictionary:
###   > pdict <- new("ULdna_PDict", dict) # takes < 4 seconds, size of pdict@ACtree slot is about 720M
###   Searching Human chr1:
###   > library(BSgenome.Hsapiens.UCSC.hg18)
###   > chr1 <- BString(Hsapiens$chr1) # because subject(dict) is a BString object
###   > system.time(pid2matchends <- Biostrings:::.match.ULdna_PDict.exact(pdict, chr1))
###      user  system elapsed
###   105.239   0.188 105.429
###   > nmatches <- sapply(pid2matchends, length)
###   > max(nmatches) # most likely no match were found
###
### I did the same as above but with a very big random dictionary
### (dict_length <- 10^7), on george1. It took 50 seconds to generate pdict.
### The resulting ACtree (pdict@ACtree) contained about 250M nodes and its
### size was about 6.6G. Then it took about 4 minutes to search chr1 (no match
### were found.
###
.match.ULdna_PDict.exact <- function(pdict, subject)
{
    .Call("match_ULdna_exact",
          pdict@length, pdict@dups,
          pdict@ACtree@xp, pdict@AC_base_codes,
          subject@data@xp, subject@offset, subject@length,
          PACKAGE="Biostrings")
}

