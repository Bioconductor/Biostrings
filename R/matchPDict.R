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
###   length: the dictionary length L i.e. the number of patterns (e.g. L=10^6)
###
###   width: the dictionary width W i.e. the number of chars per pattern
###       (e.g. W=25)
###
###   ACtree: an external integer vector (XInteger object) for the storage
###       of the Aho-Corasick 4-ary tree built from the input dictionary.
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
        length="integer",
        width="integer",
        ACtree="XInteger",
        AC_base_codes="integer",
        dups="integer"
    )
)

debug_ULdna <- function()
{
    invisible(.Call("match_ULdna_debug", PACKAGE="Biostrings"))
}

.ULdna_PDict.postinit <- function(.Object, length, slotvals)
{
    .Object@length <- length
    .Object@width <- slotvals$width
    .Object@ACtree <- XInteger(1)
    .Object@ACtree@xp <- slotvals$ACtree_xp
    .Object@AC_base_codes <- slotvals$AC_base_codes
    .Object@dups <- slotvals$dups
    .Object
}

### 'dict' must be a string vector (aka character vector) with at least 1
### element.
.ULdna_PDict.init_with_StrVect <- function(.Object, dict)
{
    if (any(is.na(dict)))
        stop("'dict' contains NAs")
    slotvals <- .Call("ULdna_init_with_StrVect", dict, PACKAGE="Biostrings")
    .ULdna_PDict.postinit(.Object, length(dict), slotvals)
}

### 'dict' must be a list of BString objects of the same class (i.e. all
### BString instances or all DNAString instances or etc...) with at least 1
### element.
.ULdna_PDict.init_with_BStringList <- function(.Object, dict)
{
    dict0 <- lapply(dict, function(pattern) list(pattern@data@xp, pattern@offset, pattern@length))
    slotvals <- .Call("ULdna_init_with_BStringList", dict0, PACKAGE="Biostrings")
    .ULdna_PDict.postinit(.Object, length(dict), slotvals)
}

### 'dict' must be a BStringViews object.
.ULdna_PDict.init_with_BStringViews <- function(.Object, dict)
{
    if (length(dict) == 0)
        stop("'dict' has no views")
    slotvals <- .Call("ULdna_init_with_views",
                  subject(dict)@data@xp, subject(dict)@offset, subject(dict)@length,
                  start(dict), end(dict),
                  PACKAGE="Biostrings")
    .ULdna_PDict.postinit(.Object, length(dict), slotvals)
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

setMethod("length", "ULdna_PDict", function(x) x@length)

setMethod("width", "ULdna_PDict", function(x) x@width)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "PDictMatches" class.
### 
### A container for storing the matches returned by matchPDict().
###
### WARNING: This class is a WORK IN PROGRESS, it's not exported and not used
### yet! For now, matchPDict() just returns a list of integer vectors.
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
### Aho-Corasick
###
### Return a list of integer vectors.
###
### A real use-case:   
###   > library(hgu95av2probe)
###   > dict <- hgu95av2probe$sequence # the input dictionary
###   > pdict <- new("ULdna_PDict", dict)
###   > length(pdict)
###   [1] 201800
###   > width(pdict)
###   [1] 25
###   > library(BSgenome.Hsapiens.UCSC.hg18)
###   > chr1 <- Hsapiens$chr1
###   > system.time(pid2matchends <- matchPDict(pdict, chr1))
###      user  system elapsed 
###    50.663   0.000  50.763
###   > nmatches <- sapply(pid2matchends, length)
###   > table(nmatches)
###   > id0 <- which(nmatches == max(nmatches))
###   > p0 <- DNAString(dict[id0])
###   > p0
###     25-letter "DNAString" instance
###   Value: CTGTAATCCCAGCACTTTGGGAGGC
###   > subBString(chr1, pid2matchends[[id0]][1]-24, pid2matchends[[id0]][1]) == p0
###   [1] TRUE
### For a more extensive validation:
###   > pidOK <- sapply(seq_len(length(pid2matchends)),
###                     function(pid) identical(pid2matchends[[pid]],
###                                             end(matchPattern(DNAString(dict[pid]), chr1))))
###   > all(pidOK)
### but be aware that THIS WILL TAKE THE WHOLE DAY!!! (20-24 hours)
###
.match.ULdna_PDict.exact <- function(pdict, subject)
{
    ## Has pdict been generated from encoded input (DNAString views) or
    ## not (character vector or BString views)?
    is_used <- pdict@AC_base_codes != -1
    used_codes <- pdict@AC_base_codes[is_used]
    if (!all(used_codes %in% DNA_STRING_CODEC@codes)) {
        used_codes <- DNA_STRING_CODEC@enc_lkup[used_codes + 1]
        if (any(is.na(used_codes)) || any(duplicated(used_codes)))
            stop("the pattern dictionary 'pdict' is incompatible with a DNAString subject")
        pdict@AC_base_codes[is_used] <- used_codes
    }
    .Call("match_ULdna_exact",
          pdict@length, pdict@dups,
          pdict@ACtree@xp, pdict@AC_base_codes,
          subject@data@xp, subject@offset, subject@length,
          PACKAGE="Biostrings")
}

### With a big random dictionary, on george1:
###
### 1. Trying to simulate Solexa data:
###      > library(Biostrings)
###      > dict_length <- 10^6
###      > s <- DNAString(paste(sample(c("A", "C", "G", "T"), 36*dict_length, replace=TRUE), collapse=""))
###      > views_start <- (0:(dict_length-1)) * 36 + 1
###      > dict <- views(s, views_start, views_start + 35) # the input dictionary
###
### 2. Building the Aho-Corasick 4-ary tree from the input dictionary:
###      > pdict <- new("ULdna_PDict", dict)
###
### 3. Using pdict on Human chr1:
###      > library(BSgenome.Hsapiens.UCSC.hg18)
###      > chr1 <- DNAString(Hsapiens$chr1)
###      > system.time(pid2matchends <- matchPDict(pdict, chr1))
###         user  system elapsed
###      105.239   0.188 105.429
###      > nmatches <- sapply(pid2matchends, length)
###      > sum(nmatches) # most likely no match were found
###
### Results obtained with some random dictionaries on george1:
###
###     dict    dict   preprocess   pdict   searching   searching     total nb
###   length   width         time    size   chr1 time       again   of matches
###   ------   -----   ----------   -----   ---------   ---------   ----------
###       1M      36      2.5 sec    717M     106 sec     103 sec            0
###      10M      36       56 sec   6724M     351 sec     200 sec            0
###      10M      12      7.5 sec    340M     227 sec     216 sec         100M
###      30M      12       27 sec    523M     491 sec           ? 


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .matchPDict()
###

.matchPDict <- function(pdict, subject, algorithm, mismatch, fixed)
{
    if (!is(pdict, "ULdna_PDict"))
        stop("the pattern dictionary 'pdict' can only be a ULdna_PDict object for now")
    if (!is(subject, "DNAString"))
        stop("'subject' can only be a DNAString object for now")
    if (!identical(algorithm, "auto"))
        stop("'algo' can only be '\"auto\"' for now")
    if (!identical(mismatch, 0))
        stop("'mismatch' can only be set to 0 for now")
    if (!identical(fixed, TRUE))
        stop("'fixed' can only be 'TRUE' for now")
    .match.ULdna_PDict.exact(pdict, subject)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "matchPDict" generic and methods.
###

setGeneric(
    "matchPDict", signature="subject",
    function(pdict, subject, algorithm="auto", mismatch=0, fixed=TRUE)
        standardGeneric("matchPDict")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPDict", "BString",
    function(pdict, subject, algorithm, mismatch, fixed)
    {
	.matchPDict(pdict, subject, algorithm, mismatch, fixed)
    }
)

