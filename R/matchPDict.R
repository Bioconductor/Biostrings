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
### The "ACtree" class (NOT EXPORTED).
###
### A low-level container for storing the Aho-Corasick tree built from the
### input dictionary (note that this tree is actually an oriented graph if we
### consider the failure links). When the input is based on an alphabet of 4
### letters (DNA input) then the Aho-Corasick tree is a 4-ary tree and the
### base_codes slot is a vector of 4 integers that will be filled with the
### internal codes (byte values) of the unique letters found in the input.
### The number of integers needed to represent a node is the number of letters
### in the alphabet plus 3 (i.e. 7 for DNA input) hence this internal
### representation of the Aho-Corasick tree can only be used for input based
### on a very small alphabet.
###
### Slot description:
###
###   nodes: an external integer vector (XInteger object) for the storage
###       of the nodes.
### 
###   base_codes: an integer vector filled with the internal codes (byte
###       values) of the unique letters found in the input dictionary during
###       its preprocessing.

setClass("ACtree",
    representation(
        nodes="XInteger",
        base_codes="integer"
    )
)

.ACtree.ints_per_acnode <- function(x) (length(x@base_codes) + 3L)

setMethod("length", "ACtree", function(x) length(x@nodes) %/% .ACtree.ints_per_acnode(x))

setMethod("show", "ACtree",
    function(object)
    {
        cat(length(object), "-node ", length(object@base_codes), "-ary Aho-Corasick tree\n", sep="")
    }
)

### Implement the 0-based subsetting semantic (like in C).
### Typical use:
###   > pdict <- new("ULdna_PDict", c("agg", "aca", "gag", "caa", "agt", "aca"))
###   > pdict@actree[0:2] # look at the 3 first nodes
###   > pdict@actree[] # look at all the nodes
###   > flinks0 <- as.matrix(pdict@actree)[ , "flink"]
###   > flinks0 # no failure links is set yet
###   > mends <- matchPDict(pdict, DNAString("acaagagagt"))
###   > flinks1 <- as.matrix(pdict@actree)[ , "flink"]
###   > flinks1 # some failure links have been set
### As you can see the 'pdict' object "learns" from being used!
###   
setMethod("[", "ACtree",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        lx <- length(x)
        if (missing(i))
            i <- 0:(lx-1)
        else if (!is.integer(i))
            i <- as.integer(i)
        ints_per_acnode <- .ACtree.ints_per_acnode(x)
        ii <- rep(i * ints_per_acnode, each=ints_per_acnode) + seq_len(ints_per_acnode)
        ans <- matrix(x@nodes[ii], ncol=ints_per_acnode, byrow=TRUE)
        colnames(ans) <- c("parent", pdict@actree@base_codes,
                           "flink", "P_id")
        rownames(ans) <- i
        ans
    }
)

setMethod("as.matrix", "ACtree", function(x) x[])

setMethod("initialize", "ACtree",
    function(.Object, xp, base_codes)
    {
        .Object@nodes <- XInteger(1)
        .Object@nodes@xp <- xp
        .Object@base_codes <- base_codes
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "ULdna_PDict" class.
###
### A container for storing a preprocessed uniform-length dictionary (or set)
### of DNA patterns.
###
### Slot description:
###
###   length: the dictionary length L i.e. the number of patterns
###       (e.g. L=10^6).
###
###   width: the dictionary width W i.e. the number of chars per pattern
###       (e.g. W=25).
###
###   actree: the Aho-Corasick tree built from the input dictionary.
###
###   dups: an integer vector of length L containing the duplicate info.

setClass("ULdna_PDict",
    contains="PDict",
    representation(
        length="integer",
        width="integer",
        actree="ACtree",
        dups="integer"
    )
)

setMethod("length", "ULdna_PDict", function(x) x@length)

setMethod("width", "ULdna_PDict", function(x) x@width)

setMethod("show", "ULdna_PDict",
    function(object)
    {
        cat(length(object), "-pattern uniform-length \"PDict\" object of width ",
            width(object), "\n", sep="")
    }
)

debug_ULdna <- function()
{
    invisible(.Call("match_ULdna_debug", PACKAGE="Biostrings"))
}

### 'ppans' is the answer returned by the preprocessing functions (the
### ULdna_pp_*() functions).
.ULdna_PDict.postinit <- function(.Object, length, ppans)
{
    .Object@length <- length
    .Object@width <- ppans$width
    .Object@actree <- new("ACtree", ppans$actree_nodes_xp, ppans$actree_base_codes)
    .Object@dups <- ppans$dups
    .Object
}

### 'dict' must be a string vector (aka character vector) with at least 1
### element.
.ULdna_PDict.init_with_StrVect <- function(.Object, dict)
{
    if (any(is.na(dict)))
        stop("'dict' contains NAs")
    ppans <- .Call("ULdna_pp_StrVect", dict, PACKAGE="Biostrings")
    .ULdna_PDict.postinit(.Object, length(dict), ppans)
}

### 'dict' must be a list of BString objects of the same class (i.e. all
### BString instances or all DNAString instances or etc...) with at least 1
### element.
.ULdna_PDict.init_with_BStringList <- function(.Object, dict)
{
    dict0 <- lapply(dict, function(pattern) list(pattern@data@xp, pattern@offset, pattern@length))
    ppans <- .Call("ULdna_pp_BStringList", dict0, PACKAGE="Biostrings")
    .ULdna_PDict.postinit(.Object, length(dict), ppans)
}

### 'dict' must be a BStringViews object.
.ULdna_PDict.init_with_BStringViews <- function(.Object, dict)
{
    if (length(dict) == 0)
        stop("'dict' has no views")
    ppans <- .Call("ULdna_pp_views",
                  subject(dict)@data@xp, subject(dict)@offset, subject(dict)@length,
                  start(dict), end(dict),
                  PACKAGE="Biostrings")
    .ULdna_PDict.postinit(.Object, length(dict), ppans)
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
        on.exit(.Call("ULdna_free_actree_nodes_buf", PACKAGE="Biostrings"))
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

### Temporary trick
.ACtree.prepare_for_use_on_DNAString <- function(actree)
{
    ## Has 'actree' been generated from encoded input (DNAString views) or
    ## not (character vector or BString views)?
    is_used <- actree@base_codes != -1
    used_codes <- actree@base_codes[is_used]
    if (!all(used_codes %in% DNA_STRING_CODEC@codes)) {
        used_codes <- DNA_STRING_CODEC@enc_lkup[used_codes + 1]
        if (any(is.na(used_codes)) || any(duplicated(used_codes)))
            stop("the pattern dictionary 'pdict' is incompatible with a DNAString subject")
        actree@base_codes[is_used] <- used_codes
    }
    actree
}

.match.ULdna_PDict.exact <- function(pdict, subject)
{
    actree <- .ACtree.prepare_for_use_on_DNAString(pdict@actree)
    .Call("match_ULdna_exact",
          pdict@length, pdict@dups,
          actree@nodes@xp, actree@base_codes,
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

