### =========================================================================
### The matchPDict() generic & related classes and functions
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "P2Views" API.
### 
### An API for manipulating the match container returned by matchPDict().
###

setClass("P2Views", representation("VIRTUAL"))

setGeneric("pids", function(x) standardGeneric("pids"))


setGeneric("p2starts", signature="x",
    function(x, all.pids=FALSE) standardGeneric("p2starts"))

setMethod("p2starts", "P2Views",
    function(x, all.pids=FALSE)
    {
        if (!is.null(pids(x)))
            return(lapply(p2ends(x, all.pids=all.pids), function(end) { end - x@width + 1L }))
        if (!missing(all.pids))
            warning("'all.pids' is ignored when \"P2Views\" object has no pattern IDs")
        lapply(p2ends(x), function(end) { end - x@width + 1L })
    }
)

setGeneric("p2ends", signature="x",
    function(x, all.pids=FALSE) standardGeneric("p2ends"))

setGeneric("p2nmatch", signature="x",
    function(x, all.pids=FALSE) standardGeneric("p2nmatch"))

setMethod("p2nmatch", "P2Views",
    function(x, all.pids=FALSE)
    {
        if (!is.null(pids(x))) {
            ends <- p2ends(x, all.pids=all.pids)
            if (length(ends) == 0)
                return(integer(0))
            return(sapply(ends, length))
        }
        if (!missing(all.pids))
            warning("'all.pids' is ignored when \"P2Views\" object has no pattern IDs")
        ends <- p2ends(x)
        if (length(ends) == 0)
            return(integer(0))
        sapply(ends, length)
    }
)

### Return a single integer or string (not NA).
setMethod("[[", "P2Views",
    function(x, i, j, ...)
    {
        # 'x' is guaranteed to be a "P2Views" object (if it's not, then
        # the method dispatch algo would not have called this method in the
        # first place), so nargs() is guaranteed to be >= 1
        if (nargs() >= 3)
            stop("too many subscripts")
        subscripts <- list(...)
        if (!missing(i))
            subscripts$i <- i
        if (!missing(j))
            subscripts$j <- j
        # At this point, 'subscripts' should be guaranteed
        # to be of length <= 1
        if (length(subscripts) == 0)
            stop("no index specified")
        key <- subscripts[[1]]
        if (!is.character(key) && !is.numeric(key))
            stop("wrong argument for subsetting an object of class \"P2Views\"")
        if (length(key) < 1)
            stop("attempt to select less than one element")
        if (length(key) > 1)
            stop("attempt to select more than one element")
        if (is.na(key))
            stop("subsetting an object of class \"P2Views\" with NA is not supported")
        if (is.numeric(key)) {
            if (!is.integer(key))
                key <- as.integer(key)
            if (key < 1 || length(x) < key)
                stop("subscript out of bounds")
        }
        key
    }
)

setMethod("$", "P2Views", function(x, name) x[[name]])


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "P2ViewsWithoutIDs" class.
### 
### When 'pdict' is a PDict object without pattern IDs then matchPDict(pdict, ...)
### returns the matches in a P2ViewsWithoutIDs object.
###
### Note that in normal operations the user NEVER needs to create a
### P2ViewsWithoutIDs object explicitely or to modify an existing one:
### P2ViewsWithoutIDs objects are created by the matchPDict() function
### and have a read-only semantic.
###
### Slot description:
###
###   ends: a list of integer vectors.
###
###   width: temporary hack. In the future we will probably want to store the
###       starts of the matches when 'pdict' is a PDict other than an
###       ULdna_PDict object. Another solution would be to keep the width slot
###       and to make it the same length as the ends slot (it's currently of
###       length 1 only).
###

setClass("P2ViewsWithoutIDs",
    contains="P2Views",
    representation(
        ends="list",
        width="integer"
    )
)

setMethod("length", "P2ViewsWithoutIDs", function(x) length(x@ends))

setMethod("pids", "P2ViewsWithoutIDs", function(x) NULL)

setMethod("show", "P2ViewsWithoutIDs",
    function(object)
    {
        cat(length(object), "-pattern \"P2Views\" object (no pattern IDs)\n", sep="")
    }
)

setMethod("[[", "P2ViewsWithoutIDs",
    function(x, i, j, ...)
    {
        key <- callNextMethod()
        if (is.character(key))
            stop("\"P2Views\" object has no pattern IDs")
        new("Views", start=x@ends[[key]]-x@width+1L, end=x@ends[[key]], check.data=FALSE)
    }
)

### An example of a P2ViewsWithoutIDs object of length 5 where only the
### 2nd pattern has matches:
###   > ends <- rep(list(integer(0)), 5)
###   > ends[[2]] <- c(199L, 402L)
###   > pm <- new("P2ViewsWithoutIDs", ends=ends, width=10L)
###   > pm[[1]]
###   > pm[[2]]
###   > pm[[6]] # Error in pm[[6]] : subscript out of bounds
###   > p2starts(pm)
###   > p2ends(pm)
###   > p2nmatch(pm)
###
setMethod("p2ends", "P2ViewsWithoutIDs",
    function(x, all.pids=FALSE)
    {
        if (!missing(all.pids))
            warning("'all.pids' is ignored when \"P2Views\" object has no pattern IDs")
        x@ends
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "P2ViewsWithIDs" class.
### 
### When 'pdict' is a PDict object with pattern IDs then matchPDict(pdict, ...)
### returns the matches in a P2ViewsWithIDs object.
###
### Note that in normal operations the user NEVER needs to create a
### P2ViewsWithIDs object explicitely or to modify an existing one:
### P2ViewsWithIDs objects are created by the matchPDict() function
### and have a read-only semantic.
###
### Slot description:
###
###   length: the length of the input dictionary.
###
###   ends: a key-value list (environment) where the values are integer
###       vectors containing the ending positions of the input pattern whose
###       position in the input dictionary is given by the key (the keys are
###       strings representing positive integers).
###       
###   width: temporary hack. In the future we will probably want to store the
###       starts of the matches when 'pdict' is a PDict other than an
###       ULdna_PDict object. Another solution would be to keep the width slot
###       and to make it the same length as the ends slot (it's currently of
###       length 1 only).
###
###   pids: a character vector containing the unique pattern IDs.
###

setClass("P2ViewsWithIDs",
    contains="P2Views",
    representation(
        length="integer",
        ends="environment",
        width="integer",
        pids="character"
    )
)

setMethod("length", "P2ViewsWithIDs", function(x) x@length)

setMethod("pids", "P2ViewsWithIDs", function(x) x@pids)

setMethod("show", "P2ViewsWithIDs",
    function(object)
    {
        cat(length(object), "-pattern \"P2Views\" object\n", sep="")
    }
)

setMethod("[[", "P2ViewsWithIDs",
    function(x, i, j, ...)
    {
        key <- callNextMethod()
        if (is.character(key)) {
            pos <- match(key, pids(x))
            if (is.na(pos))
                stop("pattern ID \"", key, "\" not found")
            key <- pos
        } 
        end <- x@ends[[as.character(key)]]
        if (is.null(end))
            end <- integer(0)
        new("Views", start=end-x@width+1L, end=end, check.data=FALSE)
    }
)

### An example of a P2ViewsWithIDs object of length 5 where only the
### 2nd pattern has matches:
###   > ends <- new.env(hash=TRUE, parent=emptyenv())
###   > ends[["2"]] <- c(199L, 402L)
###   > pm <- new("P2ViewsWithIDs", length=5L, ends=ends, width=10L, pids=letters[1:5])
###   > pm[[1]]
###   > pm[[2]]
###   > pm[[6]] # Error in pm[[6]] : subscript out of bounds
###   > pids(pm)
###   > pm[["a"]]
###   > pm[["b"]]
###   > pm[["aa"]] # Error in pm[["aa"]] : pattern ID ‘aa’ not found
###   > p2starts(pm)
###   > p2starts(pm, all.pids=TRUE)
###   > p2ends(pm)
###   > p2ends(pm, all.pids=TRUE)
###   > p2nmatch(pm)
###   > p2nmatch(pm, all.pids=TRUE)
###
setMethod("p2ends", "P2ViewsWithIDs",
    function(x, all.pids=FALSE)
    {
        if (all.pids)
            ii <- seq_len(length(x))
        else
            ii <- sort(as.integer(ls(x@ends, all.names=TRUE)))
        ans <- lapply(ii, function(i) end(x[[i]]))
        names(ans) <- pids(x)[ii]
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "extractAllMatches" function.
###

extractAllMatches <- function(subject, pdictmatches)
{
    if (!is(subject, "BString"))
        stop("'subject' must be a BString object")
    if (!is(pdictmatches, "P2Views"))
        stop("'pdictmatches' must be a P2Views object")
    if (is.null(pids(pdictmatches)))
        stop("extractAllMatches() works only with a \"P2Views\" object that has pattern IDs")
    starts <- p2starts(pdictmatches)
    if (length(starts) == 0)
        start <- integer(0)
    else
        start <- unlist(starts, recursive=FALSE, use.names=FALSE)
    ends <- p2ends(pdictmatches)
    if (length(ends) == 0)
        end <- integer(0)
    else
        end <- unlist(ends, recursive=FALSE, use.names=FALSE)
    new("BStringViews", subject=subject, start=start, end=end,
        desc=rep(names(ends), p2nmatch(pdictmatches)), check.views=FALSE)
}


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
### The number of integers needed to represent a tree node is the number of
### letters in the alphabet plus 3 (i.e. 7 for DNA input) hence this internal
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
###   > flinks0 # no failure link is set yet
###   > ends <- p2ends(matchPDict(pdict, DNAString("acaagagagt")))
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
###   dups: an integer vector of length L and containing the duplicate info.
###       If unique pattern IDs were provided as part of the input dictionary,
###       then they are used to name the elements of this vector.

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

setMethod("pids", "ULdna_PDict", function(x) names(x@dups))

setMethod("show", "ULdna_PDict",
    function(object)
    {
        cat(length(object), "-pattern uniform-length \"PDict\" object of width ", width(object), sep="")
        if (is.null(pids(object)))
            cat(" (no pattern IDs)")
        cat("\n")
    }
)

debug_ULdna <- function()
{
    invisible(.Call("match_ULdna_debug", PACKAGE="Biostrings"))
}

.checkpids <- function(pids)
{
    if (!is.null(pids)) {
        if (any(pids %in% c("", NA)))
            stop("'dict' has invalid names")
        if (any(duplicated(pids)))
            stop("'dict' has duplicated names")
    }
}

### 'ppans' is the answer returned by the preprocessing functions (the
### ULdna_pp_*() functions).
.ULdna_PDict.postinit <- function(.Object, length, ppans, pids)
{
    .Object@length <- length
    .Object@width <- ppans$width
    .Object@actree <- new("ACtree", ppans$actree_nodes_xp, ppans$actree_base_codes)
    .Object@dups <- ppans$dups
    names(.Object@dups) <- pids
    .Object
}

### 'dict' must be a string vector (aka character vector) with at least 1
### element. If it has names, then they must be unique and are considered
### to be the pattern IDs.
.ULdna_PDict.init_with_StrVect <- function(.Object, dict)
{
    if (any(is.na(dict)))
        stop("'dict' contains NAs")
    pids <- names(dict)
    .checkpids(pids)
    ppans <- .Call("ULdna_pp_StrVect", dict, PACKAGE="Biostrings")
    .ULdna_PDict.postinit(.Object, length(dict), ppans, pids)
}

### 'dict' must be a list of BString objects of the same class (i.e. all
### BString instances or all DNAString instances or etc...) with at least 1
### element.
.ULdna_PDict.init_with_BStringList <- function(.Object, dict)
{
    pids <- names(dict)
    .checkpids(pids)
    dict0 <- lapply(dict, function(pattern) list(pattern@data@xp, pattern@offset, pattern@length))
    ppans <- .Call("ULdna_pp_BStringList", dict0, PACKAGE="Biostrings")
    .ULdna_PDict.postinit(.Object, length(dict), ppans, pids)
}

### 'dict' must be a BStringViews object.
.ULdna_PDict.init_with_BStringViews <- function(.Object, dict)
{
    if (length(dict) == 0)
        stop("'dict' has no views")
    pids <- desc(dict)
    .checkpids(pids)
    ppans <- .Call("ULdna_pp_views",
                  subject(dict)@data@xp, subject(dict)@offset, subject(dict)@length,
                  start(dict), end(dict),
                  PACKAGE="Biostrings")
    .ULdna_PDict.postinit(.Object, length(dict), ppans, pids)
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
###   > system.time(ends <- p2ends(matchPDict(pdict, chr1)))
###      user  system elapsed 
###    50.663   0.000  50.763
###   > nmatches <- sapply(ends, length)
###   > table(nmatches)
###   > id0 <- which(nmatches == max(nmatches))
###   > p0 <- DNAString(dict[id0])
###   > p0
###     25-letter "DNAString" instance
###   Value: CTGTAATCCCAGCACTTTGGGAGGC
###   > subBString(chr1, ends[[id0]][1]-24, ends[[id0]][1]) == p0
###   [1] TRUE
### For a more extensive validation:
###   > pidOK <- sapply(seq_len(length(ends)),
###                     function(pid) identical(ends[[pid]],
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
    pids <- pids(pdict)
    if (is.null(pids))
        envir <- NULL
    else
        envir <- new.env(hash=TRUE, parent=emptyenv())
    ends <- .Call("match_ULdna_exact",
          pdict@length, pdict@dups,
          actree@nodes@xp, actree@base_codes,
          subject@data@xp, subject@offset, subject@length,
          envir,
          PACKAGE="Biostrings")
    if (is.null(pids))
        new("P2ViewsWithoutIDs", ends=ends, width=width(pdict))
    else
        new("P2ViewsWithIDs", length=length(pdict), ends=ends, width=width(pdict), pids=pids)
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
###      > system.time(ends <- p2ends(matchPDict(pdict, chr1)))
###         user  system elapsed
###      105.239   0.188 105.429
###      > nmatches <- sapply(ends, length)
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

setGeneric("matchPDict", signature="subject",
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

