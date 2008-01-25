#############################################################################
###
### The matchPDict() generic & related classes and functions
### ========================================================
###
### Author: Herve Pages (hpages@fhcrc.org)
###
#############################################################################




### =========================================================================
### A. PREPROCESSING
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "PDict" class.
###
### Very general container for any kind of preprocessed pattern dictionary.
###

setClass("PDict", representation("VIRTUAL"))

setGeneric("pids", function(x) standardGeneric("pids"))


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

.ACtree.ints_per_acnode <- function(x) (length(x@base_codes) + 4L)

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
###   > p2end <- p2end(matchPDict(pdict, DNAString("acaagagagt")))
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
        colnames(ans) <- c("parent_id", "depth", pdict@actree@base_codes,
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
        dups="integer",
        stats="list"
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
    .Object@stats <- ppans$stats
    .Object
}

### 'dict' must be a string vector (aka character vector) with at least 1
### element. If it has names, then they must be unique and are considered
### to be the pattern IDs.
.ULdna_PDict.init_with_StrVect <- function(.Object, dict, width)
{
    if (any(is.na(dict)))
        stop("'dict' contains NAs")
    pids <- names(dict)
    .checkpids(pids)
    ppans <- .Call("ULdna_pp_StrVect", dict, width, PACKAGE="Biostrings")
    .ULdna_PDict.postinit(.Object, length(dict), ppans, pids)
}

### 'dict' must be a list of BString objects of the same class (i.e. all
### BString instances or all DNAString instances or etc...) with at least 1
### element.
.ULdna_PDict.init_with_BStringList <- function(.Object, dict, width)
{
    pids <- names(dict)
    .checkpids(pids)
    dict0 <- lapply(dict, function(pattern) list(pattern@data@xp, pattern@offset, pattern@length))
    ppans <- .Call("ULdna_pp_BStringList", dict0, width, PACKAGE="Biostrings")
    .ULdna_PDict.postinit(.Object, length(dict), ppans, pids)
}

### 'dict' must be a BStringViews object.
.ULdna_PDict.init_with_BStringViews <- function(.Object, dict, width)
{
    if (length(dict) == 0)
        stop("'dict' has no views")
    pids <- desc(dict)
    .checkpids(pids)
    ppans <- .Call("ULdna_pp_views",
                  subject(dict)@data@xp, subject(dict)@offset, subject(dict)@length,
                  start(dict), end(dict), width,
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
    function(.Object, dict, width=NULL)
    {
        on.exit(.Call("ULdna_free_actree_nodes_buf", PACKAGE="Biostrings"))
        if (!is.null(width)) {
            if (!is.numeric(width) || length(width) != 1 || is.na(width))
                stop("when specified, 'width' must be a single integer")
            width <- as.integer(width)
            if (width <= 0)
                stop("when specified, 'width' must be a positive integer")
        }
        if (is.character(dict)) {
            if (length(dict) == 0)
                stop("'dict' is an empty character vector")
            return(.ULdna_PDict.init_with_StrVect(.Object, dict, width))
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
                return(.ULdna_PDict.init_with_StrVect(.Object, dict0, width))
            }
            if (extends(pattern_class, "BString"))
                return(.ULdna_PDict.init_with_BStringList(.Object, dict, width))
        }
        if (is(dict, "BStringViews"))
            return(.ULdna_PDict.init_with_BStringViews(.Object, dict, width))
        stop("invalid 'dict' (type '?ULdna_PDict' for more information)")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "TailedULdna_PDict" class.
###

setClass("TailedULdna_PDict",
    contains="ULdna_PDict",
    representation(
        tail="BStringViews"
    )
)

setMethod("show", "TailedULdna_PDict",
    function(object)
    {
        cat(length(object), "-pattern \"PDict\" object with a head (of width ",
            width(object), ") and a tail", sep="")
        if (is.null(pids(object)))
            cat(" (no pattern IDs)")
        cat("\n")
    }
)

setMethod("initialize", "TailedULdna_PDict",
    function(.Object, dict, width)
    {
        if (is.null(width))
            stop("'width' must be a single integer")
        .Object <- callNextMethod(.Object, dict, width=width)
        if (.Object@stats$min.width == .Object@width)
            stop("tails with empty strings are not supported yet, ",
                 "try again with 'width' decreased by 1")
        if (is.character(dict)) {
            .Object@tail <- BStringViews(substr(dict, .Object@width + 1, .Object@stats$max.width),
                                         subjectClass="DNAString")
	} else if (is.list(dict)) {
            stop("'dict' of type list is not yet supported")
        } else if (is(dict, "BStringViews")) {
            dict@views$start <- dict@views$start + .Object@width
            .Object@tail <- dict
        }
        .Object
    }
)




### =========================================================================
### B. SEARCH RESULT MANIPULATION
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "P2Views" API.
### 
### An API for manipulating the match container returned by matchPDict().
###

setClass("P2Views", representation("VIRTUAL"))


setGeneric("p2start", signature="x",
    function(x, all.pids=FALSE) standardGeneric("p2start"))

setGeneric("p2end", signature="x",
    function(x, all.pids=FALSE) standardGeneric("p2end"))

setGeneric("p2nview", signature="x",
    function(x, all.pids=FALSE) standardGeneric("p2nview"))

setMethod("p2nview", "P2Views",
    function(x, all.pids=FALSE)
    {
        if (!is.null(pids(x))) {
            p2end <- p2end(x, all.pids=all.pids)
            if (length(p2end) == 0)
                return(integer(0))
            return(sapply(p2end, length))
        }
        if (!missing(all.pids))
            warning("'all.pids' is ignored when \"P2Views\" object has no pattern IDs")
        p2end <- p2end(x)
        if (length(p2end) == 0)
            return(integer(0))
        sapply(p2end, length)
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
###   > p2v <- new("P2ViewsWithoutIDs", ends=ends, width=10L)
###   > p2v[[1]]
###   > p2v[[2]]
###   > p2v[[6]] # Error in p2v[[6]] : subscript out of bounds
###   > p2start(p2v)
###   > p2end(p2v)
###   > p2nview(p2v)
###
setMethod("p2start", "P2ViewsWithoutIDs",
    function(x, all.pids=FALSE)
    {
        if (!missing(all.pids))
            warning("'all.pids' is ignored when \"P2Views\" object has no pattern IDs")
        .Call("shiftListOfInts", x@ends, 1L - x@width, PACKAGE="Biostrings")
    }
)
setMethod("p2end", "P2ViewsWithoutIDs",
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
###   ends_envir: a key-value list (environment) where the values are integer
###       vectors containing the ending positions of the input pattern whose
###       position in the input dictionary is given by the key (the keys are
###       strings representing positive integers).
###       
###   width: temporary hack. In the future we will probably want to store the
###       starts of the matches when 'pdict' is a PDict other than an
###       ULdna_PDict object. Another solution would be to keep the width slot
###       and to make it the same length as the ends_envir slot (it's currently
###       of length 1 only).
###
###   pids: a character vector containing the unique pattern IDs.
###

setClass("P2ViewsWithIDs",
    contains="P2Views",
    representation(
        length="integer",
        ends_envir="environment",
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
        end <- x@ends_envir[[formatC(key, width=10, format="d", flag="0")]]
        if (is.null(end))
            end <- integer(0)
        new("Views", start=end-x@width+1L, end=end, check.data=FALSE)
    }
)

### An example of a P2ViewsWithIDs object of length 5 where only the
### 2nd pattern has matches:
###   > ends_envir <- new.env(hash=TRUE, parent=emptyenv())
###   > ends_envir[['0000000002']] <- c(199L, 402L)
###   > p2v <- new("P2ViewsWithIDs", length=5L, ends_envir=ends_envir, width=10L, pids=letters[1:5])
###   > p2v[[1]]
###   > p2v[[2]]
###   > p2v[[6]] # Error in p2v[[6]] : subscript out of bounds
###   > pids(p2v)
###   > p2v[["a"]]
###   > p2v[["b"]]
###   > p2v[["aa"]] # Error in p2v[["aa"]] : pattern ID ‘aa’ not found
###   > p2start(p2v)
###   > p2start(p2v, all.pids=TRUE)
###   > p2end(p2v)
###   > p2end(p2v, all.pids=TRUE)
###   > p2nview(p2v)
###   > p2nview(p2v, all.pids=TRUE)
###
setMethod("p2start", "P2ViewsWithIDs",
    function(x, all.pids=FALSE)
    {
        if (!is.logical(all.pids) || length(all.pids) != 1 || is.na(all.pids))
            stop("'all.pids' must be 'TRUE' or 'FALSE'")
        .Call("extract_p2end", x@ends_envir, 1L - x@width, x@pids, all.pids, PACKAGE="Biostrings")
    }
)
setMethod("p2end", "P2ViewsWithIDs",
    function(x, all.pids=FALSE)
    {
        if (!is.logical(all.pids) || length(all.pids) != 1 || is.na(all.pids))
            stop("'all.pids' must be 'TRUE' or 'FALSE'")
        .Call("extract_p2end", x@ends_envir, 0L, x@pids, all.pids, PACKAGE="Biostrings")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "unlist" method and the "extractAllMatches" function.
###

setMethod("unlist", "P2Views",
    function(x, recursive=TRUE, use.names=TRUE)
    {
        p2start <- p2start(x)
        if (length(p2start) == 0)
            start <- integer(0)
        else
            start <- unlist(p2start, recursive=FALSE, use.names=FALSE)
        p2end <- p2end(x)
        if (length(p2end) == 0)
            end <- integer(0)
        else
            end <- unlist(p2end, recursive=FALSE, use.names=FALSE)
        desc <- names(p2end)
        if (!is.null(desc))
            desc <- rep.int(desc, times=sapply(p2end, length))
        new("Views", start=start, end=end, desc=desc, check.data=FALSE)
    }
)

extractAllMatches <- function(subject, p2v)
{
    if (!is(subject, "BString"))
        stop("'subject' must be a BString object")
    if (!is(p2v, "P2Views"))
        stop("'p2v' must be a P2Views object")
    if (is.null(pids(p2v)))
        stop("extractAllMatches() works only with a \"P2Views\" object that has pattern IDs")
    allviews <- unlist(p2v)
    new("BStringViews", subject=subject, start=start(allviews), end=end(allviews),
        desc=desc(allviews), check.views=FALSE)
}




### =========================================================================
### C. EXACT SEARCH
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Aho-Corasick
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
###   > system.time(p2end <- p2end(matchPDict(pdict, chr1)))
###      user  system elapsed 
###    50.663   0.000  50.763
###   > nmatches <- sapply(p2end, length)
###   > table(nmatches)
###   > id0 <- which(nmatches == max(nmatches))
###   > p0 <- DNAString(dict[id0])
###   > p0
###     25-letter "DNAString" instance
###   Value: CTGTAATCCCAGCACTTTGGGAGGC
###   > subBString(chr1, p2end[[id0]][1]-24, p2end[[id0]][1]) == p0
###   [1] TRUE
### For a more extensive validation:
###   > pidOK <- sapply(seq_len(length(p2end)),
###                     function(pid) identical(p2end[[pid]],
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

.match.ULdna_PDict.exact <- function(pdict, subject, count.only)
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
          envir, count.only,
          PACKAGE="Biostrings")
    if (is.null(pids))
        new("P2ViewsWithoutIDs", ends=ends, width=width(pdict))
    else
        new("P2ViewsWithIDs", length=length(pdict), ends_envir=ends, width=width(pdict), pids=pids)
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
###      > system.time(p2end <- p2end(matchPDict(pdict, chr1)))
###         user  system elapsed
###      105.239   0.188 105.429
###      > nmatches <- sapply(p2end, length)
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




### =========================================================================
### D. THE "matchPDict" AND "countPDict" GENERIC FUNCTION AND METHODS
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .matchPDict()
###

.matchPDict <- function(pdict, subject, algorithm, max.mismatch, fixed, count.only=FALSE)
{
    if (!is(pdict, "ULdna_PDict"))
        stop("the pattern dictionary 'pdict' can only be a ULdna_PDict object for now")
    if (!is(subject, "DNAString"))
        stop("'subject' can only be a DNAString object for now")
    if (!identical(algorithm, "auto"))
        stop("'algo' can only be '\"auto\"' for now")
    if (!identical(max.mismatch, 0))
        stop("'max.mismatch' can only be set to 0 for now")
    if (!identical(fixed, TRUE))
        stop("'fixed' can only be 'TRUE' for now")
    .match.ULdna_PDict.exact(pdict, subject, count.only)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "matchPDict" generic and methods.
###

setGeneric("matchPDict", signature="subject",
    function(pdict, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)
        standardGeneric("matchPDict")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPDict", "BString",
    function(pdict, subject, algorithm, max.mismatch, fixed)
    {
	.matchPDict(pdict, subject, algorithm, max.mismatch, fixed)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "countPDict" generic and methods.
###

setGeneric("countPDict", signature="subject",
    function(pdict, subject, algorithm="auto", max.mismatch=0, fixed=TRUE)
        standardGeneric("countPDict")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPDict", "BString",
    function(pdict, subject, algorithm, max.mismatch, fixed)
    {
	.matchPDict(pdict, subject, algorithm, max.mismatch, fixed, count.only=TRUE)
    }
)

