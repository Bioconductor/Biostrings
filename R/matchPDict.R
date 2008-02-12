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

setGeneric("patternFrequency", function(x) standardGeneric("patternFrequency"))


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
###   > end_index <- endIndex(matchPDict(pdict, DNAString("acaagagagt")))
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
###       If unique pattern names were provided as part of the input dictionary,
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

setMethod("names", "ULdna_PDict", function(x) names(x@dups))

setMethod("show", "ULdna_PDict",
    function(object)
    {
        cat(length(object), "-pattern uniform-length \"PDict\" object of width ", width(object), sep="")
        if (is.null(names(object)))
            cat(" (patterns have no names)")
        cat("\n")
    }
)

setMethod("duplicated", "ULdna_PDict",
    function(x, incomparables=FALSE, ...)
        x@dups != 0
)

setMethod("patternFrequency", "ULdna_PDict",
    function(x)
    {
        ans <- rep.int(1L, x@length)
        tb <- table(x@dups[x@dups != 0])
        for (pos in names(tb)) {
            i <- as.integer(pos)
            ans[c(i, which(x@dups == i))] <- 1L + tb[pos][[1]]
        }
        ans
    }
)

debug_ULdna <- function()
{
    invisible(.Call("match_TPdna_debug", PACKAGE="Biostrings"))
}

.check.names <- function(names)
{
    if (!is.null(names)) {
        if (any(names %in% c("", NA)))
            stop("'dict' has invalid names")
        if (any(duplicated(names)))
            stop("'dict' has duplicated names")
    }
}

### 'ppans' is the answer returned by the preprocessing functions (the
### ULdna_pp_*() functions).
.ULdna_PDict.postinit <- function(.Object, length, ppans, names)
{
    .Object@length <- length
    .Object@width <- ppans$width
    .Object@actree <- new("ACtree", ppans$actree_nodes_xp, ppans$actree_base_codes)
    .Object@dups <- ppans$dups
    names(.Object@dups) <- names
    .Object@stats <- ppans$stats
    .Object
}

### 'dict' must be a string vector (aka character vector) with at least 1
### element. If it has names, then they must be unique and are considered
### to be the pattern names.
.ULdna_PDict.init_with_StrVect <- function(.Object, dict, width)
{
    if (any(is.na(dict)))
        stop("'dict' contains NAs")
    names <- names(dict)
    .check.names(names)
    ppans <- .Call("ULdna_pp_StrVect", dict, width, PACKAGE="Biostrings")
    .ULdna_PDict.postinit(.Object, length(dict), ppans, names)
}

### 'dict' must be a list of BString objects of the same class (i.e. all
### BString instances or all DNAString instances or etc...) with at least 1
### element.
.ULdna_PDict.init_with_BStringList <- function(.Object, dict, width)
{
    names <- names(dict)
    .check.names(names)
    dict0 <- lapply(dict, function(pattern) list(pattern@data@xp, pattern@offset, pattern@length))
    ppans <- .Call("ULdna_pp_BStringList", dict0, width, PACKAGE="Biostrings")
    .ULdna_PDict.postinit(.Object, length(dict), ppans, names)
}

### 'dict' must be a BStringViews object.
.ULdna_PDict.init_with_BStringViews <- function(.Object, dict, width)
{
    if (length(dict) == 0)
        stop("'dict' has no views")
    names <- desc(dict)
    .check.names(names)
    ppans <- .Call("ULdna_pp_views",
                  subject(dict)@data@xp, subject(dict)@offset, subject(dict)@length,
                  start(dict), end(dict), width,
                  PACKAGE="Biostrings")
    .ULdna_PDict.postinit(.Object, length(dict), ppans, names)
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
            if (!isSingleNumber(width))
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
### The "TPdna_PDict" class.
###
### A container for storing a preprocessed "Trusted Prefix" dictionary (or
### set) of DNA patterns.
###

setClass("TPdna_PDict",
    contains="ULdna_PDict",
    representation(
        tail="DNAStringList"
    )
)

setMethod("tail", "TPdna_PDict", function(x, ...) x@tail)

setMethod("show", "TPdna_PDict",
    function(object)
    {
        cat(length(object), "-pattern \"PDict\" object (of width ",
            width(object), ") with a tail", sep="")
        if (is.null(names(object)))
            cat(" (patterns have no names)")
        cat("\n")
    }
)

setMethod("duplicated", "TPdna_PDict",
    function(x, incomparables=FALSE, ...)
    {
        stop("not ready yet, sorry!")
    }
)

setMethod("patternFrequency", "TPdna_PDict",
    function(x)
    {
        stop("not ready yet, sorry!")
    }
)

setMethod("initialize", "TPdna_PDict",
    function(.Object, dict, width)
    {
        if (is.null(width))
            stop("'width' must be a single integer")
        .Object <- callNextMethod(.Object, dict, width=width)
        .Object@tail <- DNAStringList(dict, start=.Object@width+1L)
        .Object
    }
)




### =========================================================================
### B. SEARCH RESULT MANIPULATION
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "ViewsIndex" API.
### 
### An API for manipulating the match container returned by matchPDict().
###

setClass("ViewsIndex", representation("VIRTUAL"))


setGeneric("startIndex", signature="x",
    function(x, all.names=FALSE) standardGeneric("startIndex"))

setGeneric("endIndex", signature="x",
    function(x, all.names=FALSE) standardGeneric("endIndex"))

setGeneric("countIndex", signature="x",
    function(x, all.names=FALSE) standardGeneric("countIndex"))

setMethod("countIndex", "ViewsIndex",
    function(x, all.names=FALSE)
    {
        if (!is.null(names(x))) {
            end_index <- endIndex(x, all.names=all.names)
            if (length(end_index) == 0)
                return(integer(0))
            return(sapply(end_index, length))
        }
        if (!missing(all.names))
            warning("'all.names' is ignored when patterns have no names")
        end_index <- endIndex(x)
        if (length(end_index) == 0)
            return(integer(0))
        sapply(end_index, length)
    }
)

### Return a single integer or string (not NA).
setMethod("[[", "ViewsIndex",
    function(x, i, j, ...)
    {
        # 'x' is guaranteed to be a "ViewsIndex" object (if it's not, then
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
            stop("wrong argument for subsetting an object of class \"ViewsIndex\"")
        if (length(key) < 1)
            stop("attempt to select less than one element")
        if (length(key) > 1)
            stop("attempt to select more than one element")
        if (is.na(key))
            stop("subsetting an object of class \"ViewsIndex\" with NA is not supported")
        if (is.numeric(key)) {
            if (!is.integer(key))
                key <- as.integer(key)
            if (key < 1 || length(x) < key)
                stop("subscript out of bounds")
        }
        key
    }
)

setMethod("$", "ViewsIndex", function(x, name) x[[name]])


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "ByPos_ViewsIndex" class.
### 
### When 'pdict' is a PDict object with no names then matchPDict(pdict, ...)
### returns the matches in a ByPos_ViewsIndex object.
###
### Note that in normal operations the user NEVER needs to create a
### ByPos_ViewsIndex object explicitely or to modify an existing one:
### ByPos_ViewsIndex objects are created by the matchPDict() function
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

setClass("ByPos_ViewsIndex",
    contains="ViewsIndex",
    representation(
        ends="list",
        width="integer"
    )
)

setMethod("length", "ByPos_ViewsIndex", function(x) length(x@ends))

setMethod("names", "ByPos_ViewsIndex", function(x) NULL)

setMethod("show", "ByPos_ViewsIndex",
    function(object)
    {
        cat(length(object), "-pattern \"ViewsIndex\" object (patterns have no names)\n", sep="")
    }
)

setMethod("[[", "ByPos_ViewsIndex",
    function(x, i, j, ...)
    {
        key <- callNextMethod()
        if (is.character(key))
            stop("\"ViewsIndex\" object has no names")
        new("Views", start=x@ends[[key]]-x@width+1L, end=x@ends[[key]], check.data=FALSE)
    }
)

### An example of a ByPos_ViewsIndex object of length 5 where only the
### 2nd pattern has matches:
###   > ends <- rep(list(integer(0)), 5)
###   > ends[[2]] <- c(199L, 402L)
###   > vindex <- new("ByPos_ViewsIndex", ends=ends, width=10L)
###   > vindex[[1]]
###   > vindex[[2]]
###   > vindex[[6]] # Error in vindex[[6]] : subscript out of bounds
###   > startIndex(vindex)
###   > endIndex(vindex)
###   > countIndex(vindex)
###
setMethod("startIndex", "ByPos_ViewsIndex",
    function(x, all.names=FALSE)
    {
        if (!missing(all.names))
            warning("'all.names' is ignored when patterns have no names")
        .Call("shiftListOfInts", x@ends, 1L - x@width, PACKAGE="Biostrings")
    }
)
setMethod("endIndex", "ByPos_ViewsIndex",
    function(x, all.names=FALSE)
    {
        if (!missing(all.names))
            warning("'all.names' is ignored when patterns have no names")
        x@ends
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "ByName_ViewsIndex" class.
### 
### When 'pdict' is a PDict object with no names then matchPDict(pdict, ...)
### returns the matches in a ByName_ViewsIndex object.
###
### Note that in normal operations the user NEVER needs to create a
### ByName_ViewsIndex object explicitely or to modify an existing one:
### ByName_ViewsIndex objects are created by the matchPDict() function
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
###   names: a character vector containing the _unique_ pattern names.
###

setClass("ByName_ViewsIndex",
    contains="ViewsIndex",
    representation(
        length="integer",
        ends_envir="environment",
        width="integer",
        NAMES="character" # R doesn't like @names !!
    )
)

setMethod("length", "ByName_ViewsIndex", function(x) x@length)

setMethod("names", "ByName_ViewsIndex", function(x) x@NAMES)

setMethod("show", "ByName_ViewsIndex",
    function(object)
    {
        cat(length(object), "-pattern \"ViewsIndex\" object\n", sep="")
    }
)

setMethod("[[", "ByName_ViewsIndex",
    function(x, i, j, ...)
    {
        key <- callNextMethod()
        if (is.character(key)) {
            pos <- match(key, names(x))
            if (is.na(pos))
                stop("pattern name \"", key, "\" not found")
            key <- pos
        } 
        end <- x@ends_envir[[formatC(key, width=10, format="d", flag="0")]]
        if (is.null(end))
            end <- integer(0)
        new("Views", start=end-x@width+1L, end=end, check.data=FALSE)
    }
)

### An example of a ByName_ViewsIndex object of length 5 where only the
### 2nd pattern has matches:
###   > ends_envir <- new.env(hash=TRUE, parent=emptyenv())
###   > ends_envir[['0000000002']] <- c(199L, 402L)
###   > vindex <- new("ByName_ViewsIndex", length=5L, ends_envir=ends_envir, width=10L, NAMES=letters[1:5])
###   > vindex[[1]]
###   > vindex[[2]]
###   > vindex[[6]] # Error in vindex[[6]] : subscript out of bounds
###   > names(vindex)
###   > vindex[["a"]]
###   > vindex[["b"]]
###   > vindex[["aa"]] # Error in vindex[["aa"]] : pattern name ‘aa’ not found
###   > startIndex(vindex)
###   > startIndex(vindex, all.names=TRUE)
###   > endIndex(vindex)
###   > endIndex(vindex, all.names=TRUE)
###   > countIndex(vindex)
###   > countIndex(vindex, all.names=TRUE)
###
setMethod("startIndex", "ByName_ViewsIndex",
    function(x, all.names=FALSE)
    {
        if (!isTRUEorFALSE(all.names))
            stop("'all.names' must be 'TRUE' or 'FALSE'")
        .Call("extract_endIndex", x@ends_envir, 1L - x@width, x@NAMES, all.names, PACKAGE="Biostrings")
    }
)
setMethod("endIndex", "ByName_ViewsIndex",
    function(x, all.names=FALSE)
    {
        if (!isTRUEorFALSE(all.names))
            stop("'all.names' must be 'TRUE' or 'FALSE'")
        .Call("extract_endIndex", x@ends_envir, 0L, x@NAMES, all.names, PACKAGE="Biostrings")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "unlist" method and the "extractAllMatches" function.
###

setMethod("unlist", "ViewsIndex",
    function(x, recursive=TRUE, use.names=TRUE)
    {
        start_index <- startIndex(x)
        if (length(start_index) == 0)
            start <- integer(0)
        else
            start <- unlist(start_index, recursive=FALSE, use.names=FALSE)
        end_index <- endIndex(x)
        if (length(end_index) == 0)
            end <- integer(0)
        else
            end <- unlist(end_index, recursive=FALSE, use.names=FALSE)
        desc <- names(end_index)
        if (!is.null(desc))
            desc <- rep.int(desc, times=sapply(end_index, length))
        new("Views", start=start, end=end, desc=desc, check.data=FALSE)
    }
)

extractAllMatches <- function(subject, vindex)
{
    if (!is(subject, "BString"))
        stop("'subject' must be a BString object")
    if (!is(vindex, "ViewsIndex"))
        stop("'vindex' must be a ViewsIndex object")
    if (is.null(names(vindex)))
        stop("extractAllMatches() works only with a \"ViewsIndex\" object with names")
    allviews <- unlist(vindex)
    new("BStringViews", subject=subject, start=start(allviews), end=end(allviews),
        desc=desc(allviews), check.views=FALSE)
}




### =========================================================================
### C. EXACT MATCHING
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
###   > system.time(end_index <- endIndex(matchPDict(pdict, chr1)))
###      user  system elapsed 
###    50.663   0.000  50.763
###   > count_index <- sapply(end_index, length)
###   > table(count_index)
###   > id0 <- which(count_index == max(count_index))
###   > p0 <- DNAString(dict[id0])
###   > p0
###     25-letter "DNAString" instance
###   Value: CTGTAATCCCAGCACTTTGGGAGGC
###   > subBString(chr1, end_index[[id0]][1]-24, end_index[[id0]][1]) == p0
###   [1] TRUE
### For a more extensive validation:
###   > pidOK <- sapply(seq_len(length(end_index)),
###                     function(pid) identical(end_index[[pid]],
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
    names <- names(pdict)
    if (is.null(names))
        envir <- NULL
    else
        envir <- new.env(hash=TRUE, parent=emptyenv())
    ans <- .Call("match_TPdna",
                 actree@nodes@xp, actree@base_codes,
                 pdict@dups, NULL,
                 subject,
                 0L, c(TRUE, TRUE),
                 count.only, envir,
                 PACKAGE="Biostrings")
    if (count.only)
        return(ans)
    if (is.null(names))
        new("ByPos_ViewsIndex", ends=ans, width=width(pdict))
    else
        new("ByName_ViewsIndex", length=length(pdict), ends_envir=ans, width=width(pdict), NAMES=names)
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
###      > system.time(end_index <- endIndex(matchPDict(pdict, chr1)))
###         user  system elapsed
###      105.239   0.188 105.429
###      > count_index <- sapply(end_index, length)
###      > sum(count_index) # most likely no match were found
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
### D. INEXACT MATCHING
### -------------------------------------------------------------------------

### Example:
###
###   pdict <- new("TPdna_PDict", c("acgt", "gt", "cgt", "ac"), 2)
###   endIndex(matchPDict(pdict, DNAString("acggaccg"), max.mismatch=0))
###   endIndex(matchPDict(pdict, DNAString("acggaccg"), max.mismatch=1))
###   endIndex(matchPDict(pdict, DNAString("acggaccg"), max.mismatch=2))
###
### At the moment, preprocessing a big TPdna_PDict object is very slow:
###   > library(drosophila2probe)
###   > dict0 <- drosophila2probe$sequence 
###   > system.time(pdict0 <- new("ULdna_PDict", dict0[1:40000]))
###      user  system elapsed
###     0.040   0.032   0.072
###   > system.time(pdict <- new("TPdna_PDict", dict0[1:40000], 10))
###      user  system elapsed
###    38.158   0.052  38.217
###
###   > library(BSgenome.Dmelanogaster.FlyBase.r51)
###   > chr3R <- Dmelanogaster[["3R"]]
###   > system.time(vindex0 <- matchPDict(pdict0, chr3R))
###      user  system elapsed
###     1.352   0.000   1.352
###   > system.time(vindex <- matchPDict(pdict, chr3R))
###      user  system elapsed
###     1.332   0.008   1.338
###   > identical(countIndex(vindex0), countIndex(vindex))
###   [1] TRUE
###
### Allowing mismatches is fast:
###   > system.time(vindex_mm6 <- matchPDict(pdict, chr3R, max.mismatch=4))
###      user  system elapsed
###     1.377   0.000   1.375
###   > vindex_mm6[[103]]
###        start      end width
###   1  9381276  9381285    10
###   2 16070100 16070109    10
###   > v <- views(chr3R, start(vindex_mm6[[103]]), end(vindex_mm6[[103]])+15)
###   > mismatch(dict0[103], v)
###   [[1]]
###   [1] 14 15 19 23 24 25
###
###   [[2]]
###   integer(0)

.match.TPdna_PDict <- function(pdict, subject, max.mismatch, fixed, count.only)
{
    actree <- .ACtree.prepare_for_use_on_DNAString(pdict@actree)
    names <- names(pdict)
    if (is.null(names))
        envir <- NULL
    else
        envir <- new.env(hash=TRUE, parent=emptyenv())
    ans <- .Call("match_TPdna",
                 actree@nodes@xp, actree@base_codes,
                 pdict@dups, pdict@tail@seqs,
                 subject,
                 max.mismatch, fixed,
                 count.only, envir,
                 PACKAGE="Biostrings")
    if (count.only)
        return(ans)
    if (is.null(names))
        new("ByPos_ViewsIndex", ends=ans, width=width(pdict))
    else
        new("ByName_ViewsIndex", length=length(pdict), ends_envir=ans, width=width(pdict), NAMES=names)
}




### =========================================================================
### E. THE "matchPDict" AND "countPDict" GENERIC FUNCTION AND METHODS
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .matchPDict()
###

.matchPDict <- function(pdict, subject, algorithm, max.mismatch, fixed, count.only=FALSE)
{
    if (!is(pdict, "ULdna_PDict"))
        stop("the pattern dictionary 'pdict' can only be a ULdna_PDict (or derived) object for now")
    if (!is(subject, "DNAString"))
        stop("'subject' can only be a DNAString object for now")
    if (!identical(algorithm, "auto"))
        stop("'algo' can only be '\"auto\"' for now")
    max.mismatch <- normalize.max.mismatch(max.mismatch)
    fixed <- normalize.fixed(fixed, class(subject))
    if (is(pdict, "TPdna_PDict"))
        return(.match.TPdna_PDict(pdict, subject, max.mismatch, fixed, count.only))
    if (max.mismatch != 0 || !all(fixed))
        stop("only TPdna_PDict dictionaries support inexact matching")
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

