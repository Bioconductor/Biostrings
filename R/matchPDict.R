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
###   > pdict <- PDict(c("agg", "aca", "gag", "caa", "agt", "aca"))
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
### The "CWdna_PDict" class.
###
### A container for storing a preprocessed constant width dictionary (or set)
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

setClass("CWdna_PDict",
    contains="PDict",
    representation(
        length="integer",
        width="integer",
        actree="ACtree",
        dups="integer",
        stats="list"
    )
)

### Not intended to be used directly by the user.
setMethod("initialize", "CWdna_PDict",
    function(.Object, length, pp_Cans, names)
    {
        .Object@length <- length
        .Object@width <- pp_Cans$width
        .Object@actree <- new("ACtree", pp_Cans$actree_nodes_xp, pp_Cans$actree_base_codes)
        .Object@dups <- pp_Cans$dups
        names(.Object@dups) <- names
        .Object@stats <- pp_Cans$stats
        .Object
    }
)

setMethod("length", "CWdna_PDict", function(x) x@length)

setMethod("width", "CWdna_PDict", function(x) x@width)

setMethod("names", "CWdna_PDict", function(x) names(x@dups))

setMethod("show", "CWdna_PDict",
    function(object)
    {
        cat(length(object), "-pattern constant width PDict object of width ",
            width(object), sep="")
        if (is.null(names(object)))
            cat(" (patterns have no names)")
        cat("\n")
    }
)

setMethod("duplicated", "CWdna_PDict",
    function(x, incomparables=FALSE, ...)
        x@dups != 0
)

setMethod("patternFrequency", "CWdna_PDict",
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "TBdna_PDict" class.
###
### A container for storing a preprocessed "Trusted Band" dictionary (or
### set) of DNA patterns.
### 2 important particular cases of "Trusted Band" dictionaries are:
###   1) "Trusted Prefix": a "Trusted Prefix" dictionary is a "Trusted Band"
###      dictionary with no head.
###   2) "Trusted Suffix": a "Trusted Suffix" dictionary is a "Trusted Band"
###      dictionary with no tail.
###

### 'head' and 'tail' must be DNAStringSet objects of the same length.
### One of them can be NULL but they can't both be NULL at the same time.
setClass("TBdna_PDict",
    contains="CWdna_PDict",
    representation(
        head="DNAStringSet", # can be NULL
        tail="DNAStringSet"  # can be NULL
    )
)

### Not intended to be used directly by the user.
setMethod("initialize", "TBdna_PDict",
    function(.Object, length, pp_Cans, names)
    {
        callNextMethod(.Object, length, pp_Cans, names)
    }
)

setMethod("head", "TBdna_PDict", function(x, ...) x@head)
setMethod("tail", "TBdna_PDict", function(x, ...) x@tail)

setMethod("show", "TBdna_PDict",
    function(object)
    {
        if (is.null(object@head))
            trusted_part <- "Prefix"
        else if (is.null(object@tail))
            trusted_part <- "Suffix"
        else
            trusted_part <- "Band"
        cat(length(object), "-pattern PDict object with a Trusted ",
            trusted_part, " of width ", width(object), sep="")
        if (is.null(names(object)))
            cat(" (patterns have no names)")
        cat("\n")
    }
)

setMethod("duplicated", "TBdna_PDict",
    function(x, incomparables=FALSE, ...)
    {
        stop("not ready yet, sorry!")
    }
)

setMethod("patternFrequency", "TBdna_PDict",
    function(x)
    {
        stop("not ready yet, sorry!")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The PDict() constructor (user-friendly).
###

debug_TBdna <- function()
{
    invisible(.Call("match_TBdna_debug", PACKAGE="Biostrings"))
}

.PDict <- function(dict, names, tb.start, tb.end,
                   drop.head, drop.tail, skip.invalid.patterns)
{
    if (!is.null(names)) {
        if (any(names %in% c("", NA)))
            stop("'dict' has invalid names")
        if (any(duplicated(names)))
            stop("'dict' has duplicated names")
    }
    if (!isSingleNumberOrNA(tb.start))
        stop("'tb.start' must be a single integer or NA")
    if (!is.integer(tb.start))
        tb.start <- as.integer(tb.start)
    if (!isSingleNumberOrNA(tb.end))
        stop("'tb.end' must be a single integer or NA")
    if (!is.integer(tb.end))
        tb.end <- as.integer(tb.end)
    if (!isTRUEorFALSE(drop.head))
        stop("'drop.head' must be 'TRUE' or 'FALSE'")
    if (!isTRUEorFALSE(drop.tail))
        stop("'drop.tail' must be 'TRUE' or 'FALSE'")
    on.exit(.Call("CWdna_free_actree_nodes_buf", PACKAGE="Biostrings"))
    if (is.character(dict)) {
        pp_Cans <- .Call("CWdna_pp_charseqs",
                         dict,
                         tb.start, tb.end,
                         PACKAGE="Biostrings")
    } else if (is(dict, "DNAStringSet")) {
        pp_Cans <- .Call("CWdna_pp_BStringSet",
                         dict,
                         tb.start, tb.end,
                         PACKAGE="Biostrings")
    } else {
        stop("unsuported 'dict' type")
    }
    if (is.na(tb.start) || tb.start == 1L)
        drop.head <- TRUE
    if (is.na(tb.end) || tb.end == -1L)
        drop.tail <- TRUE
    if (drop.head && drop.tail)
        return(new("CWdna_PDict", length(dict), pp_Cans, names))
    ans <- new("TBdna_PDict", length(dict), pp_Cans, names)
    if (!drop.head) 
        ans@head <- DNAStringSet(dict, end=tb.start-1L)
    if (!drop.tail)
        ans@tail <- DNAStringSet(dict, start=tb.end+1L)
    ans
}

### The input dictionary 'dict' must be:
###   - of length >= 1
### The supported types for the input dictionary are:
###   - character vector
###   - DNAStringSet object
###   - BStringViews object
### Typical use:
###   > library(Biostrings)
###   > dict <- c("abb", "aca", "bab", "caa", "abd", "aca")
###   > pdict <- PDict(dict)
###
setGeneric("PDict", signature="dict",
    function(dict, tb.start=1, tb.end=NA,
             drop.head=FALSE, drop.tail=FALSE, skip.invalid.patterns=FALSE)
        standardGeneric("PDict")
)

setMethod("PDict", "character",
    function(dict, tb.start=1, tb.end=NA,
             drop.head=FALSE, drop.tail=FALSE, skip.invalid.patterns=FALSE)
    {
        if (length(dict) == 0)
            stop("'dict' must contain at least one pattern")
        .PDict(dict, names(dict), tb.start, tb.end,
               drop.head, drop.tail, skip.invalid.patterns)
    }
)

setMethod("PDict", "DNAStringSet",
    function(dict, tb.start=1, tb.end=NA,
             drop.head=FALSE, drop.tail=FALSE, skip.invalid.patterns=FALSE)
    {
        if (length(dict) == 0)
            stop("'dict' must contain at least one pattern")
        .PDict(dict, names(dict), tb.start, tb.end,
               drop.head, drop.tail, skip.invalid.patterns)
    }
)

setMethod("PDict", "BStringViews",
    function(dict, tb.start=1, tb.end=NA,
             drop.head=FALSE, drop.tail=FALSE, skip.invalid.patterns=FALSE)
    {
        dict <- DNAStringSet(dict)
        PDict(dict, tb.start=tb.start, tb.end=tb.end,
                    drop.head=drop.head, drop.tail=drop.tail,
                    skip.invalid.patterns=skip.invalid.patterns)
    }
)

### Just because of those silly "AsIs" objects found in the probe packages
### (e.g. drosophila2probe$sequence)
setMethod("PDict", "AsIs",
    function(dict, tb.start=1, tb.end=NA,
             drop.head=FALSE, drop.tail=FALSE, skip.invalid.patterns=FALSE)
    {
        if (!is.character(dict))
            stop("unsuported 'dict' type")
        class(dict) <- "character" # keeps the names (unlike as.character())
        PDict(dict, tb.start=tb.start, tb.end=tb.end,
                    drop.head=drop.head, drop.tail=drop.tail,
                    skip.invalid.patterns=skip.invalid.patterns)
    }
)




### =========================================================================
### B. SEARCH RESULT MANIPULATION
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "MIndex" API.
### 
### An API for manipulating the match container returned by matchPDict().
###

setClass("MIndex", representation("VIRTUAL"))


setGeneric("startIndex", signature="x",
    function(x, all.names=FALSE) standardGeneric("startIndex"))

setGeneric("endIndex", signature="x",
    function(x, all.names=FALSE) standardGeneric("endIndex"))

setGeneric("countIndex", signature="x",
    function(x, all.names=FALSE) standardGeneric("countIndex"))

setMethod("countIndex", "MIndex",
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
setMethod("[[", "MIndex",
    function(x, i, j, ...)
    {
        # 'x' is guaranteed to be a "MIndex" object (if it's not, then
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
            stop("wrong argument for subsetting an object of class \"MIndex\"")
        if (length(key) < 1)
            stop("attempt to select less than one element")
        if (length(key) > 1)
            stop("attempt to select more than one element")
        if (is.na(key))
            stop("subsetting an object of class \"MIndex\" with NA is not supported")
        if (is.numeric(key)) {
            if (!is.integer(key))
                key <- as.integer(key)
            if (key < 1 || length(x) < key)
                stop("subscript out of bounds")
        }
        key
    }
)

setMethod("$", "MIndex", function(x, name) x[[name]])


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "ByPos_MIndex" class.
### 
### When 'pdict' is a PDict object with no names then matchPDict(pdict, ...)
### returns the matches in a ByPos_MIndex object.
###
### Note that in normal operations the user NEVER needs to create a
### ByPos_MIndex object explicitely or to modify an existing one:
### ByPos_MIndex objects are created by the matchPDict() function
### and have a read-only semantic.
###
### Slot description:
###
###   ends: a list of integer vectors.
###
###   width: temporary hack. In the future we will probably want to store the
###       starts of the matches when 'pdict' is a PDict other than an
###       CWdna_PDict object. Another solution would be to keep the width slot
###       and to make it the same length as the ends slot (it's currently of
###       length 1 only).
###

setClass("ByPos_MIndex",
    contains="MIndex",
    representation(
        ends="list",
        width="integer"
    )
)

setMethod("length", "ByPos_MIndex", function(x) length(x@ends))

setMethod("names", "ByPos_MIndex", function(x) NULL)

setMethod("show", "ByPos_MIndex",
    function(object)
    {
        cat(length(object), "-pattern \"MIndex\" object (patterns have no names)\n", sep="")
    }
)

setMethod("[[", "ByPos_MIndex",
    function(x, i, j, ...)
    {
        key <- callNextMethod()
        if (is.character(key))
            stop("\"MIndex\" object has no names")
        ans_end <- x@ends[[key]]
        ans_width <- rep.int(x@width, length(ans_end))
        ans_start <- ans_end - ans_width + 1L
        new("IRanges", start=ans_start, width=ans_width, check=FALSE)
    }
)

### An example of a ByPos_MIndex object of length 5 where only the
### 2nd pattern has matches:
###   > ends <- rep(list(integer(0)), 5)
###   > ends[[2]] <- c(199L, 402L)
###   > mindex <- new("ByPos_MIndex", ends=ends, width=10L)
###   > mindex[[1]]
###   > mindex[[2]]
###   > mindex[[6]] # Error in mindex[[6]] : subscript out of bounds
###   > startIndex(mindex)
###   > endIndex(mindex)
###   > countIndex(mindex)
###
setMethod("startIndex", "ByPos_MIndex",
    function(x, all.names=FALSE)
    {
        if (!missing(all.names))
            warning("'all.names' is ignored when patterns have no names")
        .Call("shiftListOfInts", x@ends, 1L - x@width, PACKAGE="Biostrings")
    }
)
setMethod("endIndex", "ByPos_MIndex",
    function(x, all.names=FALSE)
    {
        if (!missing(all.names))
            warning("'all.names' is ignored when patterns have no names")
        x@ends
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "ByName_MIndex" class.
### 
### When 'pdict' is a PDict object with no names then matchPDict(pdict, ...)
### returns the matches in a ByName_MIndex object.
###
### Note that in normal operations the user NEVER needs to create a
### ByName_MIndex object explicitely or to modify an existing one:
### ByName_MIndex objects are created by the matchPDict() function
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
###       CWdna_PDict object. Another solution would be to keep the width slot
###       and to make it the same length as the ends_envir slot (it's currently
###       of length 1 only).
###
###   names: a character vector containing the _unique_ pattern names.
###

setClass("ByName_MIndex",
    contains="MIndex",
    representation(
        length="integer",
        ends_envir="environment",
        width="integer",
        NAMES="character" # R doesn't like @names !!
    )
)

setMethod("length", "ByName_MIndex", function(x) x@length)

setMethod("names", "ByName_MIndex", function(x) x@NAMES)

setMethod("show", "ByName_MIndex",
    function(object)
    {
        cat(length(object), "-pattern \"MIndex\" object\n", sep="")
    }
)

setMethod("[[", "ByName_MIndex",
    function(x, i, j, ...)
    {
        key <- callNextMethod()
        if (is.character(key)) {
            pos <- match(key, names(x))
            if (is.na(pos))
                stop("pattern name \"", key, "\" not found")
            key <- pos
        } 
        ans_end <- x@ends_envir[[formatC(key, width=10, format="d", flag="0")]]
        if (is.null(ans_end))
            ans_end <- integer(0)
        ans_width <- rep.int(x@width, length(ans_end))
        ans_start <- ans_end - ans_width + 1L
        new("IRanges", start=ans_start, width=ans_width, check=FALSE)
    }
)

### An example of a ByName_MIndex object of length 5 where only the
### 2nd pattern has matches:
###   > ends_envir <- new.env(hash=TRUE, parent=emptyenv())
###   > ends_envir[['0000000002']] <- c(199L, 402L)
###   > mindex <- new("ByName_MIndex", length=5L, ends_envir=ends_envir, width=10L, NAMES=letters[1:5])
###   > mindex[[1]]
###   > mindex[[2]]
###   > mindex[[6]] # Error in mindex[[6]] : subscript out of bounds
###   > names(mindex)
###   > mindex[["a"]]
###   > mindex[["b"]]
###   > mindex[["aa"]] # Error in mindex[["aa"]] : pattern name ‘aa’ not found
###   > startIndex(mindex)
###   > startIndex(mindex, all.names=TRUE)
###   > endIndex(mindex)
###   > endIndex(mindex, all.names=TRUE)
###   > countIndex(mindex)
###   > countIndex(mindex, all.names=TRUE)
###
setMethod("startIndex", "ByName_MIndex",
    function(x, all.names=FALSE)
    {
        if (!isTRUEorFALSE(all.names))
            stop("'all.names' must be 'TRUE' or 'FALSE'")
        .Call("extract_endIndex", x@ends_envir, 1L - x@width, x@NAMES, all.names, PACKAGE="Biostrings")
    }
)
setMethod("endIndex", "ByName_MIndex",
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

setMethod("unlist", "MIndex",
    function(x, recursive=TRUE, use.names=TRUE)
    {
        use.names <- normalize.use.names(use.names)
        start_index <- startIndex(x)
        if (length(start_index) == 0)
            ans_start <- integer(0)
        else
            ans_start <- unlist(start_index, recursive=FALSE, use.names=FALSE)
        end_index <- endIndex(x)
        if (length(end_index) == 0)
            ans_end <- integer(0)
        else
            ans_end <- unlist(end_index, recursive=FALSE, use.names=FALSE)
        if (use.names) {
            ans_names <- names(end_index)
            if (!is.null(ans_names))
                ans_names <- rep.int(ans_names, times=sapply(end_index, length))
        } else {
            ans_names <- NULL
        }
        ans_width <- ans_end - ans_start + 1L
        new("IRanges", start=ans_start, width=ans_width, names=ans_names, check=FALSE)
    }
)

extractAllMatches <- function(subject, mindex)
{
    if (!is(subject, "BString"))
        stop("'subject' must be a BString object")
    if (!is(mindex, "MIndex"))
        stop("'mindex' must be an MIndex object")
    if (is.null(names(mindex)))
        stop("extractAllMatches() works only with a \"MIndex\" object with names")
    allviews <- unlist(mindex)
    new("BStringViews", subject,
        start=start(allviews), width=width(allviews),
        desc=desc(allviews), check=FALSE)
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
###   > pdict <- PDict(dict)
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

.match.CWdna_PDict.exact <- function(pdict, subject, count.only)
{
    actree <- .ACtree.prepare_for_use_on_DNAString(pdict@actree)
    names <- names(pdict)
    if (is.null(names))
        envir <- NULL
    else
        envir <- new.env(hash=TRUE, parent=emptyenv())
    C_ans <- .Call("match_TBdna",
                   actree@nodes@xp, actree@base_codes,
                   pdict@dups, NULL, NULL,
                   subject,
                   0L, c(TRUE, TRUE),
                   count.only, envir,
                   PACKAGE="Biostrings")
    if (count.only)
        return(C_ans)
    if (is.null(names))
        new("ByPos_MIndex", ends=C_ans, width=width(pdict))
    else
        new("ByName_MIndex", length=length(pdict), ends_envir=C_ans, width=width(pdict), NAMES=names)
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
###      > pdict <- PDict(dict)
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
###   pdict <- PDict(c("acgt", "gt", "cgt", "ac"), tb.end=2)
###   endIndex(matchPDict(pdict, DNAString("acggaccg"), max.mismatch=0))
###   endIndex(matchPDict(pdict, DNAString("acggaccg"), max.mismatch=1))
###   endIndex(matchPDict(pdict, DNAString("acggaccg"), max.mismatch=2))
###
### At the moment, preprocessing a big TBdna_PDict object is very slow:
###   > library(drosophila2probe)
###   > dict0 <- drosophila2probe$sequence 
###   > system.time(pdict0 <- PDict(dict0[1:40000]))
###      user  system elapsed
###     0.040   0.032   0.072
###   > system.time(pdict <- PDict(dict0[1:40000], tb.end=10))
###      user  system elapsed
###    38.158   0.052  38.217
###
###   > library(BSgenome.Dmelanogaster.FlyBase.r51)
###   > chr3R <- Dmelanogaster[["3R"]]
###   > system.time(mindex0 <- matchPDict(pdict0, chr3R))
###      user  system elapsed
###     1.352   0.000   1.352
###   > system.time(mindex <- matchPDict(pdict, chr3R))
###      user  system elapsed
###     1.332   0.008   1.338
###   > identical(countIndex(mindex0), countIndex(mindex))
###   [1] TRUE
###
### Allowing mismatches is fast:
###   > system.time(mindex_mm6 <- matchPDict(pdict, chr3R, max.mismatch=4))
###      user  system elapsed
###     1.377   0.000   1.375
###   > mindex_mm6[[103]]
###        start      end width
###   1  9381276  9381285    10
###   2 16070100 16070109    10
###   > v <- views(chr3R, start(mindex_mm6[[103]]), end(mindex_mm6[[103]])+15)
###   > mismatch(dict0[103], v)
###   [[1]]
###   [1] 14 15 19 23 24 25
###
###   [[2]]
###   integer(0)

.match.TBdna_PDict <- function(pdict, subject, max.mismatch, fixed, count.only)
{
    actree <- .ACtree.prepare_for_use_on_DNAString(pdict@actree)
    names <- names(pdict)
    if (is.null(names))
        envir <- NULL
    else
        envir <- new.env(hash=TRUE, parent=emptyenv())
    C_ans <- .Call("match_TBdna",
                   actree@nodes@xp, actree@base_codes,
                   pdict@dups, pdict@head, pdict@tail,
                   subject,
                   max.mismatch, fixed,
                   count.only, envir,
                   PACKAGE="Biostrings")
    if (count.only)
        return(C_ans)
    if (is.null(names))
        new("ByPos_MIndex", ends=C_ans, width=width(pdict))
    else
        new("ByName_MIndex", length=length(pdict), ends_envir=C_ans, width=width(pdict), NAMES=names)
}




### =========================================================================
### E. THE "matchPDict" AND "countPDict" GENERIC FUNCTION AND METHODS
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .matchPDict()
###

.matchPDict <- function(pdict, subject, algorithm, max.mismatch, fixed, count.only=FALSE)
{
    if (!is(pdict, "CWdna_PDict"))
        stop("the pattern dictionary 'pdict' can only be a CWdna_PDict (or derived) object for now")
    if (!is(subject, "DNAString"))
        stop("'subject' can only be a DNAString object for now")
    if (!identical(algorithm, "auto"))
        stop("'algo' can only be '\"auto\"' for now")
    max.mismatch <- normalize.max.mismatch(max.mismatch)
    fixed <- normalize.fixed(fixed, class(subject))
    if (is(pdict, "TBdna_PDict"))
        return(.match.TBdna_PDict(pdict, subject, max.mismatch, fixed, count.only))
    if (max.mismatch != 0 || !all(fixed))
        stop("only TBdna_PDict dictionaries support inexact matching")
    .match.CWdna_PDict.exact(pdict, subject, count.only)
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

