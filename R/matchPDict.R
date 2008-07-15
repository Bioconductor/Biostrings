#############################################################################
###
### The matchPDict() generic & related classes and functions
### ========================================================
###
### Author: Herve Pages (hpages@fhcrc.org)
###
#############################################################################




### =========================================================================
### A. SEARCH RESULT MANIPULATION
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
        # 'x' is guaranteed to be an MIndex object (if it's not, then
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
###       starts of the matches when the original dictionary has a variable
###       width. Another solution would be to keep the width slot and to make
###       it the same length as the ends slot (it's currently of length 1
###       only).
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
        cat(length(object), "-pattern MIndex object\n", sep="")
    }
)

setMethod("[[", "ByPos_MIndex",
    function(x, i, j, ...)
    {
        key <- callNextMethod()
        if (is.character(key))
            stop("MIndex object has no names")
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
            warning("'all.names' is ignored when MIndex object has no names")
        .Call("shiftListOfInts", x@ends, 1L - x@width, PACKAGE="Biostrings")
    }
)
setMethod("endIndex", "ByPos_MIndex",
    function(x, all.names=FALSE)
    {
        if (!missing(all.names))
            warning("'all.names' is ignored when MIndex object has no names")
        x@ends
    }
)

ByPos_MIndex.combine <- function(mi_list, ans_width)
{
    ends_listlist <- lapply(mi_list, function(mi) mi@ends)
    #mergeEnds <- function(...)
    #{
    #    ans <- unlist(list(...), recursive=FALSE, use.names=FALSE)
    #    if (is.null(ans))
    #        return(NULL)
    #    sort(unique(ans))
    #}
    #args <- c(list(FUN=mergeEnds), ends_listlist, list(SIMPLIFY=FALSE))
    #ans_ends <- do.call(mapply, args)
    ans_ends <- .Call("ByPos_MIndex_combine",
                      ends_listlist,
                      PACKAGE="Biostrings")
    new("ByPos_MIndex", ends=ans_ends, width=ans_width)
}

setMethod("combine", signature(x="ByPos_MIndex", y="ByPos_MIndex"),
    function(x, y, ...)
    {
        stop("not ready yet, sorry!")
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
###   length: the length of the original dictionary.
###
###   ends_envir: a key-value list (environment) where the values are integer
###       vectors containing the ending positions of the pattern whose
###       position in the original dictionary is given by the key (the keys are
###       strings representing positive integers).
###       
###   width: temporary hack. In the future we will probably want to store the
###       starts of the matches when the original dictionary has a variable
###       width. Another solution would be to keep the width slot and to make
###       it the same length as the ends_envir slot (it's currently of length
###       1 only).
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
        cat(length(object), "-pattern sparse MIndex object\n", sep="")
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
        use.names <- normargUseNames(use.names)
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
    if (is(subject, "MaskedXString"))
        subject <- unmasked(subject)
    else if (!is(subject, "XString"))
        stop("'subject' must be an XString or MaskedXString object")
    if (!is(mindex, "MIndex"))
        stop("'mindex' must be an MIndex object")
    if (is.null(names(mindex)))
        stop("extractAllMatches() works only with a \"MIndex\" object with names")
    allviews <- unlist(mindex)
    new("XStringViews", subject,
        start=start(allviews), width=width(allviews),
        names=names(allviews), check=FALSE)
}




### =========================================================================
### B. MATCH FINDING
### -------------------------------------------------------------------------
###
### Some examples below with a PDict object of type "ACtree".
### TODO: All these examples need to go somewhere else!
###
### I. A real use-case
### ------------------
###   > library(hgu95av2probe)
###   > dict <- DNAStringSet(hgu95av2probe$sequence) # the original dict
###   > pdict <- PDict(dict)
###   > pdict
###   > library(BSgenome.Hsapiens.UCSC.hg18)
###   > chr1 <- Hsapiens$chr1
###   > system.time(end_index <- endIndex(matchPDict(pdict, chr1)))
###      user  system elapsed 
###    50.663   0.000  50.763
###   > count_index <- sapply(end_index, length)
###   > table(count_index)
###   > id0 <- which(count_index == max(count_index))
###   > p0 <- pdict[[id0]]
###   > p0
###     25-letter "DNAString" instance
###   Value: CTGTAATCCCAGCACTTTGGGAGGC
###   > subXString(chr1, end_index[[id0]][1]-24, end_index[[id0]][1]) == p0
###   [1] TRUE
### For a more extensive validation:
###   > pidOK <- sapply(seq_len(length(end_index)),
###                     function(pid)
###                         identical(end_index[[pid]],
###                         end(matchPattern(pdict[[pid]], chr1))))
###   > all(pidOK)
### but be aware that THIS WILL TAKE THE WHOLE DAY!!! (20-24 hours)
###
### II. With a big random dictionary, on george1
### --------------------------------------------
### 1. Trying to simulate Solexa data:
###      > library(Biostrings)
###      > dict_length <- 10^6
###      > s <- DNAString(paste(sample(c("A", "C", "G", "T"), 36*dict_length,
###                                    replace=TRUE), collapse=""))
###      > views_start <- (0:(dict_length-1)) * 36 + 1
###      > dict <- views(s, views_start, views_start + 35) # the original dict
###
### 2. Building the Aho-Corasick 4-ary tree from the original dictionary:
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
### 
### III. Inexact matching
### ---------------------
###   pdict <- PDict(c("acgt", "gt", "cgt", "ac"), tb.end=2)
###   endIndex(matchPDict(pdict, DNAString("acggaccg"), max.mismatch=0))
###   endIndex(matchPDict(pdict, DNAString("acggaccg"), max.mismatch=1))
###   endIndex(matchPDict(pdict, DNAString("acggaccg"), max.mismatch=2))
###
### TODO: Rerun the benchmarks below on the entire dict0.
###   > library(drosophila2probe)
###   > dict0 <- DNAStringSet(drosophila2probe$sequence)
###   > system.time(pdict0 <- PDict(dict0[1:40000]))
###      user  system elapsed
###     0.040   0.032   0.072
###   > system.time(pdict <- PDict(dict0[1:40000], tb.end=10))
###      user  system elapsed
###    38.158   0.052  38.217
###
###   > library(BSgenome.Dmelanogaster.UCSC.dm3)
###   > chr3R <- Dmelanogaster$chr3R
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
###   > mismatch(dict0[[103]], v)
###   [[1]]
###   [1] 14 15 19 23 24 25
###
###   [[2]]
###   integer(0)

.matchTB_PDict <- function(pdict, subject, algorithm,
                           max.mismatch, fixed, verbose, count.only)
{
    pptb <- pdict@threeparts@pptb
    pdict_type <- class(pptb)
    if (pdict_type == "Twobit") {
        pdict_pptb <- list(tb.width(pptb),
                           pptb@sign2pos@xp,
                           pptb@base_codes)
    } else if (pdict_type == "ACtree") {
        pdict_pptb <- list(pptb@nodes@xp,
                           pptb@base_codes)
    } else {
        stop(pdict_type, ": unsupported type of PDict object for 'pdict'")
    }

    pdict_head <- head(pdict)
    pdict_tail <- tail(pdict)
    if (is.null(pdict_head) && is.null(pdict_tail) && max.mismatch != 0)
        stop("'max.mismatch' must be 0 when there is no head and no tail")
    #THIS IS NOT READY YET!
    #clean_tb_unq2dup <- clean.tb.unq2dup(pdict)
    #unq2dup <- if (!is.null(clean_tb_unq2dup)) clean_tb_unq2dup else dups(pptb)@unq2dup
    unq2dup <- dups(pptb)@unq2dup
    pdict_names <- names(pdict)
    if (is.null(pdict_names))
        envir <- NULL
    else
        envir <- new.env(hash=TRUE, parent=emptyenv())
    if (is(subject, "DNAString")) {
        fixed <- normargFixed(fixed, class(subject))
        if (is.null(pdict_head) && is.null(pdict_tail) && !fixed[1])
            warning("with this value of 'fixed', IUPAC extended letters in\n",
                    "the patterns will be treated as ambiguities, even\n",
                    "though the patterns don't contain such letters")
        C_ans <- .Call("XString_match_pdict",
                       pdict_type, pdict_pptb,
                       length(pdict), tb.width(pptb),
                       unq2dup,
                       pdict_head, pdict_tail,
                       subject,
                       max.mismatch, fixed,
                       count.only, envir,
                       PACKAGE="Biostrings")
    } else if (is(subject, "XStringViews")
            && is(subject(subject), "DNAString")) {
        fixed <- normargFixed(fixed, class(subject(subject)))
        if (is.null(pdict_head) && is.null(pdict_tail) && !fixed[1])
            warning("with this value of 'fixed', IUPAC extended letters in\n",
                    "the patterns will be treated as ambiguities, even\n",
                    "though the patterns don't contain such letters")
        C_ans <- .Call("XStringViews_match_pdict",
                       pdict_type, pdict_pptb,
                       length(pdict), tb.width(pptb),
                       unq2dup,
                       pdict_head, pdict_tail,
                       subject(subject), start(subject), width(subject),
                       max.mismatch, fixed,
                       count.only, envir,
                       PACKAGE="Biostrings")
    } else {
        stop("'subject' must be a DNAString object,\n",
             "  a MaskedDNAString object,\n",
             "  or an XStringViews object with a DNAString subject")
    }
    if (count.only %in% c(TRUE, NA))
        return(C_ans)
    if (is.null(pdict_names))
        new("ByPos_MIndex", ends=C_ans, width=tb.width(pdict))
    else
        new("ByName_MIndex", length=length(pdict), ends_envir=C_ans,
                             width=tb.width(pdict), NAMES=pdict_names)
}

.matchMTB_PDict <- function(pdict, subject, algorithm,
                            max.mismatch, fixed, verbose, count.only)
{
    pdict_names <- names(pdict)
    if (!is.null(pdict_names))
        stop("MTB_PDict objects with names are not supported yet")
    tb_pdicts <- as.list(pdict)
    NTB <- length(tb_pdicts)
    max.mismatch0 <- NTB - 1L
    if (max.mismatch > max.mismatch0)
        stop("cannot use matchPDict()/countPDict()/whichPDict() ",
             "with an MTB_PDict object that was preprocessed to be used ",
             "with 'max.mismatch' <= ", max.mismatch0)
    if (max.mismatch < max.mismatch0)
        cat("WARNING: using 'max.mismatch=", max.mismatch, "' with an ",
            "MTB_PDict object that was preprocessed for allowing up to ",
            max.mismatch0, " mismatching letter(s) per match is not optimal\n",
            sep="")
    ans_parts <- lapply(seq_len(NTB),
        function(i)
        {
            if (verbose)
                cat("Getting results for TB_PDict component ",
                    i, "/", NTB, " ...\n", sep="")
            tb_pdict <- tb_pdicts[[i]]
            st <- system.time(
                {
                    ans_part <- .matchTB_PDict(tb_pdict, subject, algorithm,
                                               max.mismatch, fixed, verbose,
                                               (if (is.na(count.only)) NA else FALSE))
                }, gcFirst=TRUE)
            if (verbose) {
                print(st)
                cat(sum(countIndex(ans_part)), " match(es) found\n", sep="")
            }
            ans_part
        }
    )
    if (is.na(count.only))
        return(unique(sort(unlist(ans_parts))))
    if (is.null(pdict_names)) {
        if (verbose)
            cat("Combining the results obtained for ",
                "each TB_PDict component...\n", sep="")
        st <- system.time(
            {
                ans_width <- sum(sapply(tb_pdicts, tb.width))
                ans <- ByPos_MIndex.combine(ans_parts, ans_width)
            }, gcFirst=TRUE)
        if (verbose)
            print(st)
    } else {
        stop("not ready yet")
    }
    if (count.only)
        return(countIndex(ans))
    ans
}

.matchPDict <- function(pdict, subject, algorithm,
                        max.mismatch, fixed, verbose, count.only=FALSE)
{
    if (!identical(algorithm, "auto"))
        stop("'algorithm' can only be '\"auto\"' for now")
    max.mismatch <- normargMaxMismatch(max.mismatch)
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be 'TRUE' or 'FALSE'")
    if (is(pdict, "TB_PDict"))
        return(.matchTB_PDict(pdict, subject, algorithm,
                              max.mismatch, fixed, verbose, count.only))
    if (is(pdict, "MTB_PDict"))
        return(.matchMTB_PDict(pdict, subject, algorithm,
                               max.mismatch, fixed, verbose, count.only))
    stop("'pdict' must be a PDict object")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "matchPDict" generic and methods.
###

setGeneric("matchPDict", signature="subject",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        standardGeneric("matchPDict")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPDict", "XString",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        .matchPDict(pdict, subject, algorithm, max.mismatch, fixed, verbose)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPDict", "XStringViews",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        .matchPDict(pdict, subject, algorithm, max.mismatch, fixed, verbose)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPDict", "MaskedXString",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        matchPDict(pdict, toXStringViewsOrXString(subject), algorithm=algorithm,
                   max.mismatch=max.mismatch, fixed=fixed, verbose=verbose)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "countPDict" generic and methods.
###

setGeneric("countPDict", signature="subject",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        standardGeneric("countPDict")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPDict", "XString",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        .matchPDict(pdict, subject, algorithm,
                    max.mismatch, fixed, verbose, count.only=TRUE)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPDict", "XStringViews",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        .matchPDict(pdict, subject, algorithm,
                    max.mismatch, fixed, verbose, count.only=TRUE)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPDict", "MaskedXString",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        countPDict(pdict, toXStringViewsOrXString(subject), algorithm=algorithm,
                   max.mismatch=max.mismatch, fixed=fixed, verbose=verbose)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "whichPDict" generic and methods.
###

setGeneric("whichPDict", signature="subject",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        standardGeneric("whichPDict")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("whichPDict", "XString",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        .matchPDict(pdict, subject, algorithm,
                    max.mismatch, fixed, verbose, count.only=NA)
)

