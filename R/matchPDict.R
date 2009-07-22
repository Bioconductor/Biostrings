### =========================================================================
### The matchPDict() generic & related functions
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
###   > subseq(chr1, start=end_index[[id0]][1]-24, end=end_index[[id0]][1]) == p0
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
###      > s <- DNAString(paste(sample(DNA_BASES, 36*dict_length,
###                                    replace=TRUE), collapse=""))
###      > views_start <- (0:(dict_length-1)) * 36 + 1
###      > dict <- Views(s, start=views_start, end=views_start + 35) # the original dict
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
###   > v <- Views(chr3R, start=start(mindex_mm6[[103]]), end=end(mindex_mm6[[103]])+15)
###   > mismatch(dict0[[103]], v)
###   [[1]]
###   [1] 14 15 19 23 24 25
###
###   [[2]]
###   integer(0)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Convenience wrappers to .Call().
###

.checkUserArgsWhenTrustedBandIsFull <- function(max.mismatch, fixed)
{
    if (max.mismatch != 0)
        stop("'max.mismatch' must be 0 when there is no head and no tail")
    if (!fixed[1])
        warning("the value you specified for 'fixed' means that IUPAC extended\n",
                "letters in the patterns should be treated as ambiguities, but\n",
                "this has no effect because the patterns don't contain such\n",
                "letters")
}

### 'threeparts' is a PDict3Parts object.
.match.PDict3Parts.XString <- function(threeparts, subject,
                                       max.mismatch, fixed,
                                       matches.as, envir)
{
    fixed <- normargFixed(fixed, subject)
    if (is.null(head(threeparts)) && is.null(tail(threeparts)))
        .checkUserArgsWhenTrustedBandIsFull(max.mismatch, fixed)
    .Call("XString_match_pdict",
          threeparts@pptb, head(threeparts), tail(threeparts),
          subject,
          max.mismatch, fixed,
          matches.as, envir,
          PACKAGE="Biostrings")
}

### 'threeparts' is a PDict3Parts object.
.match.PDict3Parts.XStringViews <- function(threeparts, subject,
                                            max.mismatch, fixed,
                                            matches.as, envir)
{
    fixed <- normargFixed(fixed, subject)
    if (is.null(head(threeparts)) && is.null(tail(threeparts)))
        .checkUserArgsWhenTrustedBandIsFull(max.mismatch, fixed)
    .Call("XStringViews_match_pdict",
          threeparts@pptb, head(threeparts), tail(threeparts),
          subject(subject), start(subject), width(subject),
          max.mismatch, fixed,
          matches.as, envir,
          PACKAGE="Biostrings")
}

### 'threeparts' is a PDict3Parts object.
.vmatch.PDict3Parts.XStringSet <- function(threeparts, subject,
                                           max.mismatch, fixed,
                                           collapse, weight,
                                           matches.as, envir)
{
    fixed <- normargFixed(fixed, subject)
    if (is.null(head(threeparts)) && is.null(tail(threeparts)))
        .checkUserArgsWhenTrustedBandIsFull(max.mismatch, fixed)
    .Call("XStringSet_vmatch_pdict",
          threeparts@pptb, head(threeparts), tail(threeparts),
          subject,
          max.mismatch, fixed,
          collapse, weight,
          matches.as, envir,
          PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .matchPDict()
###

.match.TB_PDict <- function(pdict, subject, algorithm,
                            max.mismatch, fixed, verbose, matches.as)
{
    if (is(subject, "DNAString"))
        C_ans <- .match.PDict3Parts.XString(pdict@threeparts, subject,
                                            max.mismatch, fixed,
                                            matches.as, NULL)
    else if (is(subject, "XStringViews") && is(subject(subject), "DNAString"))
        C_ans <- .match.PDict3Parts.XStringViews(pdict@threeparts, subject,
                                                 max.mismatch, fixed,
                                                 matches.as, NULL)
    else
        stop("'subject' must be a DNAString object,\n",
             "  a MaskedDNAString object,\n",
             "  or an XStringViews object with a DNAString subject")
    if (matches.as != "MATCHES_AS_ENDS")
        return(C_ans)
    # matchPDict()
    new("ByPos_MIndex", width0=width(pdict), NAMES=names(pdict), ends=C_ans)
}

.match.MTB_PDict <- function(pdict, subject, algorithm,
                             max.mismatch, fixed, verbose, matches.as)
{
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
                    if (matches.as == "MATCHES_AS_COUNTS")
                        matches.as2 <- "MATCHES_AS_ENDS"
                    else
                        matches.as2 <- matches.as
                    ans_part <- .match.TB_PDict(tb_pdict, subject, algorithm,
                                                max.mismatch, fixed, verbose,
                                                matches.as2)
                }, gcFirst=TRUE)
            if (verbose) {
                print(st)
                cat(sum(countIndex(ans_part)), " match(es) found\n", sep="")
            }
            ans_part
        }
    )
    if (matches.as == "MATCHES_AS_WHICH")
        return(unique(sort(unlist(ans_parts))))
    if (verbose)
        cat("Combining the results obtained for ",
            "each TB_PDict component...\n", sep="")
    st <- system.time(ans <- ByPos_MIndex.combine(ans_parts), gcFirst=TRUE)
    if (verbose)
        print(st)
    if (matches.as == "MATCHES_AS_COUNTS")
        return(countIndex(ans))
    return(ans)
}

.matchPDict <- function(pdict, subject, algorithm,
                        max.mismatch, fixed,
                        verbose, matches.as="MATCHES_AS_ENDS")
{
    which_pp_excluded <- NULL
    dups0 <- dups(pdict)
    if (!is.null(dups0))
        which_pp_excluded <- which(duplicated(dups0))
    if (!identical(algorithm, "auto"))
        stop("'algorithm' can only be '\"auto\"' for now")
    max.mismatch <- normargMaxMismatch(max.mismatch)
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    if (is(pdict, "TB_PDict"))
        ans <- .match.TB_PDict(pdict, subject, algorithm,
                               max.mismatch, fixed, verbose, matches.as)
    else if (is(pdict, "MTB_PDict"))
        ans <- .match.MTB_PDict(pdict, subject, algorithm,
                                max.mismatch, fixed, verbose, matches.as)
    else
        stop("'pdict' must be a PDict object")
    if (length(which_pp_excluded) == 0L)
        return(ans)
    if (matches.as == "MATCHES_AS_WHICH")
        return(members(dups0, ans))
    if (matches.as == "MATCHES_AS_COUNTS") {
        ans[which_pp_excluded] <- ans[togroup(dups0, which_pp_excluded)]
        return(ans)
    }
    if (is(ans, "ByPos_MIndex")) {
        ans@dups0 <- dups0
    } else {
        stop("don't know how to store the dup info in a ",
             class(ans), " object")
    }
    return(ans)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .vmatchPDict()
###

.vmatchPDict <- function(pdict, subject, algorithm,
                         max.mismatch, fixed,
                         collapse, weight,
                         verbose, matches.as="MATCHES_AS_ENDS")
{
    which_pp_excluded <- NULL
    dups0 <- dups(pdict)
    if (!is.null(dups0))
        which_pp_excluded <- which(duplicated(dups0))
    if (!is(subject, "DNAStringSet"))
        stop("'subject' must be a DNAStringSet object")
    if (!identical(algorithm, "auto"))
        stop("'algorithm' can only be '\"auto\"' for now")
    max.mismatch <- normargMaxMismatch(max.mismatch)
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    if (matches.as == "MATCHES_AS_WHICH") {
        ## vwhichPDict()
    } else if (matches.as == "MATCHES_AS_COUNTS") {
        ## vcountPDict()
        collapse <- normargCollapse(collapse)
        if (collapse) {
            if (collapse == 1L) {
                weight <- normargWeight(weight, length(subject))
            } else {
                weight <- normargWeight(weight, length(pdict))
                if (length(which_pp_excluded) != 0L) {
                    ## Collapse weights of duplicates.
                    ## TODO: Implement this in C.
                    which_is_not_unique <- which(!sapply(low2high(dups0), is.null))
                    weight[which_is_not_unique] <-
                        weight[which_is_not_unique] +
                        sapply(which_is_not_unique,
                               function(i) sum(weight[low2high(dups0)[[i]]]))
                }
            }
        } else {
            if (!identical(weight, 1L))
                warning("'weight' is ignored when 'collapse=FALSE'")
        }
    } else {
        ## vmatchPDict()
        stop("vmatchPDict() is not supported yet, sorry")
    }
    if (is(pdict, "TB_PDict"))
        ans <- .vmatch.PDict3Parts.XStringSet(pdict@threeparts, subject,
                                              max.mismatch, fixed,
                                              collapse, weight, matches.as, NULL)
    else if (is(pdict, "MTB_PDict"))
        stop("MTB_PDict objects are not supported yet, sorry")
    else
        stop("'pdict' must be a PDict object")
    if (length(which_pp_excluded) == 0L)
        return(ans)
    if (matches.as == "MATCHES_AS_WHICH") {
        ## vwhichPDict()
        return(vmembers(dups0, ans))
    }
    if (matches.as == "MATCHES_AS_COUNTS") {
        ## vcountPDict()
        if (collapse == 0L) {
            ans[which_pp_excluded, ] <- ans[togroup(dups0, which_pp_excluded), ]
        } else if (collapse == 1L) {
            ans[which_pp_excluded] <- ans[togroup(dups0, which_pp_excluded)]
        }
        return(ans)
    }
    ## vmatchPDict()
    stop("vmatchPDict() is not supported yet, sorry")
    return(ans)
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
setMethod("matchPDict", "XStringSet",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        stop("please use vmatchPDict() when 'subject' is an XStringSet object (multiple sequence)")
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
                    max.mismatch, fixed, verbose, matches.as="MATCHES_AS_COUNTS")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPDict", "XStringSet",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        stop("please use vcountPDict() when 'subject' is an XStringSet object (multiple sequence)")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPDict", "XStringViews",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        .matchPDict(pdict, subject, algorithm,
                    max.mismatch, fixed, verbose, matches.as="MATCHES_AS_COUNTS")
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
                    max.mismatch, fixed, verbose, matches.as="MATCHES_AS_WHICH")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "vmatchPDict", "vcountPDict" and "vwhichPDict" generic and methods.
###
### These are vectorized versions of matchPDict(), countPDict() and
### whichPDict(). If M denotes the number of patterns in 'pdict' and N the
### number of sequences in 'subject', then:
###   o vwhichPDict(): 'vwhichPDict(pdict, subject, ...)' is equivalent to
###     'lapply(subject, function(s) whichPDict(pdict, s, ...))'.
###      The returned object is a list of length N.
###   o vcountPDict(): returns an M x N matrix of integers.
###   o vmatchPDict(): not supported yet!
###

setGeneric("vmatchPDict", signature="subject",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        standardGeneric("vmatchPDict")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vmatchPDict", "ANY",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        stop("vmatchPDict() is not ready yet, sorry")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vmatchPDict", "XString",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        stop("please use matchPDict() when 'subject' is an XString object (single sequence)")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vmatchPDict", "MaskedXString",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        stop("please use matchPDict() when 'subject' is a MaskedXString object (single sequence)")
)

setGeneric("vcountPDict", signature="subject",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE,
             collapse=FALSE, weight=1L, verbose=FALSE)
        standardGeneric("vcountPDict")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vcountPDict", "XString",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE,
             collapse=FALSE, weight=1L, verbose=FALSE)
        stop("please use countPDict() when 'subject' is an XString object (single sequence)")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vcountPDict", "XStringSet",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE,
             collapse=FALSE, weight=1L, verbose=FALSE)
        .vmatchPDict(pdict, subject, algorithm,
                     max.mismatch, fixed,
                     collapse, weight, verbose, matches.as="MATCHES_AS_COUNTS")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vcountPDict", "XStringViews",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE,
             collapse=FALSE, weight=1L, verbose=FALSE)
        vcountPDict(pdict, XStringViewsToSet(subject, FALSE, verbose=FALSE),
                    algorithm=algorithm,
                    max.mismatch=max.mismatch, fixed=fixed,
                    collapse=collapse, weight=weight, verbose=verbose)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vcountPDict", "MaskedXString",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE,
             collapse=FALSE, weight=1L, verbose=FALSE)
        stop("please use countPDict() when 'subject' is a MaskedXString object (single sequence)")
)

setGeneric("vwhichPDict", signature="subject",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        standardGeneric("vwhichPDict")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vwhichPDict", "XString",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        stop("please use whichPDict() when 'subject' is an XString object (single sequence)")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vwhichPDict", "XStringSet",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        .vmatchPDict(pdict, subject, algorithm,
                     max.mismatch, fixed,
                     0L, 1L, verbose, matches.as="MATCHES_AS_WHICH")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vwhichPDict", "XStringViews",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        vwhichPDict(pdict, XStringViewsToSet(subject, FALSE, verbose=FALSE),
                    algorithm=algorithm,
                    max.mismatch=max.mismatch, fixed=fixed, verbose=verbose)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vwhichPDict", "MaskedXString",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        stop("please use whichPDict() when 'subject' is a MaskedXString object (single sequence)")
)

