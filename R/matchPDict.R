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

.XString.matchPDict <- function(pdict, subject,
                                max.mismatch, fixed,
                                count.only, envir)
{
    fixed <- normargFixed(fixed, subject)
    if (is.null(head(pdict)) && is.null(tail(pdict)))
        .checkUserArgsWhenTrustedBandIsFull(max.mismatch, fixed)
    .Call("XString_match_pdict",
          pdict@threeparts@pptb, head(pdict), tail(pdict),
          subject,
          max.mismatch, fixed,
          count.only, envir,
          PACKAGE="Biostrings")
}

.XStringViews.matchPDict <- function(pdict, subject,
                                     max.mismatch, fixed,
                                     count.only, envir)
{
    fixed <- normargFixed(fixed, subject)
    if (is.null(head(pdict)) && is.null(tail(pdict)))
        .checkUserArgsWhenTrustedBandIsFull(max.mismatch, fixed)
    .Call("XStringViews_match_pdict",
          pdict@threeparts@pptb, head(pdict), tail(pdict),
          subject(subject), start(subject), width(subject),
          max.mismatch, fixed,
          count.only, envir,
          PACKAGE="Biostrings")
}

.XStringSet.vmatchPDict <- function(pdict, subject,
                                    max.mismatch, fixed,
                                    count.only, envir)
{
    fixed <- normargFixed(fixed, subject)
    if (is.null(head(pdict)) && is.null(tail(pdict)))
        .checkUserArgsWhenTrustedBandIsFull(max.mismatch, fixed)
    .Call("XStringSet_vmatch_pdict",
          pdict@threeparts@pptb, head(pdict), tail(pdict),
          subject,
          max.mismatch, fixed,
          count.only, envir,
          PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .matchPDict()
###

.matchTB_PDict <- function(pdict, subject, algorithm,
                           max.mismatch, fixed, verbose, count.only)
{
    if (is(subject, "DNAString"))
        C_ans <- .XString.matchPDict(pdict, subject,
                                     max.mismatch, fixed,
                                     count.only, NULL)
    else if (is(subject, "XStringViews") && is(subject(subject), "DNAString"))
        C_ans <- .XStringViews.matchPDict(pdict, subject,
                                          max.mismatch, fixed,
                                          count.only, NULL)
    else
        stop("'subject' must be a DNAString object,\n",
             "  a MaskedDNAString object,\n",
             "  or an XStringViews object with a DNAString subject")
    if (count.only %in% c(TRUE, NA))  # whichPDict() or countPDict()
        return(C_ans)
    # matchPDict()
    new("ByPos_MIndex", ends=C_ans, width=width(pdict))
}

.matchMTB_PDict <- function(pdict, subject, algorithm,
                            max.mismatch, fixed, verbose, count.only)
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
    if (is.na(count.only)) # whichPDict()
        return(unique(sort(unlist(ans_parts))))
    if (verbose)
        cat("Combining the results obtained for ",
            "each TB_PDict component...\n", sep="")
    st <- system.time(ans <- ByPos_MIndex.combine(ans_parts), gcFirst=TRUE)
    if (verbose)
        print(st)
    if (count.only) # countPDict()
        return(countIndex(ans))
    # matchPDict()
    ans
}

.matchPDict <- function(pdict, subject, algorithm,
                        max.mismatch, fixed, verbose, count.only=FALSE)
{
    if (!identical(algorithm, "auto"))
        stop("'algorithm' can only be '\"auto\"' for now")
    max.mismatch <- normargMaxMismatch(max.mismatch)
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    if (is(pdict, "TB_PDict"))
        ans <- .matchTB_PDict(pdict, subject, algorithm,
                              max.mismatch, fixed, verbose, count.only)
    else if (is(pdict, "MTB_PDict"))
        ans <- .matchMTB_PDict(pdict, subject, algorithm,
                               max.mismatch, fixed, verbose, count.only)
    else
        stop("'pdict' must be a PDict object")
    if (is.null(dups(pdict)))
        return(ans)
    if (is.na(count.only)) # whichPDict()
        return(sort(c(ans, unlist(dups(pdict)@unq2dup[ans]))))
    if (count.only) { # countPDict()
        which_is_dup <- which(duplicated(pdict))
        ans[which_is_dup] <- ans[dups(pdict)@dup2unq[which_is_dup]]
        return(ans)
    }
    # matchPDict()
    if (is(ans, "ByPos_MIndex")) {
        ans@dups0 <- dups(pdict)
    } else {
        stop("don't know how to store the dup info in a ",
             class(ans), " object yet")
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .vmatchPDict()
###

.vmatchTB_PDict <- function(pdict, subject, algorithm,
                            max.mismatch, fixed, verbose, count.only)
{
    if (!count.only)
        stop("'count.only' must be TRUE (other values are not supported yet, sorry)")
    if (is(subject, "DNAStringSet"))
        C_ans <- .XStringSet.vmatchPDict(pdict, subject,
                                         max.mismatch, fixed,
                                         count.only, NULL)
    else
        stop("'subject' must be a DNAStringSet object")
    C_ans
}

.vmatchPDict <- function(pdict, subject, algorithm,
                         max.mismatch, fixed, verbose, count.only=FALSE)
{
    if (!count.only)
        stop("'count.only' must be TRUE (other values are not supported yet, sorry)")
    if (!identical(algorithm, "auto"))
        stop("'algorithm' can only be '\"auto\"' for now")
    max.mismatch <- normargMaxMismatch(max.mismatch)
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    if (is(pdict, "TB_PDict"))
        ans <- .vmatchTB_PDict(pdict, subject, algorithm,
                               max.mismatch, fixed, verbose, count.only)
    else if (is(pdict, "MTB_PDict"))
        stop("MTB_PDict objects are not supported yet, sorry")
    else
        stop("'pdict' must be a PDict object")
    if (is.null(dups(pdict)))
        return(ans)
    if (count.only) { # vcountPDict()
        which_is_dup <- which(duplicated(pdict))
        ans[which_is_dup, ] <- ans[dups(pdict)@dup2unq[which_is_dup], ]
        return(ans)
    }
    stop("internal problem")
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
                    max.mismatch, fixed, verbose, count.only=TRUE)
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "vmatchPDict", "vcountPDict" and "vwhichPDict" generic and methods.
###
### These are vectorized versions of matchPDict(), countPDict() and
### whichPDict().
### ONLY vcountPDict() is supported for now. It returns an M x N matrix of
### integers where M is the number of patterns (length(pdict)) and N the
### number of sequences in the subject (length(subject)).
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
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        standardGeneric("vcountPDict")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vcountPDict", "XString",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        stop("please use countPDict() when 'subject' is an XString object (single sequence)")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vcountPDict", "XStringSet",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        .vmatchPDict(pdict, subject, algorithm,
                     max.mismatch, fixed, verbose, count.only=TRUE)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vcountPDict", "XStringViews",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        vcountPDict(pdict, XStringViewsToSet(subject, FALSE, verbose=FALSE),
                    algorithm=algorithm,
                    max.mismatch=max.mismatch, fixed=fixed,
                    verbose=verbose)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vcountPDict", "MaskedXString",
    function(pdict, subject, algorithm="auto",
             max.mismatch=0, fixed=TRUE, verbose=FALSE)
        stop("please use countPDict() when 'subject' is a MaskedXString object (single sequence)")
)

