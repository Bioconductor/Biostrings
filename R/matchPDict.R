### =========================================================================
### The matchPDict() generic & related functions
### -------------------------------------------------------------------------
###
### Some examples below with a PDict object of type "ACtree".
### TODO: Update the timings obtained when using the ACtree2 algo. All these
### examples need to go somewhere else!
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
###   > identical(elementNROWS(mindex0), elementNROWS(mindex))
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
### Checking the user-supplied arguments.
###

.checkUserArgsWhenTrustedBandIsFull <- function(max.mismatch, fixed)
{
    if (max.mismatch != 0L)
        stop("'max.mismatch' must be 0 when there is no head and no tail")
    if (!fixed[1L])
        warning("the value you specified for 'fixed' means that IUPAC\n",
                "extended letters in the patterns should be treated as\n",
                "ambiguities, but this has no effect because the patterns\n",
                "don't contain such letters")
}

.checkMaxMismatch <- function(max.mismatch, NTB)
{
    max.mismatch0 <- NTB - 1L
    if (max.mismatch > max.mismatch0)
        stop("cannot use vmatchPDict()/vcountPDict()/vwhichPDict() ",
             "with an MTB_PDict object that was preprocessed to be used ",
             "with 'max.mismatch' <= ", max.mismatch0)
    if (max.mismatch < max.mismatch0)
        cat("WARNING: using 'max.mismatch=", max.mismatch, "' with an ",
            "MTB_PDict object that was preprocessed for allowing up to ",
            max.mismatch0, " mismatching letter(s) per match is not optimal\n",
            sep="")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Convenience wrappers to .Call2().
###

### 'threeparts' is a PDict3Parts object.
.match.PDict3Parts.XString <- function(threeparts, subject,
                max.mismatch, min.mismatch, with.indels, fixed,
                algorithm, matches.as, envir)
{
    fixed <- normargFixed(fixed, subject)
    with.indels <- normargWithIndels(with.indels)
    if (with.indels)
        stop("at the moment, matchPDict() and family only support indels ",
             "on a non-preprocessed pattern dictionary, sorry")
    if (!identical(algorithm, "auto"))
        warning("'algorithm' is ignored when 'pdict' is a PDict object")
    if (is.null(head(threeparts)) && is.null(tail(threeparts)))
        .checkUserArgsWhenTrustedBandIsFull(max.mismatch, fixed)
    .Call2("match_PDict3Parts_XString",
          threeparts@pptb, head(threeparts), tail(threeparts),
          subject,
          max.mismatch, min.mismatch, fixed,
          matches.as, envir,
          PACKAGE="Biostrings")
}

### 'pattern' is an XStringSet object.
.match.XStringSet.XString <- function(pattern, subject,
                max.mismatch, min.mismatch, with.indels, fixed,
                algorithm, matches.as, envir)
{
    fixed <- normargFixed(fixed, subject)
    with.indels <- normargWithIndels(with.indels)
    if (with.indels &&
        !(matches.as %in% c("MATCHES_AS_WHICH", "MATCHES_AS_COUNTS")))
        stop("at the moment, within the matchPDict family, only ",
             "countPDict(), whichPDict(), vcountPDict() and vwhichPDict() ",
             "support indels")
    algo <- normargAlgorithm(algorithm)
    algo <- selectAlgo(algo, pattern, max.mismatch, min.mismatch,
                       with.indels, fixed)
    .Call2("match_XStringSet_XString",
          pattern,
          subject,
          max.mismatch, min.mismatch, with.indels, fixed,
          algo, matches.as, envir,
          PACKAGE="Biostrings")
}

### 'threeparts' is a PDict3Parts object.
.match.PDict3Parts.XStringViews <- function(threeparts, subject,
                max.mismatch, min.mismatch, with.indels, fixed,
                algorithm, matches.as, envir)
{
    fixed <- normargFixed(fixed, subject)
    with.indels <- normargWithIndels(with.indels)
    if (with.indels)
        stop("at the moment, matchPDict() and family only support indels ",
             "on a non-preprocessed pattern dictionary, sorry")
    if (!identical(algorithm, "auto"))
        warning("'algorithm' is ignored when 'pdict' is a PDict object")
    if (is.null(head(threeparts)) && is.null(tail(threeparts)))
        .checkUserArgsWhenTrustedBandIsFull(max.mismatch, fixed)
    .Call2("match_PDict3Parts_XStringViews",
          threeparts@pptb, head(threeparts), tail(threeparts),
          subject(subject), start(subject), width(subject),
          max.mismatch, min.mismatch, fixed,
          matches.as, envir,
          PACKAGE="Biostrings")
}

### 'pattern' is an XStringSet object.
.match.XStringSet.XStringViews <- function(pattern, subject,
                max.mismatch, min.mismatch, with.indels, fixed,
                algorithm, matches.as, envir)
{
    fixed <- normargFixed(fixed, subject)
    with.indels <- normargWithIndels(with.indels)
    if (with.indels &&
        !(matches.as %in% c("MATCHES_AS_WHICH", "MATCHES_AS_COUNTS")))
        stop("at the moment, within the matchPDict family, only ",
             "countPDict(), whichPDict(), vcountPDict() and vwhichPDict() ",
             "support indels")
    algo <- normargAlgorithm(algorithm)
    algo <- selectAlgo(algo, pattern, max.mismatch, min.mismatch,
                       with.indels, fixed)
    .Call2("match_XStringSet_XStringViews",
          pattern,
          subject(subject), start(subject), width(subject),
          max.mismatch, min.mismatch, with.indels, fixed,
          algo, matches.as, envir,
          PACKAGE="Biostrings")
}

### 'threeparts' is a PDict3Parts object.
.vmatch.PDict3Parts.XStringSet <- function(threeparts, subject,
                max.mismatch, min.mismatch, with.indels, fixed,
                algorithm, collapse, weight,
                matches.as, envir)
{
    fixed <- normargFixed(fixed, subject)
    with.indels <- normargWithIndels(with.indels)
    if (with.indels)
        stop("at the moment, matchPDict() and family only support indels ",
             "on a non-preprocessed pattern dictionary, sorry")
    if (!identical(algorithm, "auto"))
        warning("'algorithm' is ignored when 'pdict' is a PDict object")
    if (is.null(head(threeparts)) && is.null(tail(threeparts)))
        .checkUserArgsWhenTrustedBandIsFull(max.mismatch, fixed)
    .Call2("vmatch_PDict3Parts_XStringSet",
          threeparts@pptb, head(threeparts), tail(threeparts),
          subject,
          max.mismatch, min.mismatch, fixed,
          collapse, weight,
          matches.as, envir,
          PACKAGE="Biostrings")
}

### 'pattern' is an XStringSet object.
.vmatch.XStringSet.XStringSet <- function(pattern, subject,
                max.mismatch, min.mismatch, with.indels, fixed,
                algorithm, collapse, weight,
                matches.as, envir)
{
    fixed <- normargFixed(fixed, subject)
    with.indels <- normargWithIndels(with.indels)
    if (with.indels &&
        !(matches.as %in% c("MATCHES_AS_WHICH", "MATCHES_AS_COUNTS")))
        stop("at the moment, within the matchPDict family, only ",
             "countPDict(), whichPDict(), vcountPDict() and vwhichPDict() ",
             "support indels")
    algo <- normargAlgorithm(algorithm)
    algo <- selectAlgo(algo, pattern, max.mismatch, min.mismatch,
                       with.indels, fixed)
    .Call2("vmatch_XStringSet_XStringSet",
          pattern,
          subject,
          max.mismatch, min.mismatch, with.indels, fixed,
          algo, collapse, weight,
          matches.as, envir,
          PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining the results obtained for each component of an MTB_PDict object.
###

.combine.which.compons <- function(ans_compons)
    unique(sort(unlist(ans_compons)))

.combine.vwhich.compons <- function(ans_compons)
{
    lapply(seq_len(length(ans_compons[[1L]])),
        function(i)
        {
            ans_compons_elts <- lapply(ans_compons, "[[", i)
            .combine.which.compons(ans_compons_elts)
        }
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .matchPDict()
###

### 'pdict' is a TB_PDict object.
.match.TB_PDict <- function(pdict, subject,
                            max.mismatch, min.mismatch, with.indels, fixed,
                            algorithm, verbose, matches.as)
{
    if (is(subject, "DNAString"))
        C_ans <- .match.PDict3Parts.XString(pdict@threeparts, subject,
                     max.mismatch, min.mismatch, with.indels, fixed,
                     algorithm, matches.as, NULL)
    else if (is(subject, "XStringViews") && is(subject(subject), "DNAString"))
        C_ans <- .match.PDict3Parts.XStringViews(pdict@threeparts, subject,
                     max.mismatch, min.mismatch, with.indels, fixed,
                     algorithm, matches.as, NULL)
    else
        stop("'subject' must be a DNAString object,\n",
             "  a MaskedDNAString object,\n",
             "  or an XStringViews object with a DNAString subject")
    if (matches.as != "MATCHES_AS_ENDS")
        return(C_ans)
    # matchPDict()
    new("ByPos_MIndex", width0=width(pdict), NAMES=names(pdict), ends=C_ans)
}

### 'pdict' is an MTB_PDict object.
.match.MTB_PDict <- function(pdict, subject,
                             max.mismatch, min.mismatch, with.indels, fixed,
                             algorithm, verbose, matches.as)
{
    tb_pdicts <- as.list(pdict)
    NTB <- length(tb_pdicts)
    .checkMaxMismatch(max.mismatch, NTB)
    if (matches.as == "MATCHES_AS_COUNTS")
        matches.as2 <- "MATCHES_AS_ENDS"
    else
        matches.as2 <- matches.as
    ans_compons <- lapply(seq_len(NTB),
        function(i)
        {
            if (verbose)
                cat("Getting results for TB_PDict component ",
                    i, "/", NTB, " ...\n", sep="")
            tb_pdict <- tb_pdicts[[i]]
            st <- system.time({
                ans_compon <- .match.TB_PDict(tb_pdict, subject,
                                max.mismatch, min.mismatch, with.indels, fixed,
                                algorithm, verbose, matches.as2)
                  }, gcFirst=TRUE)
            if (verbose) {
                print(st)
                cat(sum(elementNROWS(ans_compon)), " match(es) found\n",
                    sep="")
            }
            ans_compon
        }
    )
    if (verbose)
        cat("Combining the results obtained for ",
            "each TB_PDict component...\n", sep="")
    if (matches.as == "MATCHES_AS_WHICH")
        return(.combine.which.compons(ans_compons))
    st <- system.time(ans <- ByPos_MIndex.combine(ans_compons), gcFirst=TRUE)
    if (verbose)
        print(st)
    if (matches.as == "MATCHES_AS_COUNTS")
        return(elementNROWS(ans))
    return(ans)
}

### 'pattern' is an XStringSet object.
.match.XStringSet <- function(pattern, subject,
                              max.mismatch, min.mismatch, with.indels, fixed,
                              algorithm, verbose, matches.as)
{
    if (is(subject, "XString"))
        C_ans <- .match.XStringSet.XString(pattern, subject,
                     max.mismatch, min.mismatch, with.indels, fixed,
                     algorithm, matches.as, NULL)
    else if (is(subject, "XStringViews"))
        C_ans <- .match.XStringSet.XStringViews(pattern, subject,
                     max.mismatch, min.mismatch, with.indels, fixed,
                     algorithm, matches.as, NULL)
    else
        stop("unsupported 'subject' type")
    if (matches.as != "MATCHES_AS_ENDS")
        return(C_ans)
    # matchPDict()
    new("ByPos_MIndex", width0=width(pattern), NAMES=names(pattern), ends=C_ans)
}

.matchPDict <- function(pdict, subject,
                        max.mismatch, min.mismatch, with.indels, fixed,
                        algorithm, verbose, matches.as="MATCHES_AS_ENDS")
{
    which_pp_excluded <- NULL
    if (is(pdict, "PDict")) {
        if (seqtype(subject) != "DNA")
            stop("'subject' must be DNA")
        dups0 <- dups(pdict)
        if (!is.null(dups0))
            which_pp_excluded <- which(duplicated(dups0))
    } else if (is(pdict, "XStringSet")) {
        if (seqtype(pdict) != seqtype(subject))
            stop("'pdict' and 'subject' must contain ",
                 "sequences of the same type")
    } else {
        stop("'pdict' must be a PDict or XStringSet object")
    }
    max.mismatch <- normargMaxMismatch(max.mismatch)
    min.mismatch <- normargMinMismatch(min.mismatch, max.mismatch)
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    ## We are doing our own dispatch here, based on the type of 'pdict'.
    ## TODO: Revisit this. Would probably be a better design to use a
    ## generic/methods approach and rely on the standard dispatch mechanism.
    ## Those low-level generic/methods wouldn't need to be exported.
    if (is(pdict, "TB_PDict"))
        ans <- .match.TB_PDict(pdict, subject,
                       max.mismatch, min.mismatch, with.indels, fixed,
                       algorithm, verbose, matches.as)
    else if (is(pdict, "MTB_PDict"))
        ans <- .match.MTB_PDict(pdict, subject,
                       max.mismatch, min.mismatch, with.indels, fixed,
                       algorithm, verbose, matches.as)
    else
        ans <- .match.XStringSet(pdict, subject,
                       max.mismatch, min.mismatch, with.indels, fixed,
                       algorithm, verbose, matches.as)
    if (is.null(which_pp_excluded))
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

### 'pdict' is an TB_PDict object.
.vmatch.TB_PDict <- function(pdict, subject,
                max.mismatch, min.mismatch, with.indels, fixed,
                algorithm, collapse, weight,
                verbose, matches.as)
{
    .vmatch.PDict3Parts.XStringSet(pdict@threeparts, subject,
                    max.mismatch, min.mismatch, with.indels, fixed,
                    algorithm, collapse, weight,
                    matches.as, NULL)
}

### 'pdict' is an MTB_PDict object.
.vmatch.MTB_PDict <- function(pdict, subject,
                max.mismatch, min.mismatch, with.indels, fixed,
                algorithm, collapse, weight,
                verbose, matches.as)
{
    tb_pdicts <- as.list(pdict)
    NTB <- length(tb_pdicts)
    .checkMaxMismatch(max.mismatch, NTB)
    if (matches.as != "MATCHES_AS_WHICH")
        stop("MTB_PDict objects are not supported yet, sorry")
    ans_compons <- lapply(seq_len(NTB),
        function(i)
        {
            if (verbose)
                cat("Getting results for TB_PDict component ",
                    i, "/", NTB, " ...\n", sep="")
            tb_pdict <- tb_pdicts[[i]]
            .vmatch.TB_PDict(tb_pdict, subject,
                             max.mismatch, min.mismatch, with.indels, fixed,
                             algorithm, collapse, weight,
                             verbose, matches.as)
        }
    )
    if (verbose)
        cat("Combining the results obtained for ",
            "each TB_PDict component...\n", sep="")
    .combine.vwhich.compons(ans_compons)
}

### 'pattern' is an XStringSet object.
.vmatch.XStringSet <- function(pattern, subject,
                max.mismatch, min.mismatch, with.indels, fixed,
                algorithm, collapse, weight,
                verbose, matches.as)
{
    .vmatch.XStringSet.XStringSet(pattern, subject,
                    max.mismatch, min.mismatch, with.indels, fixed,
                    algorithm, collapse, weight,
                    matches.as, NULL)
}

.vmatchPDict <- function(pdict, subject,
                         max.mismatch, min.mismatch, with.indels, fixed,
                         algorithm, collapse, weight,
                         verbose, matches.as="MATCHES_AS_ENDS")
{
    which_pp_excluded <- NULL
    if (is(pdict, "PDict")) {
        if (seqtype(subject) != "DNA")
            stop("'subject' must be DNA")
        dups0 <- dups(pdict)
        if (!is.null(dups0))
            which_pp_excluded <- which(duplicated(dups0))
    } else if (is(pdict, "XStringSet")) {
        if (seqtype(pdict) != seqtype(subject))
            stop("'pdict' and 'subject' must contain ",
                 "sequences of the same type")
    } else {
        stop("'pdict' must be a PDict or XStringSet object")
    }
    max.mismatch <- normargMaxMismatch(max.mismatch)
    min.mismatch <- normargMinMismatch(min.mismatch, max.mismatch)
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    if (matches.as == "MATCHES_AS_WHICH") {
        ## vwhichPDict()
    } else if (matches.as == "MATCHES_AS_COUNTS") {
        ## vcountPDict()
        collapse <- normargCollapse(collapse)
        if (collapse) {
            if (collapse == 1L) {
                weight <- recycleNumericArg(weight, "weight", length(subject))
            } else {
                weight <- recycleNumericArg(weight, "weight", length(pdict))
                if (!is.null(which_pp_excluded)) {
                    ## Collapse weights of duplicates.
                    ## TODO: Implement this in C.
                    which_is_not_unique <-
                        which(!sapply(low2high(dups0), is.null))
                    if (length(which_is_not_unique) != 0L) {
                        weight[which_is_not_unique] <-
                            weight[which_is_not_unique] +
                            sapply(which_is_not_unique,
                                function(i) sum(weight[low2high(dups0)[[i]]]))
                    }
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
    ## We are doing our own dispatch here, based on the type of 'pdict'.
    ## TODO: Revisit this. Would probably be a better design to use a
    ## generic/methods approach and rely on the standard dispatch mechanism.
    ## Those low-level generic/methods wouldn't need to be exported.
    if (is(pdict, "TB_PDict"))
        ans <- .vmatch.TB_PDict(pdict, subject,
                       max.mismatch, min.mismatch, with.indels, fixed,
                       algorithm, collapse, weight,
                       verbose, matches.as)
    else if (is(pdict, "MTB_PDict"))
        ans <- .vmatch.MTB_PDict(pdict, subject,
                       max.mismatch, min.mismatch, with.indels, fixed,
                       algorithm, collapse, weight,
                       verbose, matches.as)
    else
        ans <- .vmatch.XStringSet(pdict, subject,
                       max.mismatch, min.mismatch, with.indels, fixed,
                       algorithm, collapse, weight,
                       verbose, matches.as)
    if (is.null(which_pp_excluded))
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
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        standardGeneric("matchPDict")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPDict", "XString",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        .matchPDict(pdict, subject,
                    max.mismatch, min.mismatch, with.indels, fixed,
                    algorithm, verbose)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPDict", "XStringSet",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        stop("please use vmatchPDict() when 'subject' is an XStringSet ",
             "object (multiple sequence)")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPDict", "XStringViews",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        .matchPDict(pdict, subject,
                    max.mismatch, min.mismatch, with.indels, fixed,
                    algorithm, verbose)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPDict", "MaskedXString",
    function(pdict, subject, 
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        matchPDict(pdict, toXStringViewsOrXString(subject),
                   max.mismatch=max.mismatch, min.mismatch=min.mismatch,
                   with.indels=with.indels, fixed=fixed,
                   algorithm=algorithm, verbose=verbose)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "countPDict" generic and methods.
###

setGeneric("countPDict", signature="subject",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        standardGeneric("countPDict")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPDict", "XString",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        .matchPDict(pdict, subject,
                    max.mismatch, min.mismatch, with.indels, fixed,
                    algorithm, verbose, matches.as="MATCHES_AS_COUNTS")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPDict", "XStringSet",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        stop("please use vcountPDict() when 'subject' is an XStringSet ",
             "object (multiple sequence)")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPDict", "XStringViews",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        .matchPDict(pdict, subject,
                    max.mismatch, min.mismatch, with.indels, fixed,
                    algorithm, verbose, matches.as="MATCHES_AS_COUNTS")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPDict", "MaskedXString",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        countPDict(pdict, toXStringViewsOrXString(subject),
                   max.mismatch=max.mismatch, min.mismatch=min.mismatch,
                   with.indels=with.indels, fixed=fixed,
                   algorithm=algorithm, verbose=verbose)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "whichPDict" generic and methods.
###

setGeneric("whichPDict", signature="subject",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        standardGeneric("whichPDict")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("whichPDict", "XString",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        .matchPDict(pdict, subject,
                    max.mismatch, min.mismatch, with.indels, fixed,
                    algorithm, verbose, matches.as="MATCHES_AS_WHICH")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("whichPDict", "XStringSet",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        stop("please use vwhichPDict() when 'subject' is an XStringSet ",
             "object (multiple sequence)")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("whichPDict", "XStringViews",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        .matchPDict(pdict, subject,
                    max.mismatch, min.mismatch, with.indels, fixed,
                    algorithm, verbose, matches.as="MATCHES_AS_WHICH")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("whichPDict", "MaskedXString",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        whichPDict(pdict, toXStringViewsOrXString(subject),
                   max.mismatch=max.mismatch, min.mismatch=min.mismatch,
                   with.indels=with.indels, fixed=fixed,
                   algorithm=algorithm, verbose=verbose)
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
###   o vmatchPDict(): not supported yet! (first we need a container to
###     store the results)
###

setGeneric("vmatchPDict", signature="subject",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE, ...)
        standardGeneric("vmatchPDict")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vmatchPDict", "ANY",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        stop("vmatchPDict() is not ready yet, sorry")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vmatchPDict", "XString",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        stop("please use matchPDict() when 'subject' is an XString ",
             "object (single sequence)")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vmatchPDict", "MaskedXString",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        stop("please use matchPDict() when 'subject' is a MaskedXString ",
             "object (single sequence)")
)

setGeneric("vcountPDict", signature="subject",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", collapse=FALSE, weight=1L, verbose=FALSE, ...)
        standardGeneric("vcountPDict")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vcountPDict", "XString",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", collapse=FALSE, weight=1L, verbose=FALSE)
        stop("please use countPDict() when 'subject' is an XString ",
             "object (single sequence)")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vcountPDict", "XStringSet",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", collapse=FALSE, weight=1L, verbose=FALSE)
        .vmatchPDict(pdict, subject,
                     max.mismatch, min.mismatch, with.indels, fixed,
                     algorithm, collapse, weight,
                     verbose, matches.as="MATCHES_AS_COUNTS")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vcountPDict", "XStringViews",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", collapse=FALSE, weight=1L, verbose=FALSE)
        vcountPDict(pdict, fromXStringViewsToStringSet(subject),
                    max.mismatch=max.mismatch, min.mismatch=min.mismatch,
                    with.indels=with.indels, fixed=fixed,
                    algorithm=algorithm, collapse=collapse, weight=weight,
                    verbose=verbose)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vcountPDict", "MaskedXString",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", collapse=FALSE, weight=1L, verbose=FALSE)
        stop("please use countPDict() when 'subject' is a MaskedXString ",
             "object (single sequence)")
)

setGeneric("vwhichPDict", signature="subject",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        standardGeneric("vwhichPDict")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vwhichPDict", "XString",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        stop("please use whichPDict() when 'subject' is an XString ",
             "object (single sequence)")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vwhichPDict", "XStringSet",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        .vmatchPDict(pdict, subject,
                     max.mismatch, min.mismatch, with.indels, fixed,
                     algorithm, 0L, 1L,
                     verbose, matches.as="MATCHES_AS_WHICH")
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vwhichPDict", "XStringViews",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        vwhichPDict(pdict, fromXStringViewsToStringSet(subject),
                    max.mismatch=max.mismatch, min.mismatch=min.mismatch,
                    with.indels=with.indels, fixed=fixed,
                    algorithm=algorithm, verbose=verbose)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("vwhichPDict", "MaskedXString",
    function(pdict, subject,
             max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE,
             algorithm="auto", verbose=FALSE)
        stop("please use whichPDict() when 'subject' is a MaskedXString ",
             "object (single sequence)")
)

