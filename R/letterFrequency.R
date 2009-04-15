### =========================================================================
### alphabetFrequency()
### strrev()
### mkAllStrings()
### oligonucleotideFrequency()
### dinucleotideFrequency()
### trinucleotideFrequency()
### consensusMatrix() and consensusString()
### -------------------------------------------------------------------------


.normargFreq <- function(freq)
{
    if (!isTRUEorFALSE(freq))
        stop("'freq' must be TRUE or FALSE")
    freq
}

.set.collapse.default <- function(collapse.default)
{
    assign("collapse.default", collapse.default, envir=RTobjs)
}

.set.collapse.default(FALSE)

.normargCollapse <- function(collapse)
{
    if (is.null(collapse))
        return(get("collapse.default", envir=RTobjs))
    if (!isTRUEorFALSE(collapse))
        stop("'collapse' must be TRUE or FALSE")
    collapse
}

.normargWidth <- function(width)
{
    if (!is.numeric(width) || length(width) != 1 || is.na(width))
        stop("'width' must be a single integer")
    width <- as.integer(width)
    if (width < 0)
        stop("'width' must be a non-negative integer")
    width
}

.normargFastMovingSide <- function(fast.moving.side)
{
    if (!is.character(fast.moving.side)
     || length(fast.moving.side) != 1
     || is.na(fast.moving.side))
        stop("'fast.moving.side' must be a single string")
    match.arg(fast.moving.side, c("left", "right"))
}

.normargAs <- function(as, fast.moving.side)
{
    if (!isSingleString(as))
        stop("'as' must be a single string")
    as <- match.arg(as, c("atomic", "list", "array"))
    if (as == "array" && fast.moving.side != "left")
        stop("'fast.moving.side' must be \"left\" when 'as' is \"array\"")
    as
}

.normargWithLabels <- function(with.labels)
{
    if (!isTRUEorFALSE(with.labels))
        stop("'with.labels' must be TRUE or FALSE")
    with.labels
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "alphabetFrequency" generic and methods.
###
### sum(alphabetFrequency(x)) should always be exactly nchar(x)
###

.XString.letter_frequency <- function(x, freq)
{
    freq <- .normargFreq(freq)
    ans <- .Call("XString_letter_frequency",
                 x, NULL, FALSE,
                 PACKAGE="Biostrings")
    if (freq)
        ans <- ans / nchar(x) # nchar(x) is sum(ans) in faster
    ans
}

.XString.code_frequency <- function(x, baseOnly, freq)
{
    codes <- xscodes(x, baseOnly=baseOnly)
    freq <- .normargFreq(freq)
    ans <- .Call("XString_letter_frequency",
                 x, codes, baseOnly,
                 PACKAGE="Biostrings")
    if (freq)
        ans <- ans / nchar(x) # nchar(x) is sum(ans) in faster
    ans
}

.XStringSet.letter_frequency <- function(x, freq, ...)
{
    collapse <- .normargCollapse(list(...)$collapse)
    ## NO, we cannot use this shortcut when 'collapse' is TRUE because
    ## there is no guarantee that the elements in x cover super(x) entirely
    ## and don't overlap!
    #if (collapse)
    #    return(.XString.letter_frequency(super(x), freq))
    freq <- .normargFreq(freq)
    ans <- .Call("XStringSet_letter_frequency",
                 x, NULL, FALSE, collapse,
                 PACKAGE="Biostrings")
    if (freq) {
        if (collapse)
            ans <- ans / sum(ans)
        else
            ans <- ans / nchar(x)
    }
    ans
}

.XStringSet.code_frequency <- function(x, baseOnly, freq, ...)
{
    collapse <- .normargCollapse(list(...)$collapse)
    ## NO, we cannot use this shortcut when 'collapse' is TRUE because
    ## there is no guarantee that the elements in x cover super(x) entirely
    ## and don't overlap!
    #if (collapse)
    #    return(.XString.code_frequency(super(x), baseOnly, freq))
    codes <- xscodes(x, baseOnly=baseOnly)
    freq <- .normargFreq(freq)
    ans <- .Call("XStringSet_letter_frequency",
                 x, codes, baseOnly, collapse,
                 PACKAGE="Biostrings")
    if (freq) {
        if (collapse)
            ans <- ans / sum(ans)
        else
            ans <- ans / nchar(x)
    }
    ans
}

setGeneric("alphabetFrequency", signature="x",
    function(x, baseOnly=FALSE, freq=FALSE, ...)
        standardGeneric("alphabetFrequency")
)

setMethod("alphabetFrequency", "XString",
    function(x, baseOnly=FALSE, freq=FALSE, ...)
    {
        if (!missing(baseOnly))
            warning("'baseOnly' is ignored for a non DNA or RNA sequence")
        .XString.letter_frequency(x, freq)
    }
)

setMethod("alphabetFrequency", "DNAString",
    function(x, baseOnly=FALSE, freq=FALSE, ...)
        .XString.code_frequency(x, baseOnly, freq)
)

setMethod("alphabetFrequency", "RNAString",
    function(x, baseOnly=FALSE, freq=FALSE, ...)
        .XString.code_frequency(x, baseOnly, freq)
)

setMethod("alphabetFrequency", "XStringSet",
    function(x, baseOnly=FALSE, freq=FALSE, ...)
    {
        if (!missing(baseOnly))
            warning("'baseOnly' is ignored for a non DNA or RNA sequence")
        .XStringSet.letter_frequency(x, freq, ...)
    }
)

setMethod("alphabetFrequency", "DNAStringSet",
    function(x, baseOnly=FALSE, freq=FALSE, ...)
        .XStringSet.code_frequency(x, baseOnly, freq, ...)
)

setMethod("alphabetFrequency", "RNAStringSet",
    function(x, baseOnly=FALSE, freq=FALSE, ...)
        .XStringSet.code_frequency(x, baseOnly, freq, ...)
)

### library(drosophila2probe)
### dict0 <- drosophila2probe$sequence
### x <- XStringViews(as.character(dict0[1:2000]), subjectClass="DNAString")
### alphabetFrequency(x, baseOnly=TRUE)
### y <- DNAStringSet(x)
### alphabetFrequency(y, baseOnly=TRUE)
setMethod("alphabetFrequency", "XStringViews",
    function(x, baseOnly=FALSE, freq=FALSE, ...)
    {
        y <- XStringViewsToSet(x, use.names=FALSE, verbose=FALSE)
        alphabetFrequency(y, baseOnly=baseOnly, freq=freq, ...)
    }
)

setMethod("alphabetFrequency", "MaskedXString",
    function(x, baseOnly=FALSE, freq=FALSE, ...)
    {
        y <- toXStringViewsOrXString(x)
        .set.collapse.default(TRUE)
        ans <- alphabetFrequency(y, baseOnly=baseOnly, freq=freq, ...)
        .set.collapse.default(FALSE)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "hasOnlyBaseLetters" generic and methods.
###

setGeneric("hasOnlyBaseLetters", function(x) standardGeneric("hasOnlyBaseLetters"))

setMethod("hasOnlyBaseLetters", "DNAString",
    function(x)
        alphabetFrequency(x, baseOnly=TRUE)[["other"]] == 0
)

setMethod("hasOnlyBaseLetters", "RNAString",
    function(x) hasOnlyBaseLetters(DNAString(x))
)

setMethod("hasOnlyBaseLetters", "DNAStringSet",
    function(x)
        alphabetFrequency(x, baseOnly=TRUE, collapse=TRUE)[["other"]] == 0
)

setMethod("hasOnlyBaseLetters", "RNAStringSet",
    function(x) hasOnlyBaseLetters(DNAStringSet(x))
)

setMethod("hasOnlyBaseLetters", "XStringViews",
    function(x)
    {
        y <- XStringViewsToSet(x, use.names=FALSE, verbose=FALSE)
        hasOnlyBaseLetters(y)
    }
)

setMethod("hasOnlyBaseLetters", "MaskedDNAString",
    function(x)
    {
        y <- toXStringViewsOrXString(x)
        hasOnlyBaseLetters(y)
    }
)

setMethod("hasOnlyBaseLetters", "MaskedRNAString",
    function(x)
    {
        y <- toXStringViewsOrXString(x)
        hasOnlyBaseLetters(y)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "uniqueLetters" generic and methods.
###

.alphabetFrequencyToUniqueLetters <- function(x_af, x_codes)
{
    if (!is.null(names(x_af)))
        return(names(x_af)[x_af != 0])
    if (!identical(x_codes, 0:255))
        stop("Biostrings internal anomaly: cannot infer names of ",
             "vector returned by 'alphabetFrequency(x)'")
    x_codes <- x_codes[x_af != 0]
    if (min(x_codes) == 0)
        warning("'x' contains embedded nuls")
    intToUtf8(x_codes, multiple=TRUE)
}

setGeneric("uniqueLetters", function(x) standardGeneric("uniqueLetters"))

setMethod("uniqueLetters", "XString",
    function(x)
    {
        x_af <- alphabetFrequency(x)
        .alphabetFrequencyToUniqueLetters(x_af, xscodes(x))
    }
)

setMethod("uniqueLetters", "XStringSet",
    function(x)
    {
        x_af <- alphabetFrequency(x, collapse=TRUE)
        .alphabetFrequencyToUniqueLetters(x_af, xscodes(x))
    }
)

setMethod("uniqueLetters", "XStringViews",
    function(x)
    {
        y <- XStringViewsToSet(x, use.names=FALSE, verbose=FALSE)
        uniqueLetters(y)
    }
)

setMethod("uniqueLetters", "MaskedXString",
    function(x)
    {
        y <- toXStringViewsOrXString(x)
        uniqueLetters(y)
    }
)

### We need to be able to map *any* character whose UTF8 code is between 0 and
### 255 to its code, even the nul character.
### 'x' represents the set of characters to map: it must be a vector of 1-letter
### or empty strings, the empty string being used to represent the nul character.
### Typically, 'x' will be what was returned by uniqueLetters().
### For internal use only (not exported).
safeLettersToInt <- function(x, letters.as.names=FALSE)
{
    ii <- which(x == "")
    ans <- utf8ToInt(paste(x, collapse=""))
    for (i in ii)
        ans <- append(ans, 0, after=i-1)
    if (letters.as.names)
        names(ans) <- x
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "strrev" function.
###

strrev <- function(x)
{
    if (length(x) == 0)
        return(x)
    sapply(strsplit(x, NULL, fixed=TRUE),
           function(xx) paste(rev(xx), collapse=""))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "mkAllStrings" function.
###

.mkAllStringsR <- function(alphabet, width)
{
    if (width == 0)
        return("")
    ansR <- .mkAllStringsR(alphabet, width - 1)
    unlist(lapply(alphabet, function(l) paste(l, ansR, sep="")))
}

.mkAllStringsL <- function(alphabet, width)
{
    if (width == 0)
        return("")
    ansL <- .mkAllStringsL(alphabet, width - 1)
    unlist(lapply(alphabet, function(l) paste(ansL, l, sep="")))
}

.mkAllStrings <- function(alphabet, width, fast.moving.side)
{
    if (fast.moving.side == "right")
        .mkAllStringsR(alphabet, width)
    else
        .mkAllStringsL(alphabet, width)
}

mkAllStrings <- function(alphabet, width, fast.moving.side="right")
{
    if (!is.character(alphabet))
        stop("'alphabet' must be a character vector")
    .mkAllStrings(alphabet,
                  .normargWidth(width),
                  .normargFastMovingSide(fast.moving.side))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "oligonucleotideFrequency" generic and methods.
###
### Except for the 'other' element, oligonucleotideFrequency(x, 1L)
### should be the same as alphabetFrequency(x, baseOnly=TRUE).
###
### Something else worth checking:
###   library(BSgenome.Dmelanogaster.UCSC.dm3)
###   chr3R <- Dmelanogaster$chr3R
###   width <- 12
###   dict0 <- mkAllStrings(names(Biostrings:::DNAcodes(TRUE)), width)
###   names(dict0) <- dict0
###   pdict <- PDict(dict0)
###   system.time(c1 <- countPDict(pdict, chr3R))
###   system.time(c2 <- oligonucleotideFrequency(chr3R, width, with.labels=FALSE))
###   identical(c1, c2) # must be TRUE
### Then try for other values of 'width' (1 <= width <= 12).
### Of course oligonucleotideFrequency() is much better: it is >10x faster, does
### not require preprocessing, and uses much less memory.
###

.formatFreqVector <- function(ans, alphabet, width, freq,
                              fast.moving.side, as, with.labels)
{
    if (freq)
        ans <- ans / sum(ans)
    if (as == "array") {
        if (with.labels)
            dimnames <- rep(list(alphabet), each=width)
        else
            dimnames <- NULL
        ans <- array(ans, dim=rep.int(4L, width), dimnames=dimnames)
    } else if (with.labels) {
        names(ans) <- .mkAllStrings(alphabet, width, fast.moving.side)
    }
    ans
}

.formatFreqMatrix <- function(ans, alphabet, width, freq,
                              fast.moving.side, with.labels)
{
    if (freq)
        ans <- ans / rowSums(ans)
    if (with.labels) 
        colnames(ans) <- .mkAllStrings(alphabet, width, fast.moving.side)
    ans
}

.oligonucleotideFrequency <- function(x, width, freq,
                                      fast.moving.side, as, with.labels)
{
    width <- .normargWidth(width)
    freq <- .normargFreq(freq)
    fast.moving.side <- .normargFastMovingSide(fast.moving.side)
    as <- .normargAs(as, fast.moving.side)
    if (as == "list")
        stop("'as' cannot be \"list\" when 'x' is a DNAString or RNAString object")
    with.labels <- .normargWithLabels(with.labels)
    base_codes <- xscodes(x, baseOnly=TRUE)
    ans <- .Call("XString_oligonucleotide_frequency",
                 x, base_codes, width, fast.moving.side,
                 PACKAGE="Biostrings")
    .formatFreqVector(ans, names(base_codes), width, freq,
                      fast.moving.side, as, with.labels)
}

setGeneric("oligonucleotideFrequency", signature="x",
    function(x, width, freq=FALSE,
             fast.moving.side="right", as="atomic", with.labels=TRUE, ...)
        standardGeneric("oligonucleotideFrequency")
)

setMethod("oligonucleotideFrequency", "DNAString",
    function(x, width, freq=FALSE,
             fast.moving.side="right", as="atomic", with.labels=TRUE, ...)
    {
        .oligonucleotideFrequency(x, width, freq,
                                  fast.moving.side, as, with.labels)
    }
)

setMethod("oligonucleotideFrequency", "RNAString",
    function(x, width, freq=FALSE,
             fast.moving.side="right", as="atomic", with.labels=TRUE, ...)
    {
        .oligonucleotideFrequency(x, width, freq,
                                  fast.moving.side, as, with.labels)
    }
)

setMethod("oligonucleotideFrequency", "XStringSet",
    function(x, width, freq=FALSE,
             fast.moving.side="right", as="atomic", with.labels=TRUE, ...)
    {
        if (!(xsbasetype(x) %in% c("DNA", "RNA")))
            stop("'x' must be of base type DNA or RNA")
        width <- .normargWidth(width)
        freq <- .normargFreq(freq)
        fast.moving.side <- .normargFastMovingSide(fast.moving.side)
        as <- .normargAs(as, fast.moving.side)
        with.labels <- .normargWithLabels(with.labels)
        collapse <- .normargCollapse(list(...)$collapse)
        if (collapse && as == "list")
            stop("'as' cannot be \"list\" when 'collapse' is TRUE")
        base_codes <- xscodes(x, baseOnly=TRUE)
        if (as != "atomic" && !collapse) {
            if (as == "list")
                as <- "atomic"
            ans <- rep.int(list(NULL), length(x))
            for (i in seq_len(length(x))) {
                xx <- x[[i]]
                tmp <- .Call("XString_oligonucleotide_frequency",
                             xx, base_codes, width, fast.moving.side,
                             PACKAGE="Biostrings")
                ans[[i]] <- .formatFreqVector(tmp, names(base_codes), width,
                                              freq, fast.moving.side, as,
                                              with.labels)
            }
            return(ans)
        }
        C_ans <- .Call("XStringSet_oligonucleotide_frequency",
                       x, base_codes, width, fast.moving.side, collapse,
                       PACKAGE="Biostrings")
        if (collapse)
            ans <- .formatFreqVector(C_ans, names(base_codes), width, freq,
                                     fast.moving.side, as, with.labels)
        else
            ans <- .formatFreqMatrix(C_ans, names(base_codes), width, freq,
                                     fast.moving.side, with.labels)
        ans
    }
)

setMethod("oligonucleotideFrequency", "XStringViews",
    function(x, width, freq=FALSE,
             fast.moving.side="right", as="atomic", with.labels=TRUE, ...)
    {
        y <- XStringViewsToSet(x, use.names=FALSE, verbose=FALSE)
        oligonucleotideFrequency(y, width, freq=freq,
                                 fast.moving.side=fast.moving.side,
                                 as=as, with.labels=with.labels, ...)
    }
)

setMethod("oligonucleotideFrequency", "MaskedXString",
    function(x, width, freq=FALSE,
             fast.moving.side="right", as="atomic", with.labels=TRUE, ...)
    {
        y <- toXStringViewsOrXString(x)
        .set.collapse.default(TRUE)
        ans <- oligonucleotideFrequency(y, width, freq=freq,
                                        fast.moving.side=fast.moving.side,
                                        as=as, with.labels=with.labels, ...)
        .set.collapse.default(FALSE)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "dinucleotideFrequency", "trinucleotideFrequency", and
### "oligonucleotideTransitions" convenience wrappers.
###

dinucleotideFrequency <- function(x, freq=FALSE,
                                  fast.moving.side="right",
                                  as.matrix=FALSE, with.labels=TRUE, ...)
{
    if (isTRUE(as.matrix))
        as <- "array"
    else
        as <- "atomic"
    oligonucleotideFrequency(x, 2, freq=freq,
                             fast.moving.side=fast.moving.side,
                             as=as, with.labels=with.labels, ...)
}

trinucleotideFrequency <- function(x, freq=FALSE,
                                   fast.moving.side="right",
                                   as.array=FALSE, with.labels=TRUE, ...)
{
    if (isTRUE(as.array))
        as <- "array"
    else
        as <- "atomic"
    oligonucleotideFrequency(x, 3, freq=freq,
                             fast.moving.side=fast.moving.side,
                             as=as,
                             with.labels=with.labels, ...)
}

oligonucleotideTransitions <- function(x, left=1, right=1, freq=FALSE)
{
    frequencies <- oligonucleotideFrequency(x, width = left + right, freq = freq)
    transitions <-
      matrix(frequencies, nrow = 4 ^ left, ncol = 4 ^ right, byrow = TRUE,
             dimnames = list(unique(substring(names(frequencies), 1, left)),
                             unique(substring(names(frequencies), left + 1, left + right))))
    if (freq)
        transitions <- transitions / rowSums(transitions)
    transitions
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### consensusMatrix() and consensusString().
###

setGeneric("consensusMatrix", function(x, ...) standardGeneric("consensusMatrix"))

setMethod("consensusMatrix", "character",
    function(x, freq=FALSE)
        consensusMatrix(BStringSet(x), freq=freq)
)

setMethod("consensusMatrix", "matrix",
    function(x, freq=FALSE)
        consensusMatrix(BStringSet(apply(x, 1, paste, collapse="")), freq=freq)
)

### 'x' must be a list of FASTA records as one returned by readFASTA()
setMethod("consensusMatrix", "list",
    function(x, freq=FALSE)
        consensusMatrix(BStringSet(FASTArecordsToCharacter(x, use.names=FALSE)), freq=freq)
)

setMethod("consensusMatrix", "XStringSet",
    function(x, baseOnly=FALSE, freq=FALSE, shift=0L, width=NULL)
    {
        if (!is.integer(shift))
            shift <- as.integer(shift)
        if (length(x) != 0 && length(shift) > length(x))
            stop("'shift' has more elements than 'x'")
        if (!is.null(width)) {
            if (!isSingleNumber(width) || width < 0)
                stop("'width' must be NULL or a single non-negative integer")
            if (!is.integer(width))
                width <- as.integer(width)
        }
        codes <- xscodes(x, baseOnly=baseOnly)
        if (is.null(names(codes))) {
            names(codes) <- intToUtf8(codes, multiple = TRUE)
            removeUnused <- TRUE
        } else {
            removeUnused <- FALSE
        }
        freq <- .normargFreq(freq)
        ans <- .Call("XStringSet_letter_frequency_by_pos",
                     x, codes, baseOnly, shift, width,
                     PACKAGE="Biostrings")
        if (removeUnused) {
            ans <- ans[rowSums(ans) > 0, , drop=FALSE]
        }
        if (freq) {
            col_sums <- colSums(ans)
            col_sums[col_sums == 0] <- 1  # to avoid division by 0
            ans <- ans / rep(col_sums, each=nrow(ans))
        }
        ans
    }
)

setMethod("consensusMatrix", "XStringViews",
    function(x, baseOnly=FALSE, freq=FALSE, shift=0L, width=NULL)
    {
        y <- XStringViewsToSet(x, use.names=FALSE, verbose=FALSE)
        consensusMatrix(y, baseOnly=baseOnly, freq=freq, shift=shift, width=width)
    }
)

setGeneric("consensusString", function(x, ...) standardGeneric("consensusString"))

setMethod("consensusString", "matrix",
    function(x)
    {
        err_msg <- c("Please make sure 'x' was obtained by a ",
                     "call to consensusMatrix(..., freq=TRUE)")
        all_letters <- rownames(x)
        if (is.null(all_letters))
            stop("invalid consensus matrix 'x' (has no row names).\n",
                 "  ", err_msg)
        if (!all(nchar(all_letters) == 1))
            stop("invalid consensus matrix 'x' (row names must be single letters).\n",
                 "  ", err_msg)
        if (is.integer(x)) {
            col_sums <- colSums(x)
            col_sums[col_sums == 0] <- 1  # to avoid division by 0
            x <- x / rep(col_sums, each=nrow(x))
        }
        consensusLetter <- function(col)
        {
            i <- which.max(col)
            if (length(i) == 0L)  # 'i' is an integer(0)
                stop("invalid consensus matrix 'x' (has 0 rows or contains NAs/NaNs).\n",
                     "  ", err_msg)
            if (col[i] > 1)
                stop("invalid consensus matrix 'x' (contains values > 1).\n",
                     "  ", err_msg)
            if (col[i] > 0.5)
                return(all_letters[i])
            return("?")
        }
        paste(apply(x, 2, consensusLetter), collapse="")
    }
)

setMethod("consensusString", "XStringSet",
    function(x, shift=0L, width=NULL)
        consensusString(consensusMatrix(x, freq=TRUE, shift=shift, width=width))
)

setMethod("consensusString", "XStringViews",
    function(x, shift=0L, width=NULL)
        consensusString(consensusMatrix(x, freq=TRUE, shift=shift, width=width))
)

setMethod("consensusString", "ANY",
    function(x) consensusString(consensusMatrix(x, freq=TRUE))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Old stuff (Defunct or Deprecated).
###

setGeneric("consmat", function(x, ...) standardGeneric("consmat"))

setMethod("consmat", "ANY",
    function(x, ...)
    {
        .Deprecated("consensusMatrix")
        consensusMatrix(x, ...)
    }
)

