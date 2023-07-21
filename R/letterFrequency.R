### =========================================================================
### alphabetFrequency()
### hasOnlyBaseLetters()
### uniqueLetters()
### letterFrequencyInSlidingView()
### letterFrequency()
### mkAllStrings()
### oligonucleotideFrequency()
### dinucleotideFrequency()
### trinucleotideFrequency()
### oligonucleotideTransitions()
### nucleotideFrequencyAt()
### consensusMatrix() and consensusString()
### -------------------------------------------------------------------------

.normargWidth <- function(width, argname="width", min.width=1L)
{
    if (!isSingleNumber(width))
        stop("'", argname, "' must be a single integer")
    width <- as.integer(width)
    if (width < min.width)
        stop("'", argname, "' must be an integer >= ", min.width)
    width
}

.normargStep <- function(step)
{
    if (!isSingleNumber(step))
        stop("'step' must be a single integer")
    step <- as.integer(step)
    if (step <= 0L)
        stop("'step' must be a positive integer")
    step
}

.normargCollapse <- function(collapse)
{
    if (!isTRUEorFALSE(collapse))
        stop("'collapse' must be TRUE or FALSE")
    collapse
}

.normargAsArray <- function(as.array)
{
    if (!isTRUEorFALSE(as.array))
        stop("'as.array' must be TRUE or FALSE")
    as.array
}

.normargFastMovingSide <- function(fast.moving.side, as.array=FALSE)
{
    if (as.array)
        return("left")
    if (!isSingleString(fast.moving.side))
        stop("'fast.moving.side' must be a single string")
    match.arg(fast.moving.side, c("left", "right"))
}

.normargWithLabels <- function(with.labels)
{
    if (!isTRUEorFALSE(with.labels))
        stop("'with.labels' must be TRUE or FALSE")
    with.labels
}

.normargSimplifyAs <- function(simplify.as, as.array)
{
    if (!isSingleString(simplify.as))
        stop("'simplify.as' must be a single string")
    simplify.as <- match.arg(simplify.as, c("matrix", "list", "collapsed"))
    if (simplify.as == "matrix" && as.array)
        stop("'as.array' cannot be TRUE when 'simplify.as' is \"matrix\"")
    simplify.as
}

### Author: HJ
.normargLetters <- function(letters, alphabet)
{
    if (!is.character(letters) || any(is.na(letters)))
        stop("'letters' must be a character vector with no NAs")
    if (any(nchar(letters) == 0L))
        stop("'letters' cannot contain empty strings")
    single_letters <- unlist(strsplit(letters, NULL, fixed=TRUE))
    if (any(duplicated(single_letters)))
        stop("letters in 'letters' must be unique")
    if (!is.null(alphabet)) {
        is_valid_letter <- single_letters %in% alphabet
        if (!all(is_valid_letter))
            stop("invalid letter(s): ",
                 paste(single_letters[!is_valid_letter], collapse=", "))
    }
    single_letters
}

.normargOR <- function(OR)
{
    if (!isSingleString(OR) && !(isSingleNumber(OR) && OR == 0))
        stop("'OR' must be a single string or 0")
    OR
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "alphabetFrequency" generic and methods.
###
### sum(alphabetFrequency(x)) should always be exactly nchar(x)
###

.XString.letter_frequency <- function(x, as.prob)
{
    if (!isTRUEorFALSE(as.prob))
        stop("'as.prob' must be TRUE or FALSE")
    ans <- .Call2("XString_letter_frequency",
                 x, NULL, FALSE,
                 PACKAGE="Biostrings")
    if (as.prob)
        ans <- ans / nchar(x) # nchar(x) is sum(ans) but faster
    ans
}

.XString.nucleotide_frequency <- function(x, as.prob, baseOnly)
{
    if (!isTRUEorFALSE(as.prob))
        stop("'as.prob' must be TRUE or FALSE")
    codes <- xscodes(x, baseOnly=baseOnly)
    ans <- .Call2("XString_letter_frequency",
                 x, codes, baseOnly,
                 PACKAGE="Biostrings")
    if (as.prob)
        ans <- ans / nchar(x) # nchar(x) is sum(ans) but faster
    ans
}

.XString.amino_acid_frequency <- function(x, as.prob, baseOnly)
{
    if (!isTRUEorFALSE(as.prob))
        stop("'as.prob' must be TRUE or FALSE")
    codes <- xscodes(x, baseOnly=baseOnly)
    ans <- .Call2("XString_letter_frequency",
                 x, codes, baseOnly,
                 PACKAGE="Biostrings")
    if (as.prob)
        ans <- ans / nchar(x) # nchar(x) is sum(ans) but faster
    ans
}

.XStringSet.letter_frequency <- function(x, as.prob, collapse)
{
    if (!isTRUEorFALSE(as.prob))
        stop("'as.prob' must be TRUE or FALSE")
    collapse <- .normargCollapse(collapse)
    ans <- .Call2("XStringSet_letter_frequency",
                 x, collapse, NULL, FALSE,
                 PACKAGE="Biostrings")
    if (as.prob) {
        if (collapse)
            ans <- ans / sum(ans)
        else
            ans <- ans / nchar(x)
    }
    ans
}

.XStringSet.nucleotide_frequency <- function(x, as.prob, collapse, baseOnly)
{
    if (!isTRUEorFALSE(as.prob))
        stop("'as.prob' must be TRUE or FALSE")
    collapse <- .normargCollapse(collapse)
    codes <- xscodes(x, baseOnly=baseOnly)
    ans <- .Call2("XStringSet_letter_frequency",
                 x, collapse, codes, baseOnly,
                 PACKAGE="Biostrings")
    if (as.prob) {
        if (collapse)
            ans <- ans / sum(ans)
        else
            ans <- ans / nchar(x)
    }
    ans
}

.XStringSet.amino_acid_frequency <- function(x, as.prob, collapse, baseOnly)
{
    if (!isTRUEorFALSE(as.prob))
        stop("'as.prob' must be TRUE or FALSE")
    collapse <- .normargCollapse(collapse)
    codes <- xscodes(x, baseOnly=baseOnly)
    ans <- .Call2("XStringSet_letter_frequency",
                 x, collapse, codes, baseOnly,
                 PACKAGE="Biostrings")
    if (as.prob) {
        if (collapse)
            ans <- ans / sum(ans)
        else
            ans <- ans / nchar(x)
    }
    ans
}

setGeneric("alphabetFrequency", signature="x",
    function(x, as.prob=FALSE, ...) standardGeneric("alphabetFrequency")
)

setMethod("alphabetFrequency", "XString",
    function(x, as.prob=FALSE)
        .XString.letter_frequency(x, as.prob)
)

setMethod("alphabetFrequency", "DNAString",
    function(x, as.prob=FALSE, baseOnly=FALSE)
        .XString.nucleotide_frequency(x, as.prob, baseOnly)
)

setMethod("alphabetFrequency", "RNAString",
    function(x, as.prob=FALSE, baseOnly=FALSE)
        .XString.nucleotide_frequency(x, as.prob, baseOnly)
)

setMethod("alphabetFrequency", "AAString",
    function(x, as.prob=FALSE, baseOnly=FALSE)
        .XString.amino_acid_frequency(x, as.prob, baseOnly)
)

setMethod("alphabetFrequency", "XStringSet",
    function(x, as.prob=FALSE, collapse=FALSE)
        .XStringSet.letter_frequency(x, as.prob, collapse)
)

setMethod("alphabetFrequency", "DNAStringSet",
    function(x, as.prob=FALSE, collapse=FALSE, baseOnly=FALSE)
        .XStringSet.nucleotide_frequency(x, as.prob, collapse, baseOnly)
)

setMethod("alphabetFrequency", "RNAStringSet",
    function(x, as.prob=FALSE, collapse=FALSE, baseOnly=FALSE)
        .XStringSet.nucleotide_frequency(x, as.prob, collapse, baseOnly)
)

setMethod("alphabetFrequency", "AAStringSet",
    function(x, as.prob=FALSE, collapse=FALSE, baseOnly=FALSE)
        .XStringSet.amino_acid_frequency(x, as.prob, collapse, baseOnly)
)

### library(drosophila2probe)
### dict0 <- drosophila2probe$sequence
### x <- DNAStringSet(dict0[1:2000])
### alphabetFrequency(x, baseOnly=TRUE)
### y <- as(x, "Views")
### alphabetFrequency(y, baseOnly=TRUE)
setMethod("alphabetFrequency", "XStringViews",
    function(x, as.prob=FALSE, ...)
    {
        y <- fromXStringViewsToStringSet(x)
        alphabetFrequency(y, as.prob=as.prob, ...)
    }
)

setMethod("alphabetFrequency", "MaskedXString",
    function(x, as.prob=FALSE, ...)
    {
        y <- as(x, "XStringViews")
        alphabetFrequency(y, as.prob=as.prob, collapse=TRUE, ...)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "hasOnlyBaseLetters" generic and methods.
###

setGeneric("hasOnlyBaseLetters", function(x) standardGeneric("hasOnlyBaseLetters"))

setMethod("hasOnlyBaseLetters", "DNAString",
    function(x)
        alphabetFrequency(x, baseOnly=TRUE)[["other"]] == 0L
)

setMethod("hasOnlyBaseLetters", "RNAString",
    function(x) hasOnlyBaseLetters(DNAString(x))
)

setMethod("hasOnlyBaseLetters", "AAString",
    function(x)
        alphabetFrequency(x, baseOnly=TRUE)[["other"]] == 0L
)

setMethod("hasOnlyBaseLetters", "DNAStringSet",
    function(x)
        alphabetFrequency(x, collapse=TRUE, baseOnly=TRUE)[["other"]] == 0L
)

setMethod("hasOnlyBaseLetters", "RNAStringSet",
    function(x) hasOnlyBaseLetters(DNAStringSet(x))
)

setMethod("hasOnlyBaseLetters", "AAStringSet",
    function(x)
        alphabetFrequency(x, collapse=TRUE, baseOnly=TRUE)[["other"]] == 0L
)

setMethod("hasOnlyBaseLetters", "XStringViews",
    function(x)
    {
        y <- fromXStringViewsToStringSet(x)
        hasOnlyBaseLetters(y)
    }
)

setMethod("hasOnlyBaseLetters", "MaskedXString",
    function(x)
    {
        y <- as(x, "XStringViews")
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
    if (length(x_codes) != 0 && min(x_codes) == 0)
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
        y <- fromXStringViewsToStringSet(x)
        uniqueLetters(y)
    }
)

setMethod("uniqueLetters", "MaskedXString",
    function(x)
    {
        y <- as(x, "XStringViews")
        uniqueLetters(y)
    }
)

### We need to be able to map *any* character whose UTF8 code is between 0 and
### 255 to its code, even the nul character.
### 'x' represents the set of characters to map: it must be a vector of
### 1-letter or empty strings, the empty string being used to represent the nul
### character. Typically, 'x' will be what was returned by uniqueLetters().
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
### The "letterFrequency" and "letterFrequencyInSlidingView"
### generics and methods.
### Author: HJ
###

### joint C interface
### The value is a matrix for letterFrequencyInSlidingView
### and a matrix for letterFrequency unless collapse=TRUE.
.letterFrequency <- function(x, view.width, letters, OR, collapse=FALSE)
{
    ## letterFrequency / letterFrequencyInSlidingView switch
    is_sliding <- !is.na(view.width)
    single_letters <- .normargLetters(letters, alphabet(x))
    OR <- .normargOR(OR)
    codes <- xscodes(x)
    if (is.null(names(codes)))
        single_codes <- as.integer(BString(paste(single_letters, collapse="")))
    else
        single_codes <- codes[single_letters]
    ## Unless 'OR == 0', letters in multi-character elements of
    ## 'letters' are to be grouped (i.e. tabulated in common).
    ## We send a vector indicating the column (1-based) into which each
    ## letter in 'letters' should be tabulated.  For example, for
    ## 'letters = c("CG", "AT")' and 'OR != 0', we send 'c(1,1,2,2)'.
    ## The columns of the result are named accordingly using the OR symbol.
    nc <- nchar(letters)
    if (all(nc == 1L) || OR == 0) {
        colmap <- NULL
        colnames <- single_letters
    } else {
        colmap <- rep.int(seq_len(length(letters)), nc)
        colnames <- sapply(strsplit(letters, NULL, fixed=TRUE),
                           function(z) paste(z, collapse=OR))
    }
    if (is_sliding)
        .Call2("XString_letterFrequencyInSlidingView",
              x, view.width, single_codes, colmap, colnames,
              PACKAGE="Biostrings")
    else
        .Call2("XStringSet_letterFrequency",
              x, single_codes, colmap, colnames, collapse,
              PACKAGE="Biostrings")
}

### letterFrequencyInSlidingView
setGeneric("letterFrequencyInSlidingView", signature="x",
    function(x, view.width, letters, OR="|", as.prob=FALSE)
        standardGeneric("letterFrequencyInSlidingView")
)

### Ensure 'view.width' is not NA
setMethod("letterFrequencyInSlidingView", "XString",
    function(x, view.width, letters, OR="|", as.prob=FALSE)
    {
        view.width <- .normargWidth(view.width, "view.width")
        if (!isTRUEorFALSE(as.prob))
            stop("'as.prob' must be TRUE or FALSE")
        ans <- .letterFrequency(x, view.width, letters=letters, OR=OR)
        if (as.prob)
            ans <- ans / view.width
        ans
    }
)

### letterFrequency
setGeneric("letterFrequency", signature="x",
    function(x, letters, OR="|", as.prob=FALSE, ...)
        standardGeneric("letterFrequency")
)

setMethod("letterFrequency", "XStringSet",
    function(x, letters, OR="|", as.prob=FALSE, collapse=FALSE)
    {
        if (!isTRUEorFALSE(as.prob))
            stop("'as.prob' must be TRUE or FALSE")
        if (!isTRUEorFALSE(collapse))
            stop("'collapse' must be TRUE or FALSE")
        ans <- .letterFrequency(x, NA, letters=letters, OR=OR,
                                collapse=collapse)
        if (as.prob) {
            nc <- nchar(x)
            if (collapse)
                nc <- sum(nc)
            ans <- ans / nc
        }
        ans
    }
)

setMethod("letterFrequency", "XString",
    function(x, letters, OR="|", as.prob=FALSE)
        letterFrequency(as(x, "XStringSet"),
            letters=letters, OR=OR, as.prob=as.prob, collapse=TRUE)
)

setMethod("letterFrequency", "XStringViews",
    function(x, letters, OR="|", as.prob=FALSE, ...)
        letterFrequency(as(x, "XStringSet"),
            letters=letters, OR=OR, as.prob=as.prob, ...)
)

setMethod("letterFrequency", "MaskedXString",
    function(x, letters, OR="|", as.prob=FALSE)
        letterFrequency(as(x, "XStringViews"),
            letters=letters, OR=OR, as.prob=as.prob, collapse=TRUE)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "mkAllStrings" function.
###

.mkAllStringsR <- function(alphabet, width)
{
    if (width == 0L)
        return("")
    ansR <- .mkAllStringsR(alphabet, width - 1L)
    unlist(lapply(alphabet, function(l) paste(l, ansR, sep="")))
}

.mkAllStringsL <- function(alphabet, width)
{
    if (width == 0L)
        return("")
    ansL <- .mkAllStringsL(alphabet, width - 1L)
    unlist(lapply(alphabet, function(l) paste(ansL, l, sep="")))
}

.mkAllStrings <- function(alphabet, width, fast.moving.side)
{
    if (length(alphabet) == 0L && width != 0L)
        return(character(0))
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
                  .normargWidth(width, min.width=0L),
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

setGeneric("oligonucleotideFrequency", signature="x",
    function(x, width, step=1,
             as.prob=FALSE, as.array=FALSE,
             fast.moving.side="right", with.labels=TRUE, ...)
        standardGeneric("oligonucleotideFrequency")
)

setMethod("oligonucleotideFrequency", "XString",
    function(x, width, step=1,
             as.prob=FALSE, as.array=FALSE,
             fast.moving.side="right", with.labels=TRUE)
    {
        if (!(seqtype(x) %in% c("DNA", "RNA")))
            stop("'x' must contain sequences of type DNA or RNA")
        width <- .normargWidth(width)
        step <- .normargStep(step)
        if (!isTRUEorFALSE(as.prob))
            stop("'as.prob' must be TRUE or FALSE")
        as.array <- .normargAsArray(as.array)
        fast.moving.side <- .normargFastMovingSide(fast.moving.side, as.array)
        with.labels <- .normargWithLabels(with.labels)
        base_codes <- xscodes(x, baseOnly=TRUE)
        .Call2("XString_oligo_frequency",
               x, width, step,
               as.prob, as.array,
               fast.moving.side, with.labels,
               base_codes,
               PACKAGE="Biostrings")
    }
)

setMethod("oligonucleotideFrequency", "XStringSet",
    function(x, width, step=1,
             as.prob=FALSE, as.array=FALSE,
             fast.moving.side="right", with.labels=TRUE,
             simplify.as="matrix")
    {
        if (!(seqtype(x) %in% c("DNA", "RNA")))
            stop("'x' must contain sequences of type DNA or RNA")
        width <- .normargWidth(width)
        step <- .normargStep(step)
        if (!isTRUEorFALSE(as.prob))
            stop("'as.prob' must be TRUE or FALSE")
        as.array <- .normargAsArray(as.array)
        fast.moving.side <- .normargFastMovingSide(fast.moving.side, as.array)
        with.labels <- .normargWithLabels(with.labels)
        simplify.as <- .normargSimplifyAs(simplify.as, as.array)
        base_codes <- xscodes(x, baseOnly=TRUE)
        .Call2("XStringSet_oligo_frequency",
               x, width, step,
               as.prob, as.array,
               fast.moving.side, with.labels, simplify.as,
               base_codes,
               PACKAGE="Biostrings")
    }
)

setMethod("oligonucleotideFrequency", "XStringViews",
    function(x, width, step=1,
             as.prob=FALSE, as.array=FALSE,
             fast.moving.side="right", with.labels=TRUE, ...)
    {
        y <- fromXStringViewsToStringSet(x)
        oligonucleotideFrequency(y, width, step=step,
                                 as.prob=as.prob, as.array=as.array,
                                 fast.moving.side=fast.moving.side,
                                 with.labels=with.labels, ...)
    }
)

setMethod("oligonucleotideFrequency", "MaskedXString",
    function(x, width, step=1,
             as.prob=FALSE, as.array=FALSE,
             fast.moving.side="right", with.labels=TRUE)
    {
        y <- as(x, "XStringViews")
        oligonucleotideFrequency(y, width, step=step,
                                 as.prob=as.prob, as.array=as.array,
                                 fast.moving.side=fast.moving.side,
                                 with.labels=with.labels,
                                 simplify.as="collapsed")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "dinucleotideFrequency", "trinucleotideFrequency", and
### "oligonucleotideTransitions" convenience wrappers.
###

dinucleotideFrequency <- function(x, step=1,
                                  as.prob=FALSE, as.matrix=FALSE,
                                  fast.moving.side="right",
                                  with.labels=TRUE, ...)
{
    oligonucleotideFrequency(x, 2L, step=step,
                             as.prob=as.prob, as.array=as.matrix,
                             fast.moving.side=fast.moving.side,
                             with.labels=with.labels, ...)
}

trinucleotideFrequency <- function(x, step=1,
                                   as.prob=FALSE, as.array=FALSE,
                                   fast.moving.side="right",
                                   with.labels=TRUE, ...)
{
    oligonucleotideFrequency(x, 3L, step=step,
                             as.prob=as.prob, as.array=as.array,
                             fast.moving.side=fast.moving.side,
                             with.labels=with.labels, ...)
}

oligonucleotideTransitions <- function(x, left=1, right=1, as.prob=FALSE)
{
    if (!isSingleNumber(left))
        stop("'left' must be a single integer value")
    if (!is.integer(left))
        left <- as.integer(left)
    if (left < 1L)
        stop("'left' must be >= 1")
    if (!isSingleNumber(right))
        stop("'right' must be a single integer value")
    if (!is.integer(right))
        right <- as.integer(right)
    if (right < 1L)
        stop("'right' must be >= 1")
    if (!isTRUEorFALSE(as.prob))
        stop("'as.prob' must be TRUE or FALSE")
    if (is(x, "XStringSet") || is(x, "XStringViews")) {
        freqs <- oligonucleotideFrequency(x, width=left+right, as.prob=FALSE,
                                          simplify.as="collapsed")
    } else {
        freqs <- oligonucleotideFrequency(x, width=left+right, as.prob=FALSE)
    }
    ans_rownames <- mkAllStrings(DNA_BASES, left)
    ans_colnames <- mkAllStrings(DNA_BASES, right)
    ## A sanity check
    expected_freqs_names <- paste(rep(ans_rownames, each=length(ans_colnames)),
                                  ans_colnames, sep="")
    if (!identical(names(freqs), expected_freqs_names))
        stop("Biostrings internal error in oligonucleotideTransitions(): ",
             "'freqs' has unexpected or missing names")
    ans <- matrix(freqs, ncol=length(ans_colnames), byrow=TRUE,
                  dimnames=list(ans_rownames, ans_colnames))
    if (as.prob)
        ans <- ans / rowSums(ans)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "nucleotideFrequencyAt" generic and methods.
###

setGeneric("nucleotideFrequencyAt", signature="x",
    function(x, at, as.prob=FALSE, as.array=TRUE,
             fast.moving.side="right", with.labels=TRUE, ...)
        standardGeneric("nucleotideFrequencyAt")
)

setMethod("nucleotideFrequencyAt", "XStringSet",
    function(x, at, as.prob=FALSE, as.array=TRUE,
             fast.moving.side="right", with.labels=TRUE)
    {
        if (!(seqtype(x) %in% c("DNA", "RNA")))
            stop("'x' must contain sequences of type DNA or RNA")
        if (!is.numeric(at))
            stop("'at' must be a vector of integers")
        if (!is.integer(at))
            at <- as.integer(at)
        if (!isTRUEorFALSE(as.prob))
            stop("'as.prob' must be TRUE or FALSE")
        as.array <- .normargAsArray(as.array)
        fast.moving.side <- .normargFastMovingSide(fast.moving.side, as.array)
        with.labels <- .normargWithLabels(with.labels)
        base_codes <- xscodes(x, baseOnly=TRUE)
        .Call2("XStringSet_nucleotide_frequency_at",
              x, at, as.prob, as.array,
              fast.moving.side, with.labels,
              base_codes,
              PACKAGE="Biostrings")
    }
)

setMethod("nucleotideFrequencyAt", "XStringViews",
    function(x, at, as.prob=FALSE, as.array=TRUE,
             fast.moving.side="right", with.labels=TRUE, ...)
    {
        y <- fromXStringViewsToStringSet(x)
        if (any(width(y) < width(x)))
            stop("x contains \"out of limits\" views")
        nucleotideFrequencyAt(y, at, as.prob=as.prob, as.array=as.array,
                              fast.moving.side=fast.moving.side,
                              with.labels=with.labels, ...)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### consensusMatrix() and consensusString().
###

setGeneric("consensusMatrix", signature="x",
    function(x, as.prob=FALSE, shift=0L, width=NULL, ...)
        standardGeneric("consensusMatrix"))

setMethod("consensusMatrix", "character",
    function(x, as.prob=FALSE, shift=0L, width=NULL)
        consensusMatrix(BStringSet(x),
                        as.prob=as.prob, shift=shift, width=width)
)

setMethod("consensusMatrix", "matrix",
    function(x, as.prob=FALSE, shift=0L, width=NULL)
        consensusMatrix(BStringSet(apply(x, 1, paste, collapse="")),
                        as.prob=as.prob, shift=shift, width=width)
)

setMethod("consensusMatrix", "XStringSet",
    function(x, as.prob=FALSE, shift=0L, width=NULL, baseOnly=FALSE)
    {
        if (!isTRUEorFALSE(as.prob))
            stop("'as.prob' must be TRUE or FALSE")
        if (!is.integer(shift))
            shift <- as.integer(shift)
        if (length(x) != 0 && length(shift) > length(x))
            stop("'shift' has more elements than 'x'")
        if (!is.null(width)) {
            if (!isSingleNumber(width) || width < 0L)
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
        ans <- .Call2("XStringSet_consensus_matrix",
                     x, shift, width, baseOnly, codes,
                     PACKAGE="Biostrings")
        if (removeUnused) {
            ans <- ans[rowSums(ans) > 0, , drop=FALSE]
        }
        if (as.prob) {
            col_sums <- colSums(ans)
            col_sums[col_sums == 0] <- 1  # to avoid division by 0
            ans <- ans / rep(col_sums, each=nrow(ans))
        }
        ans
    }
)

setMethod("consensusMatrix", "XStringViews",
    function(x, as.prob=FALSE, shift=0L, width=NULL, ...)
    {
        y <- fromXStringViewsToStringSet(x)
        consensusMatrix(y, as.prob=as.prob, shift=shift, width=width, ...)
    }
)

setGeneric("consensusString", function(x, ...) standardGeneric("consensusString"))

setMethod("consensusString", "matrix",
    function(x, ambiguityMap="?", threshold=0.5)
    {
        x <- x[rowSums(x, na.rm=TRUE) > 0, , drop=FALSE]
        k <- nrow(x)
        if (k == 0)
            return(character())
        err_msg <- c("Please make sure 'x' was obtained by a ",
                     "call to consensusMatrix(..., as.prob=TRUE)")
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
            x <- x / rep(col_sums, each=k)
        }
        if (any(is.na(x) | x < 0 | x > 1))
            stop("invalid consensus matrix 'x' ",
                 "(contains NAs/NaNs or values outside [0, 1]).\n",
                 "  ", err_msg)
        if (any(abs(colSums(x) - 1) > .Machine$double.eps ^ 0.5))
            stop("invalid consensus matrix 'x' ",
                 "(some columns do not sum to 1).\n",
                 "  ", err_msg)
        if (isSingleString(ambiguityMap)) {
            if (nchar(ambiguityMap) != 1)
                stop("'ambiguityMap' must be a single character or a map ",
                     "(e.g. IUPAC_CODE_MAP)")
            if (!isSingleNumber(threshold) || threshold <= 0 || threshold > 1)
                stop("'threshold' must be a numeric in (0, 1]")
            consensusLetter <- function(col)
            {
                i <- which(col >= threshold)
                if (length(i) == 1 && col[i] >= threshold)
                    all_letters[i]
                else
                    ambiguityMap
            }
        } else {
            if (!is.character(ambiguityMap) || is.null(names(ambiguityMap)))
                stop("'ambiguityMap' must be a named character vector")
            if (!all(rownames(x) %in% names(ambiguityMap)))
                stop("'ambiguityMap' does not contain the complete alphabet")
            alphabet <- unname(ambiguityMap[nchar(ambiguityMap) == 1])
            if (!isSingleNumber(threshold) || threshold <= 0 ||
                (threshold - .Machine$double.eps ^ 0.5) > 1/length(alphabet))
                stop("'threshold' must be a numeric in ",
                     "(0, 1/sum(nchar(ambiguityMap) == 1)]")
            P <-
              sapply(strsplit(ambiguityMap[rownames(x)], ""),
                     function(y) {z <- alphabet %in% y;z/sum(z)})
            x <- P %*% x
            consensusLetter <- function(col)
            {
                i <- paste(alphabet[col >= threshold], collapse = "")
                ans <- names(ambiguityMap[ambiguityMap == i])
                if (length(ans) == 0)
                    stop("'ambiguityMap' is missing some combinations of row names")
                ans
            }
        }
        paste(apply(x, 2, consensusLetter), collapse="")
    }
)

setMethod("consensusString", "BStringSet",
    function(x, ambiguityMap="?", threshold=0.5, shift=0L, width=NULL)
        consensusString(consensusMatrix(x, as.prob=TRUE, shift=shift, width=width),
                        ambiguityMap=ambiguityMap, threshold=threshold)
)

setMethod("consensusString", "DNAStringSet",
    function(x, ambiguityMap=IUPAC_CODE_MAP, threshold=0.25, shift=0L, width=NULL)
        consensusString(consensusMatrix(x, as.prob=TRUE, shift=shift, width=width),
                        ambiguityMap=ambiguityMap, threshold=threshold)
)

setMethod("consensusString", "RNAStringSet",
    function(x,
             ambiguityMap=
             structure(as.character(RNAStringSet(DNAStringSet(IUPAC_CODE_MAP))),
                       names=
                       as.character(RNAStringSet(DNAStringSet(names(IUPAC_CODE_MAP))))),
             threshold=0.25, shift=0L, width=NULL)
        consensusString(consensusMatrix(x, as.prob=TRUE, shift=shift, width=width),
                        ambiguityMap=ambiguityMap, threshold=threshold)
)

setMethod("consensusString", "XStringViews",
    function(x, ambiguityMap, threshold, shift=0L, width=NULL)
    {
        x <- as(x, "XStringSet")
        callGeneric()
    }
)

setMethod("consensusString", "ANY",
    function(x, ambiguityMap="?", threshold=0.5)
        consensusString(consensusMatrix(x, as.prob=TRUE),
                        ambiguityMap=ambiguityMap, threshold=threshold)
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "twoWayAlphabetFrequency" generic and methods.
###
### Calculates two-way table of letter frequencies for two same-length inputs
###

.XString.two_way_code_frequency <- function(x, y, as.prob, baseOnly)
{
  if (!isTRUEorFALSE(as.prob))
    stop("'as.prob' must be TRUE or FALSE")
  x.codes <- xscodes(x, baseOnly=baseOnly)
  y.codes <- xscodes(y, baseOnly=baseOnly)
  ans <- .Call2("XString_two_way_letter_frequency",
                x, y, x.codes, y.codes, baseOnly,
                PACKAGE="Biostrings")
  if (as.prob)
    ans <- ans / nchar(x) # nchar(x) is sum(ans) but faster
  ans
}

.XStringSet.two_way_code_frequency <- function(x, y, as.prob, collapse,
                                               baseOnly)
{
  if (!isTRUEorFALSE(as.prob))
    stop("'as.prob' must be TRUE or FALSE")
  collapse <- .normargCollapse(collapse)
  x.codes <- xscodes(x, baseOnly=baseOnly)
  y.codes <- xscodes(y, baseOnly=baseOnly)
  ans <- .Call2("XStringSet_two_way_letter_frequency",
                x, y, collapse, x.codes, y.codes, baseOnly,
                PACKAGE="Biostrings")
  if (as.prob) {
    if (collapse)
      ans <- ans / sum(as.numeric(ans))
    else
      ans <- ans / rep(nchar(x), each = prod(dim(ans)[1:2]))
  }
  ans
}

setGeneric("twoWayAlphabetFrequency", signature=c("x", "y"),
           function(x, y, as.prob=FALSE, ...)
           standardGeneric("twoWayAlphabetFrequency")
           )

setMethod("twoWayAlphabetFrequency", c("XString", "XString"),
          function(x, y, as.prob=FALSE, baseOnly=FALSE)
          .XString.two_way_code_frequency(x, y, as.prob, baseOnly)
          )

setMethod("twoWayAlphabetFrequency", c("XString", "XStringSet"),
          function(x, y, as.prob=FALSE, collapse=FALSE, baseOnly=FALSE) {
          x <- rep(as(x, "XStringSet"), length(y))
          .XStringSet.two_way_code_frequency(x, y, as.prob, collapse, baseOnly)
          }
          )

setMethod("twoWayAlphabetFrequency", c("XStringSet", "XString"),
          function(x, y, as.prob=FALSE, collapse=FALSE, baseOnly=FALSE) {
          y <- rep(as(y, "XStringSet"), length(x))
          .XStringSet.two_way_code_frequency(x, y, as.prob, collapse, baseOnly)
          }
          )

setMethod("twoWayAlphabetFrequency", c("XStringSet", "XStringSet"),
          function(x, y, as.prob=FALSE, collapse=FALSE, baseOnly=FALSE)
          .XStringSet.two_way_code_frequency(x, y, as.prob, collapse, baseOnly)
          )


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "twoWayAlphabetFrequencyByQuality" generic and methods.
###
### Not sure about this yet. We may want a threeWayAlphabetFrequency,
### where the qualities are combined apriori.
###

.XStringSet.two_way_code_frequency_by_quality <-
  function(x, y, x.quality, y.quality, as.prob, baseOnly)
{
  if (!isTRUEorFALSE(as.prob))
    stop("'as.prob' must be TRUE or FALSE")
  codes <- xscodes(x, baseOnly=baseOnly)
  ans <- .Call2("XStringSet_two_way_letter_frequency_by_quality",
                x, y, x.quality, y.quality, collapse, codes, baseOnly,
                PACKAGE="Biostrings")
  if (as.prob) {
    ans <- ans / rep(nchar(x), each = prod(dim(ans)[1:2]))
  }
  ans
}

setGeneric("twoWayAlphabetFrequencyByQuality", signature=c("x", "y"),
           function(x, y, as.prob=FALSE, ...)
           standardGeneric("twoWayAlphabetFrequencyByQuality")
           )

setMethod("twoWayAlphabetFrequencyByQuality", "QualityScaledXStringSet",
          function(x, y, as.prob=FALSE, baseOnly=FALSE)
          .XStringSet.two_way_code_frequency_by_quality(x, y,
                                                        quality(x),
                                                        quality(y),
                                                        as.prob,
                                                        baseOnly)
          )
