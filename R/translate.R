### =========================================================================
### DNA/RNA transcription/translation & related functions
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some simple functions for different kinds of DNA <-> RNA transformations.
###

transcribe <- function(x)
{
    if (!is(x, "DNAString")) stop("transcribe() only works on DNA input")
    RNAString(complement(x))
}

cDNA <- function(x)
{
    if (!is(x, "RNAString")) stop("cDNA() only works on RNA input")
    DNAString(complement(x))
}

dna2rna <- function(x)
{
    if (!is(x, "DNAString")) stop("dna2rna() only works on DNA input")
    RNAString(x)
}

rna2dna <- function(x)
{
    if (!is(x, "RNAString")) stop("rna2dna() only works on RNA input")
    DNAString(x)
}



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "codons" generic and methods.
###

setGeneric("codons", signature="x",
    function(x) standardGeneric("codons")
)


### 'x' must be a DNAString or RNAString object
.XString.codons <- function(x)
{
    if (length(x) %% 3L != 0L)
        warning("the number of nucleotides in 'x' is not a multiple of 3")
    ans <- successiveViews(x, rep.int(3L, length(x) %/% 3L))
    if (alphabetFrequency(ans, baseOnly=TRUE, collapse=TRUE)[["other"]] != 0)
        stop("some trinucleotides in 'x' contain non-base letters")
    ans
}

### 'x' is assumed to in the coding strand of the DNA
setMethod("codons", "DNAString", function(x) .XString.codons(x))

setMethod("codons", "RNAString", function(x) .XString.codons(x))


### 'x' must be a MaskedDNAString or MaskedRNAString object.
### Return 1 view per codon. Each view is guaranteed to contain exactly
### 3 base letters plus eventually some '+' letters (removed during
### translation).
.MaskedXString.codons <- function(x)
{
    if (nchar(x) %% 3L != 0L)
        warning("the number of unmasked nucleotides in 'x' is not a multiple of 3")
    if (alphabetFrequency(x, baseOnly=TRUE)[["other"]] != 0)
        stop("some trinucleotides in 'x' contain non-base letters")
    ans_length <- nchar(x) %/% 3L
    ans_start <- integer(ans_length)
    ans_end <- integer(ans_length)
    x0 <- injectHardMask(x)
    codon_start <- 1L
    for (i in seq_len(ans_length)) {
        while (letter(x0, codon_start) == '+')
            codon_start <- codon_start + 1L
        codon_end <- codon_start
        codon_nchar <- 1L
        while (codon_nchar < 3L) {
            codon_end <- codon_end + 1L
            if (letter(x0, codon_end) != '+')
                codon_nchar <- codon_nchar + 1L
        }
        ans_start[i] <- codon_start
        ans_end[i] <- codon_end
        codon_start <- codon_end + 1L
    }
    Views(x0, start=ans_start, end=ans_end)
}

setMethod("codons", "MaskedDNAString", function(x) .MaskedXString.codons(x))

setMethod("codons", "MaskedRNAString", function(x) .MaskedXString.codons(x))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "translate" generic and methods.
###

setGeneric("translate", function(x) standardGeneric("translate"))

### 'x' must be an XStringViews object as returned by .XString.codons()
### or .MaskedXString.codons()
.translate.codons.as.character <- function(x)
{
    x3 <- gsub("+", "", as.character(x), fixed=TRUE)
    paste(GENETIC_CODE[x3], collapse="")
}

setMethod("translate", "DNAString",
    function(x) AAString(.translate.codons.as.character(.XString.codons(x)))
)

setMethod("translate", "RNAString",
    function(x) translate(DNAString(x))
)

setMethod("translate", "DNAStringSet",
    function(x)
    {
        if (any(width(x) %% 3L != 0L))
            warning("the number of nucleotides in some element of 'x' is not a multiple of 3")
        AAStringSet(
            sapply(
                seq_len(x),
                function(i) .translate.codons.as.character(.XString.codons(x[[i]]))
            )
        )
    }
)

setMethod("translate", "RNAStringSet", function(x) translate(DNAStringSet(x)))

setMethod("translate", "MaskedDNAString",
    function(x) AAString(.translate.codons.as.character(.MaskedXString.codons(x)))
)

setMethod("translate", "MaskedRNAString",
    function(x)
    {
        ## FIXME: Workaround until as(x, "MaskedDNAString") is available
        y <- new("MaskedDNAString", unmasked=DNAString(unmasked(x)), masks=masks(x))
        translate(y)
    }
)

