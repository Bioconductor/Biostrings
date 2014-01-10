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

setGeneric("translate", signature="x",
    function(x, if.fuzzy.codon="error") standardGeneric("translate")
)

.makeFuzzyGeneticCode <- function(keep.ambiguous.codons=FALSE)
{
    if (!isTRUEorFALSE(keep.ambiguous.codons))
        stop("'keep.ambiguous.codons' must be TRUE or FALSE")
    iupac_codes <- names(IUPAC_CODE_MAP)
    fuzzy_codons <- mkAllStrings(iupac_codes, 3L)
    fixed_codons <- DNAStringSet(names(GENETIC_CODE))
    pdict <- PDict(fixed_codons)
    fuzzy2fixed <- vwhichPDict(pdict, DNAStringSet(fuzzy_codons),
                               fixed="pattern")
    fuzzy2AAs <- relist(unname(GENETIC_CODE)[unlist(fuzzy2fixed)],
                              PartitioningByEnd(fuzzy2fixed))
    fuzzy2AAs <- unique(fuzzy2AAs)
    nAAs <- elementLengths(fuzzy2AAs)
    names(fuzzy2AAs) <- fuzzy_codons
    stopifnot(all(nAAs >= 1L))
    if (keep.ambiguous.codons) {
        idx <- which(nAAs != 1L)
        fuzzy2AAs[idx] <- "X"
    } else {
        idx <- which(nAAs == 1L)
        fuzzy2AAs <- fuzzy2AAs[idx]
    }
    unlist(fuzzy2AAs)
}

.makeTranslationLkup <- function(codon_alphabet, genetic_code)
{
    codons <- mkAllStrings(codon_alphabet, 3)
    i <- match(codons, names(genetic_code))
    if (any(is.na(i)))
        stop("some codons are not in 'genetic_code'")
    as.integer(charToRaw(paste0(genetic_code[i], collapse="")))
}

setMethod("translate", "DNAString",
    function(x, if.fuzzy.codon="error")
        translate(DNAStringSet(x), if.fuzzy.codon=if.fuzzy.codon)[[1L]]
)

setMethod("translate", "RNAString",
    function(x, if.fuzzy.codon="error")
        translate(DNAString(x), if.fuzzy.codon=if.fuzzy.codon)
)

setMethod("translate", "DNAStringSet",
    function(x, if.fuzzy.codon="error")
    {
        if.fuzzy.codon <- match.arg(if.fuzzy.codon, c("error", "translate"))
        if (if.fuzzy.codon == "error") {
            codon_alphabet <- DNA_BASES
            genetic_code <- GENETIC_CODE
        } else {
            codon_alphabet <- names(IUPAC_CODE_MAP)
            genetic_code <- .makeFuzzyGeneticCode(keep.ambiguous.codons=TRUE)
        }
        dna_codes <- DNAcodes(baseOnly=FALSE)
        skipcode <- dna_codes[["+"]]
        lkup <- .makeTranslationLkup(codon_alphabet, genetic_code)
        .Call2("DNAStringSet_translate",
              x, skipcode, dna_codes[codon_alphabet], lkup,
              PACKAGE="Biostrings")
    }
)

setMethod("translate", "RNAStringSet",
    function(x, if.fuzzy.codon="error")
        translate(DNAStringSet(x), if.fuzzy.codon=if.fuzzy.codon)
)

setMethod("translate", "MaskedDNAString",
    function(x, if.fuzzy.codon="error")
        translate(injectHardMask(x), if.fuzzy.codon=if.fuzzy.codon)
)

setMethod("translate", "MaskedRNAString",
    function(x, if.fuzzy.codon="error")
    {
        ## FIXME: Workaround until as(x, "MaskedDNAString") is available
        y <- new("MaskedDNAString", unmasked=DNAString(unmasked(x)), masks=masks(x))
        translate(y, if.fuzzy.codon=if.fuzzy.codon)
    }
)

