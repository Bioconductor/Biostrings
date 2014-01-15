### =========================================================================
### DNA/RNA transcription/translation & related functions
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some simple functions for different kinds of DNA <-> RNA transformations.
###

transcribe <- function(x)
{
    msg <- c("  transcribe() is deprecated. ",
             "Please use 'RNAString(complement(x))' instead\n",
             "  (which is how 'transcribe(x)' is implemented).")
    .Deprecated(msg=paste0(msg, collapse=""))
    if (!is(x, "DNAString")) stop("transcribe() only works on DNA input")
    RNAString(complement(x))
}

cDNA <- function(x)
{
    msg <- c("  cDNA() is deprecated. ",
             "Please use 'DNAString(complement(x))' instead\n",
             "  (which is how 'cDNA(x)' is implemented).")
    .Deprecated(msg=paste0(msg, collapse=""))
    if (!is(x, "RNAString")) stop("cDNA() only works on RNA input")
    DNAString(complement(x))
}

dna2rna <- function(x)
{
    msg <- "  dna2rna() is deprecated. Please use RNAString() instead."
    .Deprecated(msg=msg)
    if (!is(x, "DNAString")) stop("dna2rna() only works on DNA input")
    RNAString(x)
}

rna2dna <- function(x)
{
    msg <- "  rna2dna() is deprecated. Please use DNAString() instead."
    .Deprecated(msg=msg)
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

### Returns a character vector of length 2.
.normarg_if.fuzzy.codon <- function(if.fuzzy.codon)
{
    if (!is.character(if.fuzzy.codon)
     || !(length(if.fuzzy.codon) == 1L || length(if.fuzzy.codon) == 2L))
        stop("'if.fuzzy.codon' must be a single string or ",
             "a character vector of length 2")
    if (length(if.fuzzy.codon) == 1L) {
        choices <- c("error", "solve", "error.if.X", "X")
        if.fuzzy.codon <- match.arg(if.fuzzy.codon, choices)
        if.fuzzy.codon <- switch(if.fuzzy.codon,
                                 error=c("error", "error"),
                                 solve=c("solve", "X"),
                                 error.if.X=c("solve", "error"),
                                 X=c("X", "X"))
        return(if.fuzzy.codon)
    }
    ## From now on, 'if.fuzzy.codon' is guaranteed to have length 2.
    if.non.ambig <- if.fuzzy.codon[[1L]]
    choices1 <- c("error", "solve", "X")
    if.non.ambig <- match.arg(if.non.ambig, choices1)

    if.ambig <- if.fuzzy.codon[[2L]]
    choices2 <- c("error", "X")
    if.ambig <- match.arg(if.ambig, choices2)
    c(if.non.ambig, if.ambig)
}

.makeFuzzyGeneticCode <- function(keep.ambig.codons=FALSE)
{
    if (!isTRUEorFALSE(keep.ambig.codons))
        stop("'keep.ambig.codons' must be TRUE or FALSE")
    iupac_codes <- names(IUPAC_CODE_MAP)
    fuzzy_codons <- mkAllStrings(iupac_codes, 3L)
    nonfuzzy_codons <- DNAStringSet(names(GENETIC_CODE))
    fuzzy2nonfuzzy <- vwhichPDict(nonfuzzy_codons, DNAStringSet(fuzzy_codons),
                                  fixed="pattern")
    fuzzy2AAs <- relist(unname(GENETIC_CODE)[unlist(fuzzy2nonfuzzy)],
                        PartitioningByEnd(fuzzy2nonfuzzy))
    fuzzy2AAs <- unique(fuzzy2AAs)
    nAAs <- elementLengths(fuzzy2AAs)
    names(fuzzy2AAs) <- fuzzy_codons
    stopifnot(all(nAAs >= 1L))
    if (keep.ambig.codons) {
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

setMethod("translate", "DNAStringSet",
    function(x, if.fuzzy.codon="error")
    {
        if.fuzzy.codon <- .normarg_if.fuzzy.codon(if.fuzzy.codon)
        if.non.ambig <- if.fuzzy.codon[[1L]]
        if.ambig <- if.fuzzy.codon[[2L]]
        if (if.non.ambig == "error" && if.ambig == "error") {
            codon_alphabet <- DNA_BASES
            genetic_code <- GENETIC_CODE
        } else {
            codon_alphabet <- names(IUPAC_CODE_MAP)
            genetic_code <- .makeFuzzyGeneticCode(keep.ambig.codons=TRUE)
        }
        dna_codes <- DNAcodes(baseOnly=FALSE)
        skip_code <- dna_codes[["+"]]
        lkup <- .makeTranslationLkup(codon_alphabet, genetic_code)
        .Call2("DNAStringSet_translate",
              x, skip_code, dna_codes[codon_alphabet], lkup,
              if.non.ambig, if.ambig,
              PACKAGE="Biostrings")
    }
)

setMethod("translate", "RNAStringSet",
    function(x, if.fuzzy.codon="error")
        translate(DNAStringSet(x), if.fuzzy.codon=if.fuzzy.codon)
)

setMethod("translate", "DNAString",
    function(x, if.fuzzy.codon="error")
        translate(DNAStringSet(x), if.fuzzy.codon=if.fuzzy.codon)[[1L]]
)

setMethod("translate", "RNAString",
    function(x, if.fuzzy.codon="error")
        translate(DNAString(x), if.fuzzy.codon=if.fuzzy.codon)
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

