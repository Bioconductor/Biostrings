### =========================================================================
### Translating DNA/RNA sequences
### -------------------------------------------------------------------------


.normarg_genetic.code <- function(genetic.code)
{
    if (!is.character(genetic.code) || any(is.na(genetic.code)))
        stop("'genetic.code' must be a character vector with no NAs")
    if (!identical(names(genetic.code), names(GENETIC_CODE)))
        stop("'genetic.code' must have the same names as ",
             "predefined constant GENETIC_CODE")
    if (!all(nchar(genetic.code) == 1L))
        stop("'genetic.code' must contain 1-letter strings")
    ## Just a warning for now. Might become an error in the future.
    if (!all(genetic.code %in% AA_ALPHABET))
        warning("some codons in 'genetic.code' are mapped to letters ",
                "not in the Amino Acid\n  alphabet (AA_ALPHABET)")
    alt_init_codons <- attr(genetic.code, "alt_init_codons", exact=TRUE)
    if (is.null(alt_init_codons)
     || !is.character(alt_init_codons)
     || any(is.na(alt_init_codons))
     || anyDuplicated(alt_init_codons)
     || !all(alt_init_codons %in% names(genetic.code)))
        stop(wmsg("'genetic.code' must have an \"alt_init_codons\" attribute ",
                  "that lists alternative initiation codons"))
    genetic.code
}

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

.make_fuzzy_genetic_code <- function(genetic.code, keep.ambig.codons=FALSE)
{
    if (!isTRUEorFALSE(keep.ambig.codons))
        stop("'keep.ambig.codons' must be TRUE or FALSE")
    iupac_codes <- names(IUPAC_CODE_MAP)
    fuzzy_codons <- mkAllStrings(iupac_codes, 3L)
    nonfuzzy_codons <- DNAStringSet(names(genetic.code))
    fuzzy2nonfuzzy <- vwhichPDict(nonfuzzy_codons, DNAStringSet(fuzzy_codons),
                                  fixed="pattern")
    fuzzy2AAs <- relist(unname(genetic.code)[unlist(fuzzy2nonfuzzy)],
                        PartitioningByEnd(fuzzy2nonfuzzy))
    fuzzy2AAs <- unique(fuzzy2AAs)
    nAAs <- elementNROWS(fuzzy2AAs)
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

.aa2byte <- function(aa) as.integer(charToRaw(paste0(aa, collapse="")))

.make_translation_lkup <- function(codon_alphabet, genetic.code)
{
    codons <- mkAllStrings(codon_alphabet, 3)
    m <- match(codons, names(genetic.code))
    if (any(is.na(m)))
        stop("some codons are not in 'genetic.code'")
    .aa2byte(genetic.code[m])
}

.translate <- function(x, genetic.code=GENETIC_CODE, if.fuzzy.codon="error")
{
    init_genetic_code <- genetic_code <- .normarg_genetic.code(genetic.code)
    alt_init_codons <- attr(genetic_code, "alt_init_codons")
    init_genetic_code[alt_init_codons] <- "M"

    if.fuzzy.codon <- .normarg_if.fuzzy.codon(if.fuzzy.codon)
    if.non.ambig <- if.fuzzy.codon[[1L]]
    if.ambig <- if.fuzzy.codon[[2L]]

    if (if.non.ambig == "error" && if.ambig == "error") {
        codon_alphabet <- DNA_BASES
    } else {
        codon_alphabet <- names(IUPAC_CODE_MAP)
        genetic_code <- .make_fuzzy_genetic_code(genetic_code,
                                                 keep.ambig.codons=TRUE)
        init_genetic_code <- .make_fuzzy_genetic_code(init_genetic_code,
                                                      keep.ambig.codons=TRUE)
    }

    lkup <- .make_translation_lkup(codon_alphabet, genetic_code)
    init_lkup <- .make_translation_lkup(codon_alphabet, init_genetic_code)
    dna_codes <- DNAcodes(baseOnly=FALSE)
    skip_code <- dna_codes[["+"]]
    ans <- .Call2("DNAStringSet_translate",
                  x, skip_code, dna_codes[codon_alphabet],
                  lkup, init_lkup,
                  if.non.ambig, if.ambig,
                  PACKAGE="Biostrings")
    names(ans) <- names(x)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "translate" generic and methods.
###

setGeneric("translate", signature="x",
    function(x, genetic.code=GENETIC_CODE, if.fuzzy.codon="error")
        standardGeneric("translate")
)

setMethod("translate", "DNAStringSet", .translate)

setMethod("translate", "RNAStringSet", .translate)

setMethod("translate", "DNAString",
    function(x, genetic.code=GENETIC_CODE, if.fuzzy.codon="error")
        translate(DNAStringSet(x), genetic.code=genetic.code,
                  if.fuzzy.codon=if.fuzzy.codon)[[1L]]
)

setMethod("translate", "RNAString",
    function(x, genetic.code=GENETIC_CODE, if.fuzzy.codon="error")
        translate(RNAStringSet(x), genetic.code=genetic.code,
                  if.fuzzy.codon=if.fuzzy.codon)[[1L]]
)

setMethod("translate", "MaskedDNAString",
    function(x, genetic.code=GENETIC_CODE, if.fuzzy.codon="error")
        translate(injectHardMask(x), genetic.code=genetic.code,
                  if.fuzzy.codon=if.fuzzy.codon)
)

setMethod("translate", "MaskedRNAString",
    function(x, genetic.code=GENETIC_CODE, if.fuzzy.codon="error")
        translate(injectHardMask(x), genetic.code=genetic.code,
                  if.fuzzy.codon=if.fuzzy.codon)
)


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

