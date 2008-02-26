### =========================================================================
### BStringCodec objects
### --------------------
###
### A BStringCodec object allows fast mapping between letters and codes.
###
### In addition to the slots 'letters' and 'codes', it has 2 extra slots
### 'enc_lkup' and 'dec_lkup' that are lookup tables. They allow fast
### translation from letters to codes and from codes to letters.
### Those lookup tables are used at the C level for fast encoding/decoding
### of the sequence contained in a DNAString or RNAString object.
### See the buildLookupTable() function for more details about lookup tables.
###
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions.
###

### Could not find a simpler way to get the vector of ascii codes :-/
.letterAsByteVal <- function(letters)
{
    if (!all(nchar(letters) == 1))
        stop("strings in 'letters' must have one single letter")
    as.integer(charToRaw(paste(letters, collapse="")))
}

### Builds a lookup table that can be used for fast mapping from 'keys'
### (unique integers >= 0) to 'vals'.
### The returned value is a vector 'lkup' such that:
###   lkup[keys + 1] is identical to vals
### Note that if 'x' and 'y' are both integer vectors of the same length,
### then lkupx2y <- buildLookupTable(x, y) and lkupy2x <- buildLookupTable(y, x)
### are reverse lookup tables.
### The key property of reverse lookup tables is:
###   lkupy2x[lkupx2y[x + 1]] + 1 is identical to x
### More generally, if 'x', 'y1', 'y2', 'z' verify:
###   a) 'x' is a vector of unique non-negative integers
###   b) 'y1' is a vector of non-negative integers and same length as 'x'
###   c) 'y2' is a vector of unique non-negative integers
###   d) 'z' is a vector of same length as 'y2'
### and if 'lkupx2y' and 'lkupy2z' are the lookup tables from 'x' to 'y1'
### and from 'y2' to 'z' (respectively), then this table:
###   lkupx2z <- lkupy2z[lkupx2y + 1]
### is the lookup table from 'x' to the subset of 'z' defined by
###   lkupx2z[x + 1]
### Note that 'lkupx2z[x + 1]' is exactly the same as 'lkupy2z[y1 + 1]'.

buildLookupTable <- function(keys, vals)
{
    if (!is.integer(keys) || min(keys) < 0)
        stop("'keys' must be a vector of non-negative integers")
    if (any(duplicated(keys)))
        stop("'keys' are not unique")
    if (length(keys) != length(vals))
        stop("'keys' and 'vals' must have the same length")
    table <- vector(mode=typeof(vals), length=max(keys)+1)
    table[] <- NA
    for (i in 1:length(keys))
        table[keys[i] + 1] <- vals[[i]]
    table
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The BStringCodec class.
###

setClass("BStringCodec",
    representation(
        letters="character",
        codes="integer",
        enc_lkup="integer",    # Lookup table for fast encoding
        dec_lkup="integer"     # Lookup table for fast decoding
    )
)

setMethod("initialize", "BStringCodec",
    function(.Object, letters, codes, extra_letters=NULL, extra_codes=NULL)
    {
        letter_byte_vals <- .letterAsByteVal(letters)
        codes <- as.integer(codes)
        .Object@letters <- letters
        .Object@codes <- codes
        .Object@dec_lkup <- buildLookupTable(codes, letter_byte_vals)
        if (!is.null(extra_letters)) {
            letter_byte_vals <- c(letter_byte_vals, .letterAsByteVal(extra_letters))
            codes <- c(codes, as.integer(extra_codes))
        }
        .Object@enc_lkup <- buildLookupTable(letter_byte_vals, codes)
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The DNA and RNA alphabets and codecs.
###

DNA_BASE_CODES <- c(A=1L, C=2L, G=4L, T=8L)
RNA_BASE_CODES <- c(A=1L, C=2L, G=4L, U=8L)

.IUPACcodes <- function(baseCodes)
{
    baseIsU <- names(baseCodes) == "U"
    if (any(baseIsU))
        names(baseCodes)[baseIsU] <- "T"
    code_list <- strsplit(IUPAC_CODE_MAP, "", fixed=TRUE)
    codes <- sapply(code_list, function(x) sum(baseCodes[x]))
    if (any(baseIsU))
        names(codes)[names(codes) == "T"] <- "U"
    codes
}

.DNAorRNAcodes <- function(baseCodes, baseOnly)
{
    if (!isTRUEorFALSE(baseOnly))
        stop("'baseOnly' must be 'TRUE' or 'FALSE'")
    codes <- .IUPACcodes(baseCodes)
    if (baseOnly)
        codes[names(codes) %in% names(baseCodes)]
    else
        c(codes, '-'=16L)
}

DNAcodes <- function(baseOnly) .DNAorRNAcodes(DNA_BASE_CODES, baseOnly)
RNAcodes <- function(baseOnly) .DNAorRNAcodes(RNA_BASE_CODES, baseOnly)

### DNA and RNA alphabets
.DNA_CODES <- DNAcodes(FALSE)
.RNA_CODES <- RNAcodes(FALSE)
DNA_ALPHABET <- names(.DNA_CODES)
RNA_ALPHABET <- names(.RNA_CODES)


### DNA and RNA codecs
.BStringCodec.DNAorRNA <- function(codes)
{
    letters <- names(codes)
    extra_letters <- setdiff(tolower(letters), letters)
    extra_codes <- codes[toupper(extra_letters)]
    new("BStringCodec", letters, codes, extra_letters, extra_codes)
}

DNA_STRING_CODEC <- .BStringCodec.DNAorRNA(.DNA_CODES)
RNA_STRING_CODEC <- .BStringCodec.DNAorRNA(.RNA_CODES)


### Return the lookup table that transforms a DNA sequence into its
### complementary sequence.
getDNAComplementLookup <- function()
{
    complement_base_codes <- c(A=DNA_BASE_CODES["T"][[1]],
                               C=DNA_BASE_CODES["G"][[1]],
                               G=DNA_BASE_CODES["C"][[1]],
                               T=DNA_BASE_CODES["A"][[1]])
    complement_codes <- .DNAorRNAcodes(complement_base_codes, FALSE)
    complement_codec <- .BStringCodec.DNAorRNA(complement_codes)
    complement_codec@enc_lkup[DNA_STRING_CODEC@dec_lkup + 1]
}

### Return the lookup table that transforms an RNA sequence into its
### complementary sequence.
getRNAComplementLookup <- function()
{
    complement_base_codes <- c(A=RNA_BASE_CODES["U"][[1]],
                               C=RNA_BASE_CODES["G"][[1]],
                               G=RNA_BASE_CODES["C"][[1]],
                               U=RNA_BASE_CODES["A"][[1]])
    complement_codes <- .DNAorRNAcodes(complement_base_codes, FALSE)
    complement_codec <- .BStringCodec.DNAorRNA(complement_codes)
    complement_codec@enc_lkup[RNA_STRING_CODEC@dec_lkup + 1]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Add extra codecs below...


