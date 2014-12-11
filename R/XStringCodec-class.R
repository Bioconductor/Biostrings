### =========================================================================
### XStringCodec objects
### --------------------
###
### A XStringCodec object allows fast mapping between letters and codes.
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
### More generally, if 'x', 'y1', 'y2', and 'z' verify:
###   a) 'x' is a vector of unique non-negative integers
###   b) 'y1' is a vector of non-negative integers with same length as 'x'
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
    ## Checking 'keys'.
    if (!is.integer(keys))
        stop("'keys' must be a an integer vector")
    if (any(is.na(keys)))
        stop("'keys' cannot contain NAs")
    keys_len <- length(keys)
    if (keys_len != 0L && min(keys) < 0L)
        stop("'keys' cannot contain negative integers")
    if (any(duplicated(keys)))
        stop("'keys' cannot contain duplicates")

    ## Checking 'vals'.
    if (!is.atomic(vals) || length(vals) != keys_len)
        stop("'vals' must be a vector of the length of 'keys'")

    ## Build the lookup table ('ans') and return it.
    if (keys_len == 0L) {
        ans_len <- 0L  # but could be anything as long as we fill with NAs
    } else {
        ans_len <- max(keys) + 1L
    }
    ans <- vector(mode=typeof(vals), length=ans_len)
    ans[] <- NA
    ans[keys + 1L] <- vals
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The XStringCodec class.
###

setClass("XStringCodec",
    representation(
        letters="character",
        codes="integer",
        enc_lkup="integer",    # Lookup table for fast encoding
        dec_lkup="integer"     # Lookup table for fast decoding
    )
)

setMethod("initialize", "XStringCodec",
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
        stop("'baseOnly' must be TRUE or FALSE")
    codes <- .IUPACcodes(baseCodes)
    if (baseOnly) {
        codes[names(codes) %in% names(baseCodes)]
    } else {
        c(codes, `-`=16L, `+`=32L, `.`=64L)
    }
}

DNAcodes <- function(baseOnly) .DNAorRNAcodes(DNA_BASE_CODES, baseOnly)
RNAcodes <- function(baseOnly) .DNAorRNAcodes(RNA_BASE_CODES, baseOnly)

### DNA and RNA alphabets.
DNA_CODES <- DNAcodes(FALSE)
RNA_CODES <- RNAcodes(FALSE)
DNA_ALPHABET <- names(DNA_CODES)
RNA_ALPHABET <- names(RNA_CODES)

### DNA_BASES could be defined more simply as being 'names(DNA_BASE_CODES)'
### but calling DNAcodes() gives us the guarantee that the order of the
### bases will be consistent with DNA_ALPHABET.
DNA_BASES <- names(DNAcodes(TRUE))
RNA_BASES <- names(RNAcodes(TRUE))


### DNA and RNA codecs.
.XStringCodec.DNAorRNA <- function(codes)
{
    letters <- names(codes)
    extra_letters <- setdiff(tolower(letters), letters)
    extra_codes <- codes[toupper(extra_letters)]
    new("XStringCodec", letters, codes, extra_letters, extra_codes)
}

DNA_STRING_CODEC <- .XStringCodec.DNAorRNA(DNA_CODES)
RNA_STRING_CODEC <- .XStringCodec.DNAorRNA(RNA_CODES)


### Return the lookup table that transforms a DNA sequence into its
### complementary sequence.
getDNAComplementLookup <- function()
{
    complement_base_codes <- c(A=DNA_BASE_CODES["T"][[1]],
                               C=DNA_BASE_CODES["G"][[1]],
                               G=DNA_BASE_CODES["C"][[1]],
                               T=DNA_BASE_CODES["A"][[1]])
    complement_codes <- .DNAorRNAcodes(complement_base_codes, FALSE)
    complement_codec <- .XStringCodec.DNAorRNA(complement_codes)
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
    complement_codec <- .XStringCodec.DNAorRNA(complement_codes)
    complement_codec@enc_lkup[RNA_STRING_CODEC@dec_lkup + 1]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Add extra codecs below...


