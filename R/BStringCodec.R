### Could not find a simpler way to get the vector of ascii codes :-/
letterAsByteVal <- function(letters)
{
    if (!all(nchar(letters) == 1))
        stop("strings in 'letters' must have one single letter")
    as.integer(charToRaw(paste(letters, collapse="")))
}

### Builds a lookup table for 'keys' (unique integers >= 0) and 'vals'.
### The returned value is a vector 'lkup' such that:
###   lkup[keys + 1] is identical to vals
### Note that if 'x' and 'y' are both integer vectors of the same length,
### then lkupxy <- buildLookupTable(x, y) and lkupyx <- buildLookupTable(y, x)
### are reverse lookup tables.
### The key property of reverse lookup tables is:
###   lkupyx[lkupxy[x + 1]] + 1 is identical to x

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


# ===========================================================================
# The BStringCodec class
# ---------------------------------------------------------------------------

# A BStringCodec object is used to store a mapping between letters
# and their codes.
# The 'initialize' method checks that slots 'letters' and 'codes'
# (both "CharBuffer" objects) have the same length.

setClass(
    "BStringCodec",
    representation(
        letters="character",
        codes="integer",
        enc_lkup="integer",    # Lookup table for fast encoding
        dec_lkup="integer"     # Lookup table for fast decoding
    )
)

setMethod("initialize", "BStringCodec",
    function(.Object, letters, codes, extra_letters, extra_codes)
    {
        letter_byte_vals <- letterAsByteVal(letters)
        codes <- as.integer(codes)
        .Object@letters <- letters
        .Object@codes <- codes
        .Object@dec_lkup <- buildLookupTable(codes, letter_byte_vals)
        if (!missing(extra_letters)) {
            letter_byte_vals <- c(letter_byte_vals, letterAsByteVal(extra_letters))
            extra_codes <- as.integer(extra_codes)
            codes <- c(codes, extra_codes)
        }
        .Object@enc_lkup <- buildLookupTable(letter_byte_vals, codes)
        .Object
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setGeneric("alphabet", function(x) standardGeneric("alphabet"))

setMethod("alphabet", "BStringCodec", function(x) x@letters)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IUPAC extended genetic alphabet
#   R :== [GA]
#   Y :== [TC]
#   M :== [AC]
#   K :== [GT]
#   S :== [GC]
#   W :== [AT]
#   H :== [ACT]
#   B :== [GTC]
#   V :== [GCA]
#   D :== [GAT]
#   N :== [GATC]

BStringCodec.DNA <- function()
{
    letters <- strsplit("ACGTMRSVWYHKDBN-", "", fixed=TRUE)[[1]]
    codes <- c(1, 2, 4, 8, 3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16)
    extra_letters <- strsplit("acgtmrsvwyhkdbn", "", fixed=TRUE)[[1]]
    extra_codes <- c(1, 2, 4, 8, 3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15)
    new("BStringCodec", letters, codes, extra_letters, extra_codes)
}

BStringCodec.RNA <- function()
{
    letters <- strsplit("UGCAKYSBWRDMHVN-", "", fixed=TRUE)[[1]]
    codes <- c(1, 2, 4, 8, 3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16)
    extra_letters <- strsplit("ugcakysbwrdmhvn", "", fixed=TRUE)[[1]]
    extra_codes <- c(1, 2, 4, 8, 3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15)
    new("BStringCodec", letters, codes, extra_letters, extra_codes)
}

# More codecs to come below...


