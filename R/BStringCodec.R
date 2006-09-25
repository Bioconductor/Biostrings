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
        letters="PrintableCharBuffer",
        codes="CharBuffer",
        enc_hash="CharBuffer",    # Hash table for encoding
        dec_hash="CharBuffer"     # Hash table for decoding
    )
)

# The first value in a hash table is reserved: it must contain
# the "hole" value.
# The key property for the encoding and decoding hash tables is:
#   enc_hash[dec_hash[x + 2] + 2] == x
# for every x to be decoded.

buildCodecHashTable <- function(x, y, hole)
{
    if (!is.integer(x) || !is.integer(y))
        stop("'x' and 'y' must be integer vectors") 
    if (length(x) != length(y))
        stop("lengths of 'x' and 'y' are different")
    hash <- CharBuffer(max(x) + 2)
    hash[] <- hole
    for (i in 1:length(x)) {
        if (y[i] == hole)
            stop("'hole' must be different from any value in 'y'")
        hash[x[i] + 2] <- y[i]
    }
    hash
}

checkCodecLettersAndCodes <- function(letters, codes)
{
    if (!is.character(letters) || length(letters) != 1)
        stop("'letters' must be a single string")
    if (nchar(letters) == 0)
        stop("no letter provided")
    if (nchar(letters) != length(codes))
        stop("number of letters and codes are different")
    if (!is.integer(codes))
        stop("'storage.mode(codes)' must be \"integer\"")
}

# 'enc_hole' is used to mark holes in the encoding hash table. Its value must
# be different from any code.
# 'dec_hole' is used to mark holes in the decoding hash table. Its value must
# be different from any letter.
setMethod("initialize", "BStringCodec",
    function(.Object, letters, codes, enc_hole, dec_hole, extra_letters, extra_codes)
    {
        checkCodecLettersAndCodes(letters, codes)
        .Object@letters <- new("PrintableCharBuffer", nchar(letters))
        .Object@letters[] <- letters
        .Object@codes <- CharBuffer(length(codes))
        .Object@codes[] <- codes
        ord <- as.integer(.Object@letters)
        .Object@dec_hash <- buildCodecHashTable(codes, ord, dec_hole)
        if (!missing(extra_letters)) {
            checkCodecLettersAndCodes(extra_letters, extra_codes)
            tmp <- new("PrintableCharBuffer", nchar(extra_letters))
            tmp[] <- extra_letters
            ord <- c(ord, as.integer(tmp))
            codes <- c(codes, extra_codes)
        }
        .Object@enc_hash <- buildCodecHashTable(ord, codes, enc_hole)
        .Object
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setGeneric("alphabet", function(x) standardGeneric("alphabet"))

setMethod("alphabet", "BStringCodec", function(x) x@letters[])


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
    letters <- c("ACGTMRSVWYHKDBN-")
    codes <- c(1, 2, 4, 8, 3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16)
    extra_letters <- c("acgtmrsvwyhkdbn")
    extra_codes <- c(1, 2, 4, 8, 3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15)
    new("BStringCodec", letters, as.integer(codes), 99, 199,
        extra_letters, as.integer(extra_codes))
}

BStringCodec.RNA <- function()
{
    letters <- c("UGCAKYSBWRDMHVN-")
    codes <- c(1, 2, 4, 8, 3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16)
    extra_letters <- c("ugcakysbwrdmhvn")
    extra_codes <- c(1, 2, 4, 8, 3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15)
    new("BStringCodec", letters, as.integer(codes), 99, 199,
        extra_letters, as.integer(extra_codes))
}

# More codecs to come below...


