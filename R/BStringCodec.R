# ===========================================================================
# The BStringCodec class
# ---------------------------------------------------------------------------

# A BStringCodec object is used to store a mapping between letters
# and their codes.
# The 'initialize' method checks that slots 'letters' and 'codes'
# (both "bbuf" objects) have the same length.

setClass(
    "BStringCodec",
    representation(
        letters="cbuf",
        codes="bbuf",
        enc_hash="bbuf",    # Hash table for encoding
        dec_hash="bbuf"     # Hash table for decoding
    )
)

# The first value in a hash table is reserved: it must contain
# the "hole" value.
# The key property for the encoding and decoding hash tables is:
#   dec_hash[enc_hash[x + 2] + 2] == x
# for every x to be encoded.

buildHash <- function(x, y, hole)
{
    if (!is.integer(x) || !is.integer(y))
        stop("'x' and 'y' must be integer vectors") 
    if (length(x) != length(y))
        stop("lengths of 'x' and 'y' are different")
    hash <- bbuf(max(x) + 2)
    hash[] <- hole
    for (i in 1:length(x)) {
        if (y[i] == hole)
            stop("'hole' must be different from any value in 'y'")
        hash[x[i] + 2] <- y[i]
    }
    hash
}

# 'enc_hole' is used to mark holes in the encoding hash table. Its value must
# be different from any code.
# 'dec_hole' is used to mark holes in the decoding hash table. Its value must
# be different from any letter.
setMethod("initialize", "BStringCodec",
    function(.Object, letters, codes, enc_hole, dec_hole)
    {
        if (!is.character(letters) || length(letters) != 1)
            stop("'letters' must be a single string")
        if (nchar(letters) == 0)
            stop("no letter provided")
        if (nchar(letters) != length(codes))
            stop("number of letters and codes are different")
        if (!is.integer(codes))
            stop("'storage.mode(codes)' must be \"integer\"")
        .Object@letters <- new("cbuf", nchar(letters))
        .Object@letters[] <- letters
        .Object@codes <- bbuf(length(codes))
        .Object@codes[] <- codes
        ord <- as.integer(.Object@letters)
        .Object@enc_hash <- buildHash(ord, codes, enc_hole)
        .Object@dec_hash <- buildHash(codes, ord, dec_hole)
        .Object
    }
)

DNAcodec <- function()
{
        letters <- c("-TGCANBDHKMRSVWY")
        codes <- c(1, 2, 4, 8, 16, 30, 14, 22, 26, 6, 24, 20, 12, 28, 18, 10)
        new("BStringCodec", letters, as.integer(codes), 99, 199)
}

# The only difference with DNAcodec() is that letter "T" is replaced
# by letter "U". That's all folks!
RNAcodec <- function()
{
        letters <- c("-UGCANBDHKMRSVWY")
        codes <- c(1, 2, 4, 8, 16, 30, 14, 22, 26, 6, 24, 20, 12, 28, 18, 10)
        new("BStringCodec", letters, as.integer(codes), 99, 199)
}

# More codecs to come below...


