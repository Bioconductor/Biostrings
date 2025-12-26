### =========================================================================
### The seqtype() generic & related functions
### -------------------------------------------------------------------------
###
### Most sequence containers in Biostrings have a "sequence type" that
### reflects the type of sequences that it represents as well as the
### encoding that is used to store the sequence internally:
###
###   sequence  |                           |              |
###   type      | description               | alphabet     | encoded
###   ----------|---------------------------|--------------|--------
###   "B"       | general purpose string(s) | bytes 0-255  | no
###   "DNA"     | DNA sequence(s)           | DNA_ALPHABET | yes
###   "RNA"     | RNA sequence(s)           | RNA_ALPHABET | yes
###   "AA"      | amino acid sequence(s)    | AA_ALPHABET  | yes
###
### The seqtype() function returns the sequence type. For
### example 'seqtype(AAString())' returns "AA".
###
### The Modstrings package by Felix Ernst introduces two additional sequence
### types, "ModDNA" and "ModRNA", that are treated as sequence type "B" by
### compatible_seqtypes(), get_seqtype_conversion_lookup(), and
### get_seqtype_switches_before_binary_op() below.
###
### The Structstrings package by Felix Ernst introduces one additional sequence
### types, "DotBracket", that are treated as sequence type "B" by
### compatible_seqtypes(), get_seqtype_conversion_lookup(), and
### get_seqtype_switches_before_binary_op() below.
###
### Unless specified otherwise, things in this file are not exported.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "seqtype" and "seqtype<-" generics.
###
### seqtype() and `seqtype<-`() have methods defined for the 4 basic
### string containers: XString (single sequence), XStringSet (multiple
### sequences), XStringViews (multiple sequences) and MaskedXString (single
### sequence).
###

### Exported.
setGeneric("seqtype", function(x) standardGeneric("seqtype"))

### Exported.
### Use to switch the seqtype of an object (e.g. from "DNA" to "B").
setGeneric("seqtype<-", signature="x",
    function(x, value) standardGeneric("seqtype<-")
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions for which the returned value depends on 'seqtype(x)',
### not on what particular data are in 'x'. Not exported.
###

xsbaseclass <- function(x) paste(seqtype(x), "String", sep="")

setGeneric("xscodes", signature="x",
    function(x, baseOnly=FALSE, ...) standardGeneric("xscodes")
)

setMethod("xscodes", "ANY",
    function(x, baseOnly=FALSE)
    {
        if (!isTRUEorFALSE(baseOnly))
            stop("'baseOnly' must be TRUE or FALSE")
        seqtype <- seqtype(x)
        switch(seqtype,
               DNA=DNAcodes(baseOnly),
               RNA=RNAcodes(baseOnly),
               AA=AAcodes(baseOnly),
               0:255
        )
    }
)

xscodec <- function(x)
{
    switch(seqtype(x),
        DNA=DNA_STRING_CODEC,
        RNA=RNA_STRING_CODEC,
        AA=AA_STRING_CODEC,
        NULL
    )
}

xs_enc_lkup <- function(x)
{
    codec <- xscodec(x)
    if (is.null(codec)) NULL else codec@enc_lkup
}

xs_dec_lkup <- function(x)
{
    codec <- xscodec(x)
    if (is.null(codec)) NULL else codec@dec_lkup
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some restrictions apply for converting from a sequence type to another
### or for comparing XString objects of different sequence types. This is
### due to the fact that XString objects with different sequence types can
### use different encodings for their sequence data (or no encoding at all)
### or simply to the fact that the conversion or comparison doesn't make
### sense from a biological perspective.
### The helper functions below are used internally (they are NOT exported) to
### determine those restrictions.
###

### Note that sequence types "ModDNA", "ModRNA" and "DotBracket" are treated as
### sequence type "B" by compatible_seqtypes(), get_seqtype_conversion_lookup(),
### and get_seqtype_switches_before_binary_op() below.
.SUPPORTED_SEQTYPES <- c("B", "DNA", "RNA", "AA", "ModDNA", "ModRNA",
                         "DotBracket")

compatible_seqtypes <- function(seqtype1, seqtype2)
{
    stopifnot(isSingleString(seqtype1), seqtype1 %in% .SUPPORTED_SEQTYPES,
              isSingleString(seqtype2), seqtype2 %in% .SUPPORTED_SEQTYPES)
    if (seqtype1 == seqtype2 ||
        seqtype1 %in% c("B", "ModDNA", "ModRNA", "DotBracket") ||
        seqtype2 %in% c("B", "ModDNA", "ModRNA", "DotBracket"))
        return(TRUE)
    is_nucleo1 <- seqtype1 %in% c("DNA", "RNA")
    is_nucleo2 <- seqtype2 %in% c("DNA", "RNA")
    is_nucleo1 == is_nucleo2
}

### Exported.
get_seqtype_conversion_lookup <- function(from_seqtype, to_seqtype)
{
    if (!compatible_seqtypes(from_seqtype, to_seqtype))
        stop("incompatible sequence types \"",
             from_seqtype, "\" and \"", to_seqtype, "\"")
    if (from_seqtype %in% c("ModDNA", "ModRNA", "DotBracket"))
        from_seqtype <- "B"
    if (to_seqtype %in% c("ModDNA", "ModRNA", "DotBracket"))
        to_seqtype <- "B"
    if (from_seqtype == to_seqtype)
        return(NULL)
    from_is_nucleo <- from_seqtype %in% c("DNA", "RNA")
    to_is_nucleo <- to_seqtype %in% c("DNA", "RNA")
    if (from_is_nucleo && to_is_nucleo)
        return(NULL)
    if (to_seqtype == "DNA")
        return(DNA_STRING_CODEC@enc_lkup)
    if (to_seqtype == "RNA")
        return(RNA_STRING_CODEC@enc_lkup)
    if (from_seqtype == "DNA")
        return(DNA_STRING_CODEC@dec_lkup)
    if (from_seqtype == "RNA")
        return(RNA_STRING_CODEC@dec_lkup)
    if (to_seqtype == "AA")
        return(AA_STRING_CODEC@enc_lkup)
    if (from_seqtype == "AA")
        return(AA_STRING_CODEC@dec_lkup)
    stop("Biostrings internal error, please report")  # should never happen
}

### Returns a character vector of length 2 indicating the 2 target seqtypes
### two use for the seqtype switches that need to happen before proceeding
### with the binary operation.
.get_seqtype_switches_before_binary_op <- function(seqtype1, seqtype2,
                                                   what_op, class1, class2)
{
    stopifnot(isSingleString(seqtype1), isSingleString(seqtype2))
    if (seqtype1 %in% c("ModDNA", "ModRNA", "DotBracket"))
        seqtype1 <- "B"
    if (seqtype2 %in% c("ModDNA", "ModRNA", "DotBracket"))
        seqtype2 <- "B"
    if (seqtype1 == seqtype2)
        return(c(seqtype1, seqtype2))
    if (!compatible_seqtypes(seqtype1, seqtype2))
        stop(wmsg(what_op, " between a ", class1, " object ",
                  "and a ", class2, " object is not supported"))
    if (seqtype1 != "AA" && seqtype2 == "B" ||
        seqtype2 != "AA" && seqtype1 == "B")
        return(c("B", "B"))
    c(seqtype1, seqtype2)
}

get_seqtype_switches_before_binary_op <- function(x, y, what_op="comparison")
{
    seqtype1 <- try(seqtype(x), silent=TRUE)
    if (is(seqtype1, "try-error"))
        seqtype1 <- "B"
    seqtype2 <- try(seqtype(y), silent=TRUE)
    if (is(seqtype2, "try-error"))
        seqtype2 <- "B"
    .get_seqtype_switches_before_binary_op(seqtype1, seqtype2,
                                           what_op, class(x), class(y))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### alphabet()
###
### Exported.

### Could be made just a regular function but that would cause problems to
### people wanting to redefine alphabet() for their own objects (this is the
### case at least in the ShortRead package).
setGeneric("alphabet", function(x, ...) standardGeneric("alphabet"))

setMethod("alphabet", "ANY",
    function(x, baseOnly=FALSE)
    {
        if (!isTRUEorFALSE(baseOnly))
            stop("'baseOnly' must be TRUE or FALSE")
        switch(seqtype(x),
            DNA=if (baseOnly) DNA_BASES else DNA_ALPHABET,
            RNA=if (baseOnly) RNA_BASES else RNA_ALPHABET,
            AA=AA_ALPHABET,
            NULL
        )
    }
)

