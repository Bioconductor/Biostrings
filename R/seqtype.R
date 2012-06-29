### =========================================================================
### The seqtype() generic & related functions
### -------------------------------------------------------------------------
###
### Most sequence containers in Biostrings have a "sequence type" that
### indicates the nature of the sequence(s) that the container can store:
###
###   sequence  |                           |              |
###   type      | description               | alphabet     | encoded
###   ----------|---------------------------|--------------|--------
###   "B"       | general purpose string(s) | bytes 0-255  | no
###   "DNA"     | DNA sequence(s)           | DNA_ALPHABET | yes
###   "RNA"     | RNA sequence(s)           | RNA_ALPHABET | yes
###   "AA"      | amino acid sequence(s)    | AA_ALPHABET  | no
###
### seqtype() returns that sequence type. For example 'seqtype(AAString())'
### returns "AA".
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
setGeneric("seqtype<-", signature="x",
    function(x, value) standardGeneric("seqtype<-")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions for which the returned value depends on 'seqtype(x)',
### not on what particular data are in 'x'. Not exported.
### 

xsbaseclass <- function(x) paste(seqtype(x), "String", sep="")

xscodes <- function(x, baseOnly=FALSE)
{
    switch(seqtype(x),
        DNA=DNAcodes(baseOnly),
        RNA=RNAcodes(baseOnly),
        0:255
    )
}

xscodec <- function(x)
{
    switch(seqtype(x),
        DNA=DNA_STRING_CODEC,
        RNA=RNA_STRING_CODEC,
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
### Some restrictions apply for converting from a sequence type to another or
### for comparing XString objects of different sequence types. This is
### due to the fact that XString objects with different sequence types can use
### different encodings for their sequence data (or no encoding at all) or
### simply to the fact that the conversion or comparison doesn't make sense
### from a biological perspective.
### The helper functions below are used internally (they are NOT exported) to
### determine those restrictions.
###

compatible_seqtypes <- function(seqtype1, seqtype2)
{
    if (seqtype1 %in% c("DNA", "RNA"))
        return(seqtype2 != "AA")
    if (seqtype1 == "AA")
        return(!(seqtype2 %in% c("DNA", "RNA")))
    TRUE
}

### Exported.
get_seqtype_conversion_lookup <- function(from_seqtype, to_seqtype)
{
    if (!compatible_seqtypes(from_seqtype, to_seqtype))
        stop("incompatible sequence types \"",
             from_seqtype, "\" and \"", to_seqtype, "\"")
    from_nucleo <- from_seqtype %in% c("DNA", "RNA")
    to_nucleo <- to_seqtype %in% c("DNA", "RNA")
    if (from_nucleo == to_nucleo)
        return(NULL)
    if (to_seqtype == "DNA")
        return(DNA_STRING_CODEC@enc_lkup)
    if (to_seqtype == "RNA")
        return(RNA_STRING_CODEC@enc_lkup)
    if (from_seqtype == "DNA")
        return(DNA_STRING_CODEC@dec_lkup)
    if (from_seqtype == "RNA")
        return(RNA_STRING_CODEC@dec_lkup)
    stop("Biostrings internal error, please report") # should never happen
}

comparable_seqtypes <- function(seqtype1, seqtype2)
{
    is_nucleo1 <- seqtype1 %in% c("DNA", "RNA")
    is_nucleo2 <- seqtype2 %in% c("DNA", "RNA")
    is_nucleo1 == is_nucleo2
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

