### =========================================================================
### The xsbasetype() generic & related functions
### -------------------------------------------------------------------------
###
### Most sequence containers in Biostrings have an "XString base type" that
### indicates the nature of the sequence(s) that the container can store:
###
###   XString   |                           |              |
###   base type | description               | alphabet     | encoded
###   ----------|---------------------------|--------------|--------
###   "B"       | general purpose string(s) | bytes 0-255  | no
###   "DNA"     | DNA sequence(s)           | DNA_ALPHABET | yes
###   "RNA"     | RNA sequence(s)           | RNA_ALPHABET | yes
###   "AA"      | amino acid sequence(s)    | AA_ALPHABET  | no
###
### xsbasetype() returns that base type. For example 'xsbasetype(AAString())'
### returns "AA".
### Unless specified otherwise, things in this file are not exported.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "xsbasetype" and "xsbasetype<-" generics.
###
### xsbasetype() and `xsbasetype<-`() have methods defined for the 4 basic
### string containers: XString (single sequence), XStringSet (multiple
### sequences), XStringViews (multiple sequences) and MaskedXString (single
### sequence).
###   

### Exported.
setGeneric("xsbasetype",
    function(x) standardGeneric("xsbasetype")
)

### Exported.
setGeneric("xsbasetype<-", signature="x",
    function(x, value) standardGeneric("xsbasetype<-")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions for which the returned value depends on 'xsbasetype(x)',
### not on what particular data are in 'x'. Not exported.
### 

xsbaseclass <- function(x) paste(xsbasetype(x), "String", sep="")

xscodes <- function(x, baseOnly=FALSE)
{
    switch(xsbasetype(x),
        DNA=DNAcodes(baseOnly),
        RNA=RNAcodes(baseOnly),
        0:255
    )
}

xscodec <- function(x)
{
    switch(xsbasetype(x),
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
### Some restrictions apply for converting from an XString base type to
### another or for comparing XString objects of different base types. This is
### due to the fact that different XString base types can use different
### encodings for their data (or no encoding at all) or simply to the fact
### that the conversion or comparison doesn't make sense from a biological
### perspective.
### The helper functions below are used internally (they are NOT exported) to
### determine those restrictions.
###

compatible_xsbasetypes <- function(basetype1, basetype2)
{
    if (basetype1 %in% c("DNA", "RNA"))
        return(basetype2 != "AA")
    if (basetype1 == "AA")
        return(!(basetype2 %in% c("DNA", "RNA")))
    TRUE
}

### Exported.
get_xsbasetypes_conversion_lookup <- function(from_basetype, to_basetype)
{
    if (!compatible_xsbasetypes(from_basetype, to_basetype))
        stop("incompatible XString base types \"",
             from_basetype, "\" and \"", to_basetype, "\"")
    from_nucleo <- from_basetype %in% c("DNA", "RNA")
    to_nucleo <- to_basetype %in% c("DNA", "RNA")
    if (from_nucleo == to_nucleo)
        return(NULL)
    if (to_basetype == "DNA")
        return(DNA_STRING_CODEC@enc_lkup)
    if (to_basetype == "RNA")
        return(RNA_STRING_CODEC@enc_lkup)
    if (from_basetype == "DNA")
        return(DNA_STRING_CODEC@dec_lkup)
    if (from_basetype == "RNA")
        return(RNA_STRING_CODEC@dec_lkup)
    stop("Biostrings internal error, please report") # should never happen
}

comparable_xsbasetypes <- function(basetype1, basetype2)
{
    is_nucleo1 <- basetype1 %in% c("DNA", "RNA")
    is_nucleo2 <- basetype2 %in% c("DNA", "RNA")
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
        switch(xsbasetype(x),
            DNA=if (baseOnly) DNA_BASES else DNA_ALPHABET,
            RNA=if (baseOnly) RNA_BASES else RNA_ALPHABET,
            AA=AA_ALPHABET,
            NULL
        )
    }
)

