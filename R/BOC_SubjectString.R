### =========================================================================
### Preprocessed Subject Strings
### -------------------------------------------------------------------------

### The most general type of preprocessed string (not much we can say about
### it)
setClass("SubjectString", representation("VIRTUAL"))

### Now a real one!
### Note that 'base1_code', 'base2_code' and 'base3_code' must be distinct
setClass("BOC_SubjectString",
    contains="SubjectString",
    representation(
        subject="DNAString",        # TODO: support "RNAString" too
        pattern_length="integer",   # A single integer e.g. 36
        base1_code="integer",       # integer belonging to DNA_BASE_CODES
        base1_OCbuffer="XRaw",      # all buffers must be of length nchar(subject) - pattern_length + 1
        base2_code="integer",
        base2_OCbuffer="XRaw",
        base3_code="integer",
        base3_OCbuffer="XRaw",
        base4_code="integer"
    )
)

### subject <- DNAString("acag-tcgatgc-NNNDR")
### boc <- new("BOC_SubjectString", subject, 4, c("A", "C", "G"))
setMethod("initialize", "BOC_SubjectString",
    function(.Object, subject, pattern_length, base_letters)
    {
        .Object@subject <- subject
        if (!is.numeric(pattern_length) || length(pattern_length) != 1 || is.na(pattern_length))
            stop("'pattern_length' must be a single integer")
        pattern_length <- as.integer(pattern_length)
        if (pattern_length < 1L || 254L < pattern_length)
            stop("'pattern_length' must be >= 1 and <= 254")
        if (pattern_length > nchar(subject))
            stop("'pattern_length' must be <= 'nchar(subject)'")
        .Object@pattern_length <- pattern_length
        if (!is.character(base_letters) || length(base_letters) != 3
         || !all(base_letters %in% names(DNA_BASE_CODES)) || any(duplicated(base_letters)))
            stop("'base_letters' must contain 3 distinct DNA base-letters")
        buf_length <- nchar(subject) - pattern_length + 1
	code1 <- DNA_BASE_CODES[base_letters[1]]
        buf1 <- XRaw(buf_length)
	code2 <- DNA_BASE_CODES[base_letters[2]]
        buf2 <- XRaw(buf_length)
	code3 <- DNA_BASE_CODES[base_letters[3]]
        buf3 <- XRaw(buf_length)
	code4 <- DNA_BASE_CODES[setdiff(names(DNA_BASE_CODES), base_letters)]
        .Call("match_BOC_preprocess",
              subject@data@xp, subject@offset, subject@length,
              pattern_length,
              code1, buf1@xp,
              code2, buf2@xp,
              code3, buf3@xp,
              code4,
              PACKAGE="Biostrings")
        .Object@base1_code <- code1
        .Object@base1_OCbuffer <- buf1
        .Object@base2_code <- code2
        .Object@base2_OCbuffer <- buf2
        .Object@base3_code <- code3
        .Object@base3_OCbuffer <- buf3
        .Object@base4_code <- code4
        .Object
    }
)

