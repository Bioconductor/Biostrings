### =========================================================================
### Input/output of BStringViews objects
### ------------------------------------
###
### Create a BStringViews object by reading a file or write a BStringViews
### object to a file.
###
### NOTE: Only FASTA files are supported for now.
###
### Typical use:
###   file <- system.file("Exfiles", "someORF.fsa", package="Biostrings")
###   v <- read.BStringViews(file(file), "fasta", "DNAString")
###   write.BStringViews(v, file="someORF2.fsa", "fasta")
###
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Conversion between a list of FASTA records (as one returned by
### readFASTA) and a named character vector.
###

FASTArecordsToCharacter <- function(FASTArecs, use.names=TRUE)
{
    ans <- sapply(FASTArecs, function(rec) rec$seq)
    if (use.names)
        names(ans) <- sapply(FASTArecs, function(rec) rec$desc)
    ans
}

CharacterToFASTArecords <- function(x)
{
    if (!is.character(x))
        stop("'x' must be a character vector")
    lapply(seq_len(length(x)),
           function(i) list(desc=names(x)[i], seq=x[[i]]))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Conversion between a list of FASTA records (as one returned by
### readFASTA) and a BStringViews object.
###
### Note that, for any well-formed list of FASTA records 'FASTArecs',
###
###   BStringViewsToFASTArecords(FASTArecordsToBStringViews(FASTArecs))
###
### is identical to 'FASTArecs'.
### But it is NOT the case that any BStringViews object y can
### be "reconstructed" with
###
###   FASTArecordsToBStringViews(BStringViewsToFASTArecords(y))
###

FASTArecordsToBStringViews <- function(FASTArecs, subjectClass, collapse="")
{
    if (!is.character(subjectClass) || length(subjectClass) != 1 || is.na(subjectClass))
        stop("'subjectClass' must be a single string")
    if (!is.character(collapse) || length(collapse) != 1 || is.na(collapse))
        stop("'collapse' must be a single string")
    src <- FASTArecordsToCharacter(FASTArecs)
    BStringViews(src, subjectClass, collapse)
}

BStringViewsToFASTArecords <- function(x)
{
    if (!is(x, "BStringViews"))
        stop("'x' must be a BStringViews object")
    CharacterToFASTArecords(as.character(x))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "read.BStringViews" function.
###

.read.fasta <- function(file, subjectClass, collapse)
{
    FASTArecs <- readFASTA(file, strip.desc=TRUE)
    FASTArecordsToBStringViews(FASTArecs, subjectClass, collapse)
}

read.BStringViews <- function(file, format, subjectClass, collapse="")
{
    if (missing(file))
        stop("'file' must be specified")
    if (!is.character(format) || length(format) != 1 || is.na(format))
        stop("'format' must be a single string")
    format <- match.arg(tolower(format), c("fasta"))
    if (!is.character(subjectClass) || length(subjectClass) != 1 || is.na(subjectClass))
        stop("'subjectClass' must be a single string")
    if (!is.character(collapse) || length(collapse) != 1 || is.na(collapse))
        stop("'collapse' must be a single string")
    switch(format,
        "fasta"=.read.fasta(file, subjectClass, collapse)
    )
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "write.BStringViews" function.
###

.write.fasta <- function(x, file, width)
{
    FASTArecs <- BStringViewsToFASTArecords(x)
    writeFASTA(FASTArecs, file, width)
}

write.BStringViews <- function(x, file="", format, width=80)
{
    if (!is.character(format) || length(format) != 1 || is.na(format))
        stop("'format' must be a single string")
    format <- match.arg(tolower(format), c("fasta"))
    switch(format,
        "fasta"=.write.fasta(x, file, width)
    )
}

