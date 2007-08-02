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
### BStringViews <-> FASTA file
###

### Conversion between a "FASTA records object" (a list as one returned by
### readFASTA) and a BStringViews object.
###
### Note that, for any well-formed "FASTA records object" FASTArecs,
###
###   BStringViewsToFASTArecords(FASTArecordsToBStringViews(FASTArecs))
###
### is identical to FASTArecs.
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
    src <- sapply(FASTArecs, function(rec) rec$seq)
    x <- BStringViews(src, subjectClass, collapse)
    x@views$desc <- sapply(FASTArecs, function(rec) rec$desc)
    x
}
BStringViewsToFASTArecords <- function(x)
{
    if (!is(x, "BStringViews"))
        stop("'x' must be a BStringViews object")
    FASTArecs <- list()
    for (i in seq_len(length(x))) {
        FASTArecs[[i]] <- list(desc=desc(x)[i], seq=as.character(x[[i]]))
    }
    FASTArecs
}

.read.fasta <- function(file, subjectClass, collapse)
{
    FASTArecs <- readFASTA(file)
    FASTArecordsToBStringViews(FASTArecs, subjectClass, collapse)
}

.write.fasta <- function(x, file, width)
{
    FASTArecs <- BStringViewsToFASTArecords(x)
    writeFASTA(FASTArecs, file, width)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "read.BStringViews" function.
###

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

write.BStringViews <- function(x, file="", format, width=80)
{
    if (!is.character(format) || length(format) != 1 || is.na(format))
        stop("'format' must be a single string")
    format <- match.arg(tolower(format), c("fasta"))
    switch(format,
        "fasta"=.write.fasta(x, file, width)
    )
}
