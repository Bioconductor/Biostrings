### =========================================================================
### Input/output of XStringSet objects
### -------------------------------------------------------------------------


.normargFilepath <- function(filepath)
{
    if (!is.character(filepath) || any(is.na(filepath)))
        stop("'filepath' must be a character vector with no NAs")
    ## First pass: expand local paths and download any remote file.
    filepath2 <- character(length(filepath))
    for (i in seq_len(length(filepath))) {
        fp <- filepath[i]
        con <- file(fp)
        con_class <- class(con)[1L]
        close(con)
        if (con_class == "url") {
            filepath2[i] <- tempfile()
            download.file(fp, filepath2[i])
        } else {
            filepath2[i] <- path.expand(fp)
        }
    }
    ## Second pass: check the type of the local files (all files are
    ## now local).
    filetype <- character(length(filepath2))
    for (i in seq_len(length(filepath2))) {
        fp <- filepath2[i]
        con <- file(fp)
        ## Ugly trick to get the type of 'con'. Is there a better way?
        filetype[i] <- showConnections(TRUE)[as.character(con), "class"]
        close(con)
        if (filetype[i] != "file")
            stop("file \"", filepath[i], "\" ",
                 "has unsupported type: ", filetype[i])
    }
    names(filepath2) <- filetype
    filepath2
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### FASTA
###

### TODO: Rename this fasta.geometry() and deprecate fasta.info()
fasta.info <- function(filepath, use.descs=TRUE)
{
    filepath <- .normargFilepath(filepath)
    use.descs <- normargUseNames(use.descs)
    on.exit(.Call("io_cleanup", PACKAGE="Biostrings"))
    .Call("fasta_info", filepath, use.descs, PACKAGE="Biostrings")
}

.read.fasta.in.XStringSet <- function(filepath, set.names, elementType, lkup)
{
    on.exit(.Call("io_cleanup", PACKAGE="Biostrings"))
    .Call("read_fasta_in_XStringSet",
          filepath, set.names, elementType, lkup,
          PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### FASTQ
###

fastq.geometry <- function(filepath)
{
    filepath <- .normargFilepath(filepath)
    on.exit(.Call("io_cleanup", PACKAGE="Biostrings"))
    .Call("fastq_geometry", filepath, PACKAGE="Biostrings")
}

.read.fastq.in.XStringSet <- function(filepath, set.names, elementType, lkup)
{
    on.exit(.Call("io_cleanup", PACKAGE="Biostrings"))
    .Call("read_fastq_in_XStringSet",
          filepath, set.names, elementType, lkup,
          PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "read.BStringSet", "read.DNAStringSet", "read.RNAStringSet" and
### "read.AAStringSet" functions.
###

.read.XStringSet <- function(filepath, format, set.names, basetype)
{
    filepath <- .normargFilepath(filepath)
    if (!isSingleString(format))
        stop("'format' must be a single string")
    format <- match.arg(tolower(format), c("fasta", "fastq"))
    if (!isTRUEorFALSE(set.names))
        stop("'set.names' must be TRUE or FALSE")
    elementType <- paste(basetype, "String", sep="")
    lkup <- get_xsbasetypes_conversion_lookup("B", basetype)
    switch(format,
        "fasta"=.read.fasta.in.XStringSet(filepath, set.names,
                                          elementType, lkup),
        "fastq"=.read.fastq.in.XStringSet(filepath, set.names,
                                          elementType, lkup)
    )
}

read.BStringSet <- function(filepath, format="fasta")
    .read.XStringSet(filepath, format, TRUE, "B")

read.DNAStringSet <- function(filepath, format="fasta")
    .read.XStringSet(filepath, format, TRUE, "DNA")

read.RNAStringSet <- function(filepath, format="fasta")
    .read.XStringSet(filepath, format, TRUE, "RNA")

read.AAStringSet <- function(filepath, format="fasta")
    .read.XStringSet(filepath, format, TRUE, "AA")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### write.XStringSet()
###

### TODO: Implement this in C.
.write.XStringSet.to.fasta <- function(x, file, append, width)
{
    FASTArecs <- XStringSetToFASTArecords(x)
    writeFASTA(FASTArecs, file=file, append=append, width=width)
}

.write.XStringSet.to.fastq <- function(x, file, append)
    stop("writing to a FASTQ file is not supported yet, sorry!")

write.XStringSet <- function(x, file="", append=FALSE, format="fasta", width=80)
{
    if (!isSingleString(format))
        stop("'format' must be a single string")
    format <- match.arg(tolower(format), c("fasta"))
    switch(format,
        "fasta"=.write.XStringSet.to.fasta(x, file, append, width),
        "fastq"=.write.XStringSet.to.fastq(x, file, append)
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Serialization of XStringSet objects.
###

save.XStringSet <- function(x, objname, dirpath=".",
        save.dups=FALSE, verbose=TRUE)
{
    if (!is(x, "XStringSet"))
        stop("'x' must be an XStringSet object")
    if (!isSingleString(objname))
        stop("'objname' must be a single string")
    if (!isSingleString(dirpath))
        stop("'dirpath' must be a single string")
    if (!isTRUEorFALSE(save.dups))
        stop("'save.dups' must be TRUE or FALSE")
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    x_dups <- NULL
    ## Don't use 'is(x, "DNAStringSet")' here since we only want to use the
    ## "pre-compression trick" on a DNAStringSet *instance*. There is no
    ## guarantee that using this trick on an object deriving from the
    ## DNAStringSet class won't corrupt the data stored in the extra slots!
    if (class(x) == "DNAStringSet") {
        ## The "pre-compression trick" is based on the following property.
        ## If 'x_dup2unq' is an integer vector (of the same length as 'x')
        ## that maps from duplicated to unique elements in 'x', then
        ## 'x' and 'x[x_dup2unq]' contain exactly the same *values* i.e.
        ## 'all(x == x[x_dup2unq])' is TRUE. Note that other metadata
        ## attached to 'x' like the names etc could differ though!
        ## So by replacing 'x' with 'x[x_dup2unq]' below , we actually don't
        ## modify the sequences in 'x', but the internal representation of 'x'
        ## has changed. What has changed is that duplicated elements are not
        ## duplicated in memory anymore (i.e. they all point to the same place
        ## in memory) so calling compact() on 'x' will be much more efficient.
        ## Note that, for efficiency reasons, we use PDict() to extract the
        ## 'x_dup2unq' mapping. This means that the "pre-compression trick"
        ## works only if 'x' is a rectangular DNAStringSet instance with no
        ## IUPAC ambiguity codes.
        pdict <- try(PDict(x), silent=TRUE)  
        if (!is(pdict, "try-error")) {
            x_dups <- pdict@dups0
            x_names <- names(x)
            x_dup2unq <- togroup(x_dups)
            x <- x[x_dup2unq]
            names(x) <- x_names
        }
    }
    x <- compact(x)
    assign(objname, x)
    objfile <- paste(objname, ".rda", sep="")
    filepath <- file.path(dirpath, objfile)
    if (save.dups) {
        if (is.null(x_dups))
            stop("could not determine 'x_dups'")
        objname2 <- paste(objname, "_dups", sep="")
        assign(objname2, x_dups)
        objname <- c(objname, objname2)
    }
    if (verbose)
        cat("Saving ", filepath, " ... ", sep="")
    save(list=objname, file=filepath)
    if (verbose)
        cat("OK\n")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Legacy stuff.
###

### Conversion between a list of FASTA records (as one returned by
### readFASTA) and a named character vector.
FASTArecordsToCharacter <- function(FASTArecs, use.names=TRUE)
{
    use.names <- normargUseNames(use.names)
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

### Conversion between a list of FASTA records (as one returned by
### readFASTA) and an XStringViews object.
###
### Note that, for any well-formed list of FASTA records 'FASTArecs',
###
###   XStringSetToFASTArecords(BStringSet(FASTArecordsToXStringViews(FASTArecs)))
###
### is identical to 'FASTArecs'.
### But it is NOT the case that any XStringViews object y can
### be "reconstructed" with
###
###   FASTArecordsToXStringViews(XStringSetToFASTArecords(BStringSet(y)))
FASTArecordsToXStringViews <- function(FASTArecs, subjectClass, collapse="")
{
    if (!isSingleString(subjectClass))
        stop("'subjectClass' must be a single string")
    if (!isSingleString(collapse))
        stop("'collapse' must be a single string")
    x <- FASTArecordsToCharacter(FASTArecs)
    XStringViews(x, subjectClass, collapse)
}

XStringSetToFASTArecords <- function(x)
{
    if (!is(x, "XStringSet"))
        stop("'x' must be an XStringSet object")
    CharacterToFASTArecords(as.character(x))
}

.read.fasta <- function(filepath, subjectClass, collapse)
{
    FASTArecs <- readFASTA(filepath, strip.descs=TRUE)
    FASTArecordsToXStringViews(FASTArecs, subjectClass, collapse)
}

.read.fastq <- function(filepath, drop.quality=FALSE, subjectClass="DNAString")
    stop("FASTQ format temporarily unsupported in read.XStringViews(), sorry!")

read.XStringViews <- function(filepath, format="fasta",
                              subjectClass, collapse="")
{
    if (!is.character(filepath) || any(is.na(filepath)))
        stop("'filepath' must be a character vector with no NAs")
    if (!isSingleString(format))
        stop("'format' must be a single string")
    format <- match.arg(tolower(format), c("fasta", "fastq"))
    if (!isSingleString(subjectClass))
        stop("'subjectClass' must be a single string")
    if (!isSingleString(collapse))
        stop("'collapse' must be a single string")
    switch(format,
        "fasta"=.read.fasta(filepath, subjectClass, collapse),
        "fastq"=.read.fastq(filepath, drop.quality=TRUE, subjectClass)
    )
}

write.XStringViews <- function(x, file="", append=FALSE,
                               format="fasta", width=80)
{
    y <- XStringViewsToSet(x, use.names=TRUE)
    write.XStringSet(y, file=file, append=append, format, width=width)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Old stuff (Defunct or Deprecated).
###

FASTArecordsToBStringViews <- function(FASTArecs, subjectClass, collapse="")
{
    .Deprecated("FASTArecordsToXStringViews")
    FASTArecordsToXStringViews(FASTArecs, subjectClass, collapse=collapse)
}

read.BStringViews <- function(file, format, subjectClass, collapse="")
{
    .Deprecated("read.XStringViews")
    read.XStringViews(file, format, subjectClass, collapse=collapse)
}

write.BStringViews <- function(x, file="", format, width=80)
{
    .Deprecated("write.XStringViews")
    write.XStringViews(x, file=file, format=format, width=width)
}

