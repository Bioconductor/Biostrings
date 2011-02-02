### =========================================================================
### Input/output of XStringSet objects
### -------------------------------------------------------------------------


.normargInputFilepath <- function(filepath)
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

.ExternalFilePtr.close <- function(x)
    .Call("ExternalFilePtr_close", x, PACKAGE="Biostrings")

### Returns a list of "external file pointers".
.openInputFiles <- function(filepath)
{
    filepath2 <- .normargInputFilepath(filepath)
    ans <- lapply(filepath2,
           function(fp)
           {
               efp <- .Call("new_input_ExternalFilePtr", fp,
                            PACKAGE="Biostrings")
               reg.finalizer(efp, .ExternalFilePtr.close, onexit=TRUE)
               efp
           })
    names(ans) <- filepath
    ans
}

### Returns a length-1 list of "external file pointers".
.openOutputFile <- function(filepath, append)
{
    if (!isSingleString(filepath))
        stop("'filepath' must be a single string")
    if (!isTRUEorFALSE(append))
        stop("'append' must be TRUE or FALSE")
    filepath2 <- path.expand(filepath)
    efp <- .Call("new_output_ExternalFilePtr", filepath2, append,
                 PACKAGE="Biostrings")
    reg.finalizer(efp, .ExternalFilePtr.close, onexit=TRUE)
    ans <- list(efp)
    names(ans) <- filepath
    ans
}

### 'efp_list' must be a list of "external file pointers" returned by
### .openInputFiles() or .openOutputFiles().
.closeFiles <- function(efp_list)
{
    for (efp in efp_list) .ExternalFilePtr.close(efp)
}

.normargNrec <- function(nrec)
{
    if (!isSingleNumber(nrec))
        stop("'nrec' must be a single integer value")
    if (!is.integer(nrec))
        nrec <- as.integer(nrec)
    nrec
}

.normargSkip <- function(skip)
{
    if (!isSingleNumber(skip))
        stop("'skip' must be a single integer value")
    if (!is.integer(skip))
        skip <- as.integer(skip)
    if (skip < 0L)
        stop("'skip' cannot be negative")
    skip
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### FASTA
###

### TODO (maybe): Rename this fasta.geometry() and deprecate fasta.info()
fasta.info <- function(filepath, nrec=-1L, skip=0L, use.names=TRUE)
{
    efp_list <- .openInputFiles(filepath)
    on.exit(.closeFiles(efp_list))
    nrec <- .normargNrec(nrec)
    skip <- .normargSkip(skip)
    use.names <- normargUseNames(use.names)
    .Call("fasta_info",
          efp_list, nrec, skip, use.names,
          PACKAGE="Biostrings")
}

.read.fasta.in.XStringSet <- function(efp_list, nrec, skip,
                                      use.names, elementType, lkup)
{
    .Call("read_fasta_in_XStringSet",
          efp_list, nrec, skip, use.names, elementType, lkup,
          PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### FASTQ
###

fastq.geometry <- function(filepath, nrec=-1L, skip=0L)
{
    efp_list <- .openInputFiles(filepath)
    on.exit(.closeFiles(efp_list))
    nrec <- .normargNrec(nrec)
    skip <- .normargSkip(skip)
    .Call("fastq_geometry", efp_list, nrec, skip, PACKAGE="Biostrings")
}

.read.fastq.in.XStringSet <- function(efp_list, nrec, skip,
                                      use.names, elementType, lkup)
{
    .Call("read_fastq_in_XStringSet",
          efp_list, nrec, skip, use.names, elementType, lkup,
          PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "read.BStringSet", "read.DNAStringSet", "read.RNAStringSet" and
### "read.AAStringSet" functions.
###

.read.XStringSet <- function(filepath, format, nrec, skip, use.names, basetype)
{
    efp_list <- .openInputFiles(filepath)
    on.exit(.closeFiles(efp_list))
    if (!isSingleString(format))
        stop("'format' must be a single string")
    format <- match.arg(tolower(format), c("fasta", "fastq"))
    nrec <- .normargNrec(nrec)
    skip <- .normargSkip(skip)
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    elementType <- paste(basetype, "String", sep="")
    lkup <- get_xsbasetypes_conversion_lookup("B", basetype)
    switch(format,
        "fasta"=.read.fasta.in.XStringSet(efp_list, nrec, skip,
                                          use.names, elementType, lkup),
        "fastq"=.read.fastq.in.XStringSet(efp_list, nrec, skip,
                                          use.names, elementType, lkup)
    )
}

read.BStringSet <- function(filepath, format="fasta",
                            nrec=-1L, skip=0L, use.names=TRUE)
    .read.XStringSet(filepath, format, nrec, skip, use.names, "B")

read.DNAStringSet <- function(filepath, format="fasta",
                              nrec=-1L, skip=0L, use.names=TRUE)
    .read.XStringSet(filepath, format, nrec, skip, use.names, "DNA")

read.RNAStringSet <- function(filepath, format="fasta",
                              nrec=-1L, skip=0L, use.names=TRUE)
    .read.XStringSet(filepath, format, nrec, skip, use.names, "RNA")

read.AAStringSet <- function(filepath, format="fasta",
                             nrec=-1L, skip=0L, use.names=TRUE)
    .read.XStringSet(filepath, format, nrec, skip, use.names, "AA")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### write.XStringSet()
###

.write.XStringSet.to.fasta <- function(x, efp_list, width=80L)
{
    if (!isSingleNumber(width)) 
        stop("'width' must be a single integer")
    if (!is.integer(width)) 
        width <- as.integer(width)
    if (width < 1L) 
        stop("'width' must be an integer >= 1")
    lkup <- get_xsbasetypes_conversion_lookup(xsbasetype(x), "B")
    .Call("write_XStringSet_to_fasta",
          x, efp_list, width, lkup,
          PACKAGE="Biostrings")
}

.write.XStringSet.to.fastq <- function(x, efp_list, qualities)
{
    if (missing(qualities))
        stop("'qualities' is missing")
    if (!is(qualities, "BStringSet"))
        stop("'qualities' must be a BStringSet object")
    if (length(qualities) != length(x))
        stop("'x' and 'qualities' must have the same length")
    lkup <- get_xsbasetypes_conversion_lookup(xsbasetype(x), "B")
    .Call("write_XStringSet_to_fastq",
          x, efp_list, qualities, lkup,
          PACKAGE="Biostrings")
}

write.XStringSet <- function(x, filepath, append=FALSE, format="fasta", ...)
{
    if (!is(x, "XStringSet"))
        stop("'x' must be an XStringSet object")
    if (!isSingleString(format))
        stop("'format' must be a single string")
    format <- match.arg(tolower(format), c("fasta", "fastq"))
    efp_list <- .openOutputFile(filepath, append)
    res <- try(switch(format,
                   "fasta"=.write.XStringSet.to.fasta(x, efp_list, ...),
                   "fastq"=.write.XStringSet.to.fastq(x, efp_list, ...)
               ),
               silent=FALSE)
    .closeFiles(efp_list)
    if (is(res, "try-error") && !append) {
        ## Get the expamded path and remove the file.
        expath <- attr(efp_list[[1L]], "expath")
        if (!file.remove(expath))
            warning("cannot remove file '", expath, "'")
    }
    invisible(NULL)
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
### Old stuff (Defunct or Deprecated).
###

FASTArecordsToCharacter <- function(FASTArecs, use.names=TRUE)
    .Defunct(msg="FASTArecordsToCharacter() is defunct with no replacement, sorry!")

CharacterToFASTArecords <- function(x)
    .Defunct(msg="CharacterToFASTArecords() is defunct with no replacement, sorry!")

FASTArecordsToXStringViews <- function(FASTArecs, subjectClass, collapse="")
    .Defunct(msg="FASTArecordsToXStringViews() is defunct with no replacement, sorry!")

XStringSetToFASTArecords <- function(x)
    .Defunct(msg="XStringSetToFASTArecords() is defunct with no replacement, sorry!")

read.XStringViews <- function(filepath, format="fasta",
                              subjectClass, collapse="")
    .Defunct(msg="read.XStringViews() is defunct: use as(read.DNAStringSet(...), \"Views\") instead")

write.XStringViews <- function(x, file="", append=FALSE,
                               format="fasta", width=80)
{
    .Defunct(msg="write.XStringViews() is defunct: use write.XStringSet(as(x, \"XStringSet\"), ...) instead")
}

