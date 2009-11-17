### =========================================================================
### Input/output of XStringSet objects
### -------------------------------------------------------------------------


.SharedRaw.saveFASTA <- function(x, filepath, dec_lkup=NULL)
{
    stop("not ready yet, sorry!")
}

### Return a list of 4 elements (see comments for .SharedRaw_loadFASTA() in
### src/SharedRaw_utils.c for the details).
### 'filepath' must a path to an uncompressed FASTA file. Note that,
### unlike with the file() function, it cannot an URL, '""', '"stdin"'
### or '"clipboard"'.
.SharedRaw.loadFASTA <- function(x, filepath, collapse="", enc_lkup=NULL)
{
    if (!isSingleString(filepath))
        stop("'filepath' must be a single string")
    if (!isSingleString(collapse))
        stop("'collapse' must be a single string")
    if (!is.null(enc_lkup) && !is.integer(enc_lkup))
        stop("'enc_lkup' must be an integer vector")
    filepath <- path.expand(filepath)
    .Call("SharedRaw_loadFASTA",
          x@xp, filepath, collapse, enc_lkup, PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Conversion between a list of FASTA records (as one returned by
### readFASTA) and a named character vector.
###

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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
###

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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "read.XStringViews" function.
###

.read.fasta <- function(filepath, subjectClass, collapse)
{
    FASTArecs <- readFASTA(filepath, strip.descs=TRUE)
    FASTArecordsToXStringViews(FASTArecs, subjectClass, collapse)
}

### WORK IN PROGRESS!
###
### Attempt to implement a _fast_ version of the above .read.fasta() by:
###   - shortcutting the use of readFASTA()
###   - use readLines to read the FASTA file line by line
###   - call SharedRaw.write() on each line + SharedRaw() and
###     SharedVector.copy() when it's time to reallocate a biggest
###     SharedRaw object.
###
### 'file' must be a character string or connection.
### If 'file' is a connection, then 'type' is ignored.
### If it's a character string then 'type' can be "default" (i.e. 'file' is a
### path to the file to be opened or a complete URL, or '""' or '"stdin"' or
### '"clipboard"'), "gzfile" (then 'file' should be path to a file that
### is compressed by 'gzip'), "bzfile" (then 'file' should be path to a file
### that is compressed by 'bzip2') or "unz" (not supported yet).
.read.fasta2 <- function(file, subjectClass, collapse, type="default")
{
    filesize <- NA
    if (is.character(file)) {
        if (type == "default")
            filesize <- file.info(file)$size
        file <- switch(type,
                    default = file(file, "r"),
                    gzfile = gzfile(file, "r"),
                    bzfile = bzfile(file, "r"),
                    unz = stop("type \"unz\" not supported yet")
                )
        on.exit(close(file))
    } else {
        if (!inherits(file, "connection"))
            stop("'file' must be a character string or connection")
        if (!isOpen(file)) {
            open(file, "r")
            on.exit(close(file))
        }
    }
    if (is.na(filesize)) {
        datasize <- 0
        while (length(line <- readLines(file, n=1)) != 0) {
            if (substr(line, 1, 1) != ">")
                datasize <- datasize + nchar(line, type="bytes")
        }
    } else {
        datasize <- filesize # an estimate only, should be >= real datasize
    }
    datasize
}

### 'file' must be a filesystem path (character string) to an uncompressed FASTA file.
.read.uncompressed_fasta_file <- function(file, subjectClass, collapse)
{
    if (!isSingleString(file))
        stop("'file' must be a single string")
    filesize <- file.info(file)$size
    if (is.na(filesize))
        stop(file, ": file not found")
    filesize <- as.integer(filesize)
    if (is.na(filesize))
        stop(file, ": file is too big")
    file <- file(file, "r")
    on.exit(close(file))
    shared <- SharedRaw(filesize)
    subject <- new(subjectClass, shared=shared, length=length(shared))
    subject@length <- 0L

    #ans <- .SharedRaw.loadFASTA(subject@shared, file, collapse, enc_lkup=xs_enc_lkup(x))

#-- implement this in C, from here
    width <- integer(0)
    current_width <- 0L
    desc <- character(0)
    ## Even if scan() is faster than readLines() for loading all the lines of a
    ## file, the "load-all-the-lines-in-memory" solution is slower than the
    ## "load-one-line-at-a-time" solution.
    #lines <- scan(file=file, what="", sep="\n", allowEscapes=FALSE)
    #for (lineno in seq_len(length(lines))) {
    #    line <- lines[lineno]
    lineno <- 0
    while (length(line <- readLines(file, n=1)) != 0) {
        lineno <- lineno + 1
        nbytes <- nchar(line, type="bytes")
        if (nbytes == 0L)
            next
        char0 <- substr(line, 1, 1)
        if (char0 == ";")
            next
        if (char0 != ">") {
            if (length(desc) != length(width) + 1)
                stop("in file ", file, ", line ", lineno, ": ",
                     "number of FASTA sequences doesn't match number of description lines")
            subject <- XString.write(subject, value=line)
            current_width <- current_width + nbytes
            next
        }
        desc <- c(desc, substr(line, 2, nbytes))
        if (current_width == 0L)
            next
        width <- c(width, current_width)
        current_width <- 0L
        subject <- XString.write(subject, value=collapse)
    }
    if (subject@length == 0L)
        stop(file, ": file doesn't seem to be FASTA (no data in it)")
    if (current_width != 0L)
        width <- c(width, current_width)
#--- to here
    ans <- successiveViews(subject, width, gapwidth=nchar(collapse))
    names(ans) <- desc
    ans
}

fastq.geometry <- function(filepath)
{
    if (!is.character(filepath))
        stop("'filepath' must be a character vector")
    on.exit(.Call("io_cleanup", PACKAGE="Biostrings"))
    .Call("fastq_geometry", filepath, PACKAGE="Biostrings")
}

.read.fastq <- function(filepath, drop.quality=FALSE, subjectClass="DNAString")
{
    if (!isTRUEorFALSE(drop.quality))
        stop("'drop.quality' must be TRUE or FALSE")
    if (!identical(subjectClass, "DNAString"))
        stop("'subjectClass' must be \"DNAString\"")
    on.exit(.Call("io_cleanup", PACKAGE="Biostrings"))
    C_ans <- .Call("read_fastq", filepath, drop.quality, PACKAGE="Biostrings")
    views_width <- rep.int(C_ans[[1]][2], C_ans[[1]][1])
    fastq_seqs <- successiveViews(C_ans[[2]], views_width)
    if (drop.quality)
        return(fastq_seqs)
    fastq_quals <- successiveViews(C_ans[[3]], views_width)
    return(list(seq=fastq_seqs, qual=fastq_quals))
}

read.XStringViews <- function(filepath, format, subjectClass, collapse="")
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "read.BStringSet", "read.DNAStringSet", "read.RNAStringSet" and
### "read.AAStringSet" functions.
###

.read.fasta.in.XStringSet <- function(filepath, set.names, elementType, lkup)
{
    .Call("read_fasta_in_XStringSet",
          filepath, set.names, elementType, lkup,
          PACKAGE="Biostrings")
}

.read.fastq.in.XStringSet <- function(filepath, set.names, elementType, lkup)
{
    .read.fastq(filepath, drop.quality=TRUE, elementType)
}

.read.XStringSet <- function(filepath, format, set.names, basetype)
{
    if (!is.character(filepath))
        stop("'filepath' must be a character vector")
    if (!isSingleString(format))
        stop("'format' must be a single string")
    format <- match.arg(tolower(format), c("fasta", "fastq"))
    if (!isTRUEorFALSE(set.names))
        stop("'set.names' must be TRUE or FALSE")
    elementType <- paste(basetype, "String", sep="")
    lkup <- get_xsbasetypes_conversion_lookup("B", basetype)
    switch(format,
        "fasta"=.read.fasta.in.XStringSet(filepath, set.names, elementType, lkup),
        "fastq"=.read.fastq.in.XStringSet(filepath, set.names, elementType, lkup)
    )
}

read.BStringSet <- function(filepath, format)
    .read.XStringSet(filepath, format, TRUE, "B")

read.DNAStringSet <- function(filepath, format)
    .read.XStringSet(filepath, format, TRUE, "DNA")

read.RNAStringSet <- function(filepath, format)
    .read.XStringSet(filepath, format, TRUE, "RNA")

read.AAStringSet <- function(filepath, format)
    .read.XStringSet(filepath, format, TRUE, "AA")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "write.XStringSet" and "write.XStringViews" functions.
###

.write.fasta <- function(x, file, append, width)
{
    FASTArecs <- XStringSetToFASTArecords(x)
    writeFASTA(FASTArecs, file=file, append=append, width=width)
}

write.XStringSet <- function(x, file="", append=FALSE, format, width=80)
{
    if (!isSingleString(format))
        stop("'format' must be a single string")
    format <- match.arg(tolower(format), c("fasta"))
    switch(format,
        "fasta"=.write.fasta(x, file, append, width)
    )
}

write.XStringViews <- function(x, file="", append=FALSE, format, width=80)
{
    y <- XStringViewsToSet(x, use.names=TRUE)
    write.XStringSet(y, file=file, append=append, format, width=width)
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

