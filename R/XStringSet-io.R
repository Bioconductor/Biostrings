### =========================================================================
### Input/output of XStringSet objects
### ------------------------------------
###
### NOTE: Only FASTA files are supported for now.
###
### Typical use:
###   file <- system.file("extdata", "someORF.fa", package="Biostrings")
###   v <- read.DNAStringSet(file(file), "fasta", "DNAString")
###   write.XStringSet(v, file="someORF2.fa", "fasta")
###
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Conversion between a list of FASTA records (as one returned by
### readFASTA) and a named character vector.
###

FASTArecordsToCharacter <- function(FASTArecs, use.names=TRUE)
{
    use.names <- normalize.use.names(use.names)
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
###   XStringSetToFASTArecords(BStringSet(FASTArecordsToBStringViews(FASTArecs)))
###
### is identical to 'FASTArecs'.
### But it is NOT the case that any BStringViews object y can
### be "reconstructed" with
###
###   FASTArecordsToBStringViews(XStringSetToFASTArecords(BStringSet(y)))
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

XStringSetToFASTArecords <- function(x)
{
    if (!is(x, "XStringSet"))
        stop("'x' must be an XStringSet object")
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

### WORK IN PROGRESS!
###
### Attempt to implement a _fast_ version of the above .read.fasta() by:
###   - shortcutting the use of readFASTA()
###   - use readLines to read the FASTA file line by line
###   - call XRaw.write() on each line + XRaw() and XRaw.copy() when it's
###     time to reallocate a biggest XRaw object.
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
    if (!is.character(file) || length(file) != 1 || is.na(file))
        stop("'file' must be a non-NA character string")
    filesize <- file.info(file)$size
    if (is.na(filesize))
        stop(file, ": file not found")
    filesize <- as.integer(filesize)
    if (is.na(filesize))
        stop(file, ": file is too big")
    file <- file(file, "r")
    on.exit(close(file))
    data <- XRaw(filesize)
    subject <- new(subjectClass, data, 0L, length(data))
    subject@length <- 0L

    #ans <- XRaw.loadFASTA(subject@data, file, collapse, enc_lkup=enc_lkup(x))

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
    ans <- adjacentViews(subject, width, gapwidth=nchar(collapse, type="bytes"))
    desc(ans) <- desc
    ans
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
### The "read.BStringSet", "read.DNAStringSet", "read.RNAStringSet" and
### "read.AAStringSet" functions.
###

read.BStringSet <- function(file, format)
    BStringSet(read.BStringViews(file, format, "BString"))

read.DNAStringSet <- function(file, format)
    DNAStringSet(read.BStringViews(file, format, "DNAString"))

read.RNAStringSet <- function(file, format)
    RNAStringSet(read.BStringViews(file, format, "RNAString"))

read.AAStringSet <- function(file, format)
    AAStringSet(read.BStringViews(file, format, "AAString"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "write.XStringSet" and "write.BStringViews" functions.
###

.write.fasta <- function(x, file, width)
{
    FASTArecs <- XStringSetToFASTArecords(x)
    writeFASTA(FASTArecs, file, width)
}

write.XStringSet <- function(x, file="", format, width=80)
{
    if (!isSingleString(format))
        stop("'format' must be a single string")
    format <- match.arg(tolower(format), c("fasta"))
    switch(format,
        "fasta"=.write.fasta(x, file, width)
    )
}

write.BStringViews <- function(x, file="", format, width=80)
{
    y <- BStringViewsToSet(x, use.names=TRUE)
    write.XStringSet(y, file=file, format, width=width)
}

