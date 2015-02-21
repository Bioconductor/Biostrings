### =========================================================================
### Input/output of XStringSet objects
### -------------------------------------------------------------------------


.normarg_input_filepath <- function(filepath)
{
    if (!is.character(filepath) || any(is.na(filepath)))
        stop(wmsg("'filepath' must be a character vector with no NAs"))
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
        if (!(filetype[i] %in% c("file", "gzfile")))
            stop(wmsg("file \"", filepath[i], "\" ",
                      "has unsupported type: ", filetype[i]))
    }
    names(filepath2) <- filetype
    filepath2
}

.finalize_ExternalFilePtr <- function(x)
    .Call2("finalize_ExternalFilePtr", x, PACKAGE="Biostrings")

### Returns a list of "external file pointers".
.open_input_files <- function(filepath)
{
    filepath2 <- .normarg_input_filepath(filepath)
    ans <- lapply(filepath2,
           function(fp)
           {
               efp <- .Call2("new_input_ExternalFilePtr", fp,
                            PACKAGE="Biostrings")
               reg.finalizer(efp, .finalize_ExternalFilePtr, onexit=TRUE)
               efp
           })
    names(ans) <- filepath
    ans
}

.normarg_compress <- function(compress)
{
    if (isTRUEorFALSE(compress)) {
        if (compress)
            return("gzip")
        return("no")
    }
    if (isSingleString(compress)) {
        # Types of compression supported by save():
        #VALID_COMPRESS <- c("no", "gzip", "bzip2", "xz")
        VALID_COMPRESS <- c("no", "gzip")
        if (!(compress %in% VALID_COMPRESS))
            stop(wmsg("when 'compress' is a single string, it must be one of ",
                      paste(paste0("\"", VALID_COMPRESS, "\""),
                            collapse=", ")))
        return(compress)
    }
    stop(wmsg("'compress' must be TRUE or FALSE or a single string"))
}

.normarg_compression_level <- function(compression_level, compress)
{
    if (!isSingleNumberOrNA(compression_level))
        stop(wmsg("'compression_level' must be a single number or NA"))
    if (is.na(compression_level))
        return(switch(compress, no=0L, gzip=6L, bzip2=9L, xz=9L))
    if (!is.integer(compression_level))
        compression_level <- as.integer(compression_level)
    if (compression_level < 0L)
        stop(wmsg("'compression_level' cannot be negative"))
    compression_level
}

### Returns a length-1 list of "external file pointers".
.open_output_file <- function(filepath, append, compress, compression_level)
{
    if (!isSingleString(filepath))
        stop(wmsg("'filepath' must be a single string"))
    if (!isTRUEorFALSE(append))
        stop(wmsg("'append' must be TRUE or FALSE"))
    compress <- .normarg_compress(compress)
    compression_level <- .normarg_compression_level(compression_level, compress)
    filepath2 <- path.expand(filepath)
    efp <- .Call2("new_output_ExternalFilePtr",
                  filepath2, append, compress, compression_level,
                  PACKAGE="Biostrings")
    reg.finalizer(efp, .finalize_ExternalFilePtr, onexit=TRUE)
    ans <- list(efp)
    names(ans) <- filepath
    ans
}

### 'efp_list' must be a list of "external file pointers" returned by
### .open_input_files() or .open_output_files().
.finalize_efp_list <- function(efp_list)
{
    for (efp in efp_list) .finalize_ExternalFilePtr(efp)
}

.normarg_nrec <- function(nrec)
{
    if (!isSingleNumber(nrec))
        stop(wmsg("'nrec' must be a single integer value"))
    if (!is.integer(nrec))
        nrec <- as.integer(nrec)
    nrec
}

.normarg_skip <- function(skip)
{
    if (!isSingleNumber(skip))
        stop(wmsg("'skip' must be a single integer value"))
    if (!is.integer(skip))
        skip <- as.integer(skip)
    if (skip < 0L)
        stop(wmsg("'skip' cannot be negative"))
    skip
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### FASTA
###

fasta.index <- function(filepath, nrec=-1L, skip=0L, seek.first.rec=FALSE,
                        seqtype="B")
{
    efp_list <- .open_input_files(filepath)
    on.exit(.finalize_efp_list(efp_list))
    nrec <- .normarg_nrec(nrec)
    skip <- .normarg_skip(skip)
    if (!isTRUEorFALSE(seek.first.rec)) 
        stop(wmsg("'seek.first.rec' must be TRUE or FALSE"))
    seqtype <- match.arg(seqtype, c("B", "DNA", "RNA", "AA"))
    lkup <- get_seqtype_conversion_lookup("B", seqtype)
    ans <- .Call2("fasta_index",
                  efp_list, nrec, skip, seek.first.rec, lkup,
                  PACKAGE="Biostrings")
    ans$filepath <- filepath[ans[ , "fileno"]]
    ans
}

### TODO (maybe): Rename this fasta.seqlengths() and deprecate fasta.info()
fasta.info <- function(filepath, nrec=-1L, skip=0L, seek.first.rec=FALSE,
                       seqtype="B", use.names=TRUE)
{
    if (!isTRUEorFALSE(use.names))
        stop(wmsg("'use.names' must be TRUE or FALSE"))
    fai <- fasta.index(filepath, nrec=nrec, skip=skip,
                       seek.first.rec=seek.first.rec,
                       seqtype=seqtype)
    ans <- fai[ , "seqlength"]
    if (use.names)
        names(ans) <- fai[ , "desc"]
    ans
}

.check_fasta_index <- function(fai)
{
    .REQUIRED_COLS <- c("recno", "fileno", "offset",
                        "desc", "seqlength", "filepath")
    if (!all(.REQUIRED_COLS %in% colnames(fai)))
        stop(wmsg("invalid FASTA index: a FASTA index must be a data frame ",
                  "with columns: ", paste0(.REQUIRED_COLS, collapse=", ")))

    recno <- fai[ , "recno"]
    if (!is.integer(recno)
     || S4Vectors:::anyMissingOrOutside(recno, lower=1L))
        stop(wmsg("invalid FASTA index: the \"recno\" column must be ",
                  "an integer vector with no NAs and with positive values"))

    fileno <- fai[ , "fileno"]
    if (!is.integer(fileno)
     || S4Vectors:::anyMissingOrOutside(fileno, lower=1L))
        stop(wmsg("invalid FASTA index: the \"fileno\" column must be ",
                  "an integer vector with no NAs and with positive values"))

    offset <- fai[ , "offset"]
    if (!is.numeric(offset) || any(is.na(offset)) || any(offset < 0))
        stop(wmsg("invalid FASTA index: the \"offset\" column must be ",
                  "a numeric vector with no NAs and no negative values"))

    desc <- fai[ , "desc"]
    if (!is.character(desc) || any(is.na(desc)))
        stop(wmsg("invalid FASTA index: the \"desc\" column must be ",
                  "a character vector with no NAs"))

    seqlength <- fai[ , "seqlength"]
    if (!is.integer(seqlength)
     || S4Vectors:::anyMissingOrOutside(seqlength, lower=0L))
        stop(wmsg("invalid FASTA index: the \"seqlength\" column must be ",
                  "an integer vector with no NAs and no negative values"))

    filepath <- fai[ , "filepath"]
    if (!is.character(filepath) || any(is.na(filepath)))
        stop(wmsg("invalid FASTA index: the \"filepath\" column must be ",
                  "a character vector with no NAs"))
}

### "FASTA blocks" are groups of consecutive FASTA records.
### Fasta index 'norm_fai' must already be normalized (i.e. no duplicated rows
### and ordered by ascending recno). This is NOT checked!
.compute_sorted_fasta_blocks_from_normalized_fasta_index <- function(norm_fai)
{
    recno <- norm_fai[ , "recno"]
    fileno <- norm_fai[ , "fileno"]
    offset <- norm_fai[ , "offset"]
    blockid <- recno - seq_along(recno)  # this block id is unique only within
                                         # a given file
    is_first_in_block <- !S4Vectors:::duplicatedIntegerPairs(blockid, fileno)
    first_in_block_idx <- which(is_first_in_block)
    data.frame(fileno=fileno[first_in_block_idx],
               nrec=diff(c(first_in_block_idx, nrow(norm_fai) + 1L)),
               offset=offset[first_in_block_idx])
}

### Fasta index 'norm_fai' must already be normalized (i.e. no duplicated rows
### and ordered by ascending recno). This is NOT checked!
.read_XStringSet_from_normalized_fasta_index <- function(norm_fai,
                                                         elementType, lkup)
{
    ## Prepare 'nrec_list' and 'offset_list'.
    fasta_blocks <-
        .compute_sorted_fasta_blocks_from_normalized_fasta_index(norm_fai)
    nrec_list <- split(fasta_blocks[ , "nrec"], fasta_blocks[ , "fileno"],
                       drop=TRUE)
    offset_list <- split(fasta_blocks[ , "offset"], fasta_blocks[ , "fileno"],
                         drop=TRUE)

    ## Prepare 'efp_list'.
    filepath <- norm_fai[ , "filepath"]
    fileno <- norm_fai[ , "fileno"]
    used_fileno <- as.integer(names(nrec_list))
    used_filepath <- filepath[match(used_fileno, fileno)]
    efp_list <- .open_input_files(used_filepath)
    on.exit(.finalize_efp_list(efp_list))

    ## Prepare 'seqlength'.
    seqlength <- norm_fai[ , "seqlength"]

    .Call2("read_XStringSet_from_fasta_blocks",
           seqlength, efp_list, nrec_list, offset_list,
           elementType, lkup,
           PACKAGE="Biostrings")
}

.read_XStringSet_from_fasta_index <- function(fai, use.names, elementType, lkup)
{
    .check_fasta_index(fai)

    ## Normalize the FASTA index (i.e. remove duplicated rows and order by
    ## ascending recno).
    recno <- fai[ , "recno"]
    norm_recno <- sort(unique(recno))
    norm_fai <- fai[match(norm_recno, recno), , drop=FALSE]

    C_ans <- .read_XStringSet_from_normalized_fasta_index(norm_fai,
                                                          elementType, lkup)

    ## Re-order XStringSet object so it's parallel to 'recno'.
    ans <- C_ans[match(recno, norm_recno)]

    ## Set names on XStringSet object.
    if (use.names)
        names(ans) <- fai[ , "desc"]
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### FASTQ
###

fastq.geometry <- function(filepath, nrec=-1L, skip=0L, seek.first.rec=FALSE)
{
    efp_list <- .open_input_files(filepath)
    on.exit(.finalize_efp_list(efp_list))
    nrec <- .normarg_nrec(nrec)
    skip <- .normarg_skip(skip)
    if (!isTRUEorFALSE(seek.first.rec)) 
        stop(wmsg("'seek.first.rec' must be TRUE or FALSE"))
    .Call2("fastq_geometry",
           efp_list, nrec, skip, seek.first.rec,
           PACKAGE="Biostrings")
}

.read_XStringSet_from_fastq <- function(filepath, nrec, skip, seek.first.rec,
                                        use.names, elementType, lkup)
{
    efp_list <- .open_input_files(filepath)
    on.exit(.finalize_efp_list(efp_list))
    nrec <- .normarg_nrec(nrec)
    skip <- .normarg_skip(skip)
    if (!isTRUEorFALSE(seek.first.rec)) 
        stop(wmsg("'seek.first.rec' must be TRUE or FALSE"))
    .Call2("read_XStringSet_from_fastq",
           efp_list, nrec, skip, seek.first.rec,
           use.names, elementType, lkup,
           PACKAGE="Biostrings")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The readBStringSet(), readDNAStringSet(), readRNAStringSet(), and
### readAAStringSet() functions.
###

.read_XStringSet <- function(filepath, format,
                             nrec=-1L, skip=0L, seek.first.rec=FALSE,
                             use.names=TRUE, seqtype="B")
{
    if (!isSingleString(format))
        stop(wmsg("'format' must be a single string"))
    format <- match.arg(tolower(format), c("fasta", "fastq"))
    if (!isTRUEorFALSE(use.names)) 
        stop(wmsg("'use.names' must be TRUE or FALSE"))
    elementType <- paste(seqtype, "String", sep="")
    lkup <- get_seqtype_conversion_lookup("B", seqtype)

    ## Read FASTQ.
    if (format == "fastq") {
        ans <- .read_XStringSet_from_fastq(filepath,
                                           nrec, skip, seek.first.rec,
                                           use.names, elementType, lkup)
        return(ans)
    }

    ## Read FASTA.
    if (is.data.frame(filepath)) {
        if (!(identical(nrec, -1L) &&
              identical(skip, 0L) &&
              identical(seek.first.rec, FALSE)))
            warning(wmsg("'nrec', 'skip', and 'seek.first.rec' are ",
                         "ignored when 'filepath' is a data frame"))
        fai <- filepath
    } else {
        fai <- fasta.index(filepath, nrec=nrec, skip=skip,
                           seek.first.rec=seek.first.rec,
                           seqtype=seqtype)
    }
    .read_XStringSet_from_fasta_index(fai, use.names, elementType, lkup)
}

readBStringSet <- function(filepath, format="fasta",
                           nrec=-1L, skip=0L, seek.first.rec=FALSE,
                           use.names=TRUE)
    .read_XStringSet(filepath, format, nrec, skip, seek.first.rec,
                     use.names, "B")

readDNAStringSet <- function(filepath, format="fasta",
                             nrec=-1L, skip=0L, seek.first.rec=FALSE,
                             use.names=TRUE)
    .read_XStringSet(filepath, format, nrec, skip, seek.first.rec,
                     use.names, "DNA")

readRNAStringSet <- function(filepath, format="fasta",
                             nrec=-1L, skip=0L, seek.first.rec=FALSE,
                             use.names=TRUE)
    .read_XStringSet(filepath, format, nrec, skip, seek.first.rec,
                     use.names, "RNA")

readAAStringSet <- function(filepath, format="fasta",
                            nrec=-1L, skip=0L, seek.first.rec=FALSE,
                            use.names=TRUE)
    .read_XStringSet(filepath, format, nrec, skip, seek.first.rec,
                     use.names, "AA")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### writeXStringSet()
###

.write_XStringSet_to_fasta <- function(x, efp_list, width=80L)
{
    if (!isSingleNumber(width)) 
        stop(wmsg("'width' must be a single integer"))
    if (!is.integer(width)) 
        width <- as.integer(width)
    if (width < 1L) 
        stop(wmsg("'width' must be an integer >= 1"))
    lkup <- get_seqtype_conversion_lookup(seqtype(x), "B")
    .Call2("write_XStringSet_to_fasta",
          x, efp_list, width, lkup,
          PACKAGE="Biostrings")
}

.write_XStringSet_to_fastq <- function(x, efp_list, qualities=NULL)
{
    if (!is.null(qualities)) {
        if (!is(qualities, "BStringSet"))
            stop(wmsg("'qualities' must be NULL or a BStringSet object"))
        if (length(qualities) != length(x))
            stop(wmsg("'x' and 'qualities' must have the same length"))
    }
    lkup <- get_seqtype_conversion_lookup(seqtype(x), "B")
    .Call2("write_XStringSet_to_fastq",
          x, efp_list, qualities, lkup,
          PACKAGE="Biostrings")
}

writeXStringSet <- function(x, filepath, append=FALSE,
                            compress=FALSE, compression_level=NA,
                            format="fasta", ...)
{
    if (!is(x, "XStringSet"))
        stop(wmsg("'x' must be an XStringSet object"))
    if (!isSingleString(format))
        stop(wmsg("'format' must be a single string"))
    format <- match.arg(tolower(format), c("fasta", "fastq"))
    efp_list <- .open_output_file(filepath, append, compress, compression_level)
    res <- try(switch(format,
                   "fasta"=.write_XStringSet_to_fasta(x, efp_list, ...),
                   "fastq"=.write_XStringSet_to_fastq(x, efp_list, ...)
               ),
               silent=FALSE)
    .finalize_efp_list(efp_list)
    if (is(res, "try-error") && !append) {
        ## Get the expamded path and remove the file.
        expath <- attr(efp_list[[1L]], "expath")
        if (!file.remove(expath))
            warning(wmsg("cannot remove file '", expath, "'"))
    }
    invisible(NULL)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Serialization of XStringSet objects.
###

saveXStringSet <- function(x, objname, dirpath=".",
                           save.dups=FALSE, verbose=TRUE)
{
    if (!is(x, "XStringSet"))
        stop(wmsg("'x' must be an XStringSet object"))
    if (!isSingleString(objname))
        stop(wmsg("'objname' must be a single string"))
    if (!isSingleString(dirpath))
        stop(wmsg("'dirpath' must be a single string"))
    if (!isTRUEorFALSE(save.dups))
        stop(wmsg("'save.dups' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(verbose))
        stop(wmsg("'verbose' must be TRUE or FALSE"))
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
            stop(wmsg("could not determine 'x_dups'"))
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

