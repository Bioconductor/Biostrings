### =========================================================================
### Input/output of XStringSet objects
### -------------------------------------------------------------------------


### 'filexp_list' must be a list of "file external pointers" returned by
### XVector:::open_input_files() or XVector:::open_output_file().
.finalize_filexp_list <- function(filexp_list)
{
    for (filexp in filexp_list) XVector:::finalize_filexp(filexp)
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
    filexp_list <- XVector:::open_input_files(filepath)
    on.exit(.finalize_filexp_list(filexp_list))
    nrec <- .normarg_nrec(nrec)
    skip <- .normarg_skip(skip)
    if (!isTRUEorFALSE(seek.first.rec)) 
        stop(wmsg("'seek.first.rec' must be TRUE or FALSE"))
    seqtype <- match.arg(seqtype, c("B", "DNA", "RNA", "AA"))
    lkup <- get_seqtype_conversion_lookup("B", seqtype)
    ans <- .Call2("fasta_index",
                  filexp_list, nrec, skip, seek.first.rec, lkup,
                  PACKAGE="Biostrings")
    ## 'expath' will usually be the same as 'filepath', except when 'filepath'
    ## contains URLs which will be replaced by the path to the downloaded file.
    expath <- vapply(filexp_list, attr, character(1), "expath", USE.NAMES=FALSE)
    ans$filepath <- expath[ans[ , "fileno"]]
    ans
}

fasta.seqlengths <- function(filepath, nrec=-1L, skip=0L, seek.first.rec=FALSE,
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
### Fasta index 'ssorted_fai' must be strictly sorted by "recno". This is NOT
### checked!
.compute_sorted_fasta_blocks_from_ssorted_fasta_index <- function(ssorted_fai)
{
    recno <- ssorted_fai[ , "recno"]
    fileno <- ssorted_fai[ , "fileno"]
    offset <- ssorted_fai[ , "offset"]
    blockid <- recno - seq_along(recno)  # this block id is unique only within
                                         # a given file
    is_first_in_block <- !duplicatedIntegerPairs(blockid, fileno)
    first_in_block_idx <- which(is_first_in_block)
    data.frame(fileno=fileno[first_in_block_idx],
               nrec=diff(c(first_in_block_idx, nrow(ssorted_fai) + 1L)),
               offset=offset[first_in_block_idx])
}

### Fasta index 'ssorted_fai' must be strictly sorted by "recno". This is NOT
### checked!
.read_XStringSet_from_ssorted_fasta_index <- function(ssorted_fai,
                                                      elementType, lkup)
{
    ## Prepare 'nrec_list' and 'offset_list'.
    fasta_blocks <-
        .compute_sorted_fasta_blocks_from_ssorted_fasta_index(ssorted_fai)
    nrec_list <- split(fasta_blocks[ , "nrec"], fasta_blocks[ , "fileno"],
                       drop=TRUE)
    offset_list <- split(fasta_blocks[ , "offset"], fasta_blocks[ , "fileno"],
                         drop=TRUE)

    ## Prepare 'filexp_list'.
    filepath <- ssorted_fai[ , "filepath"]
    fileno <- ssorted_fai[ , "fileno"]
    used_fileno <- as.integer(names(nrec_list))
    used_filepath <- filepath[match(used_fileno, fileno)]
    filexp_list <- XVector:::open_input_files(used_filepath)
    on.exit(.finalize_filexp_list(filexp_list))

    ## Prepare 'seqlength'.
    seqlength <- ssorted_fai[ , "seqlength"]

    .Call2("read_XStringSet_from_fasta_blocks",
           seqlength, filexp_list, nrec_list, offset_list,
           elementType, lkup,
           PACKAGE="Biostrings")
}

.read_XStringSet_from_fasta_index <- function(fai, use.names, elementType, lkup)
{
    .check_fasta_index(fai)

    ## Create a "strictly sorted" version of 'fai' by removing duplicated rows
    ## and sorting the remaining rows by ascending "recno".
    recno <- fai[ , "recno"]
    ssorted_recno <- sort(unique(recno))
    ssorted_fai <- fai[match(ssorted_recno, recno), , drop=FALSE]

    C_ans <- .read_XStringSet_from_ssorted_fasta_index(ssorted_fai,
                                                       elementType, lkup)

    ## Re-order XStringSet object to make it parallel to 'recno'.
    ans <- C_ans[match(recno, ssorted_recno)]

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
    filexp_list <- XVector:::open_input_files(filepath)
    on.exit(.finalize_filexp_list(filexp_list))
    nrec <- .normarg_nrec(nrec)
    skip <- .normarg_skip(skip)
    if (!isTRUEorFALSE(seek.first.rec)) 
        stop(wmsg("'seek.first.rec' must be TRUE or FALSE"))
    .Call2("fastq_geometry",
           filexp_list, nrec, skip, seek.first.rec,
           PACKAGE="Biostrings")
}

.read_XStringSet_from_fastq <- function(filepath, nrec, skip, seek.first.rec,
                                        use.names, elementType, lkup)
{
    filexp_list <- XVector:::open_input_files(filepath)
    on.exit(.finalize_filexp_list(filexp_list))
    nrec <- .normarg_nrec(nrec)
    skip <- .normarg_skip(skip)
    if (!isTRUEorFALSE(seek.first.rec)) 
        stop(wmsg("'seek.first.rec' must be TRUE or FALSE"))
    .Call2("read_XStringSet_from_fastq",
           filexp_list, nrec, skip, seek.first.rec,
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

.write_XStringSet_to_fasta <- function(x, filexp_list, width=80L)
{
    if (!isSingleNumber(width)) 
        stop(wmsg("'width' must be a single integer"))
    if (!is.integer(width)) 
        width <- as.integer(width)
    if (width < 1L) 
        stop(wmsg("'width' must be an integer >= 1"))
    lkup <- get_seqtype_conversion_lookup(seqtype(x), "B")
    .Call2("write_XStringSet_to_fasta",
          x, filexp_list, width, lkup,
          PACKAGE="Biostrings")
}

.write_XStringSet_to_fastq <- function(x, filexp_list, qualities=NULL)
{
    if (!is.null(qualities)) {
        if (!is(qualities, "BStringSet"))
            stop(wmsg("'qualities' must be NULL or a BStringSet object"))
        if (length(qualities) != length(x))
            stop(wmsg("'x' and 'qualities' must have the same length"))
    }
    lkup <- get_seqtype_conversion_lookup(seqtype(x), "B")
    .Call2("write_XStringSet_to_fastq",
          x, filexp_list, qualities, lkup,
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
    filexp_list <- XVector:::open_output_file(filepath, append,
                                              compress, compression_level)
    res <- try(switch(format,
                   "fasta"=.write_XStringSet_to_fasta(x, filexp_list, ...),
                   "fastq"=.write_XStringSet_to_fastq(x, filexp_list, ...)
               ),
               silent=FALSE)
    .finalize_filexp_list(filexp_list)
    if (is(res, "try-error") && !append) {
        ## Get the expamded path and remove the file.
        expath <- attr(filexp_list[[1L]], "expath")
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
        if (!is(pdict, "try-error") && !is.null(pdict@dups0)) {
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

