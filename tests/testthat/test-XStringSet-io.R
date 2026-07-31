## TODO: Add in I/O tests for FASTQ files
# (will revisit after doing QualityScaledXStringSet testing)

check_rw_strings_fasta <- function(XStringSetRFUN, XStringSetMAKEFUN, X_ALPH) {
    # setup
    set.seed(123L)
    NUM_STRINGS <- 10L
    seq_lens <- sample(100L, NUM_STRINGS, replace=TRUE)
    f <- tempfile()
    if (file.exists(f)) unlink(f)

    # make some test strings
    strings <- character(NUM_STRINGS)
    for (i in seq_len(NUM_STRINGS))
        strings[i] <- paste(sample(X_ALPH, seq_lens[i], replace=TRUE),
                            collapse="")

    # make an XStringSet and name it
    ss <- XStringSetMAKEFUN(strings)
    names(ss) <- as.character(seq_len(NUM_STRINGS))

    # write, then ensure it wrote
    writeXStringSet(ss, f, format="fasta")
    expect_true(file.exists(f))

    # check that seq lengths in the fasta file are correct
    ssn <- fasta.seqlengths(f, seqtype=seqtype(ss), use.names=TRUE)
    expected_lengths <- lengths(ss)
    names(expected_lengths) <- names(ss)
    expect_equal(ssn, expected_lengths)

    # check that optional parameters work
    object <- fasta.seqlengths(f, seqtype=seqtype(ss), nrec=2, use.names=FALSE)
    expect_equal(object, seq_lens[seq_len(2)])
    object <- fasta.seqlengths(f, seqtype=seqtype(ss), skip=2, nrec=3,
                               use.names=FALSE)
    expect_equal(object, seq_lens[3:5])

    # check indexing
    indices <- fasta.index(f, seqtype=seqtype(ss))
    indices_sub <- fasta.index(f, skip=2L, nrec=3L, seqtype=seqtype(ss))

    expect_s3_class(indices, "data.frame")
    expected <- c("recno", "fileno", "offset", "desc", "seqlength", "filepath")
    expect_equal(colnames(indices), expected)
    expect_true(all(indices$filepath == f))
    expect_true(all(indices$recno == names(ss)))
    expect_equal(indices$seqlength, seq_lens)
    expect_equal(nrow(indices_sub), 3L)

    # rownames will screw this up
    rownames(indices_sub) <- 3:5
    expect_equal(indices[3:5,], indices_sub)

    # check that the actual values are correct
    ss_read <- XStringSetRFUN(f, format="fasta")
    ss_read_sub <- XStringSetRFUN(f, nrec=3, skip=2)

    # double checking with equals.XStringSet and standard character
    expect_equal(class(ss_read), class(ss))
    expect_true(all(ss_read == ss))
    expect_equal(as.character(ss_read), as.character(ss))
    expect_true(all(ss_read_sub == ss[3:5]))
    expect_equal(as.character(ss_read_sub), as.character(ss[3:5]))

    expect_equal(names(ss_read), names(ss))
    expect_equal(names(ss_read_sub), names(ss)[3:5])

    # check that we can append to file
    writeXStringSet(ss, f, format="fasta", append=TRUE)
    ss_read <- XStringSetRFUN(f, format="fasta")
    expect_true(all(ss_read == c(ss, ss)))
    expect_equal(as.character(ss_read), as.character(c(ss, ss)))

    unlink(f)
}

check_rw_strings_serial <- function(XStringSetMAKEFUN, X_ALPH) {
    # setup
    set.seed(345L)
    NUM_STRINGS <- 10L
    seq_lens <- sample(100L, NUM_STRINGS, replace=TRUE)
    d <- tempdir()
    f <- "example_xs"
    outfile <- file.path(d, paste0(f, ".rda"))
    if (file.exists(outfile)) unlink(outfile)

    # make some test strings
    strings <- character(NUM_STRINGS)
    for (i in seq_len(NUM_STRINGS))
        strings[i] <- paste(sample(X_ALPH, seq_lens[i], replace=TRUE),
                            collapse="")

    ss <- XStringSetMAKEFUN(strings)

    expect_output(saveXStringSet(ss, f, dirpath=d, save.dups=FALSE),
                  "Saving .+/example_xs.rda")
    expect_true(file.exists(outfile))

    if ("example_xs" %in% ls()) rm("example_xs")
    load(outfile)
    expect_equal(as.character(example_xs), as.character(ss))
    expect_true(all(example_xs == ss))
    unlink(outfile)

    # unsure how this is supposed to work
    #object <- saveXStringSet(XStringSetMAKEFUN(rep(strings, 2)),
    #                         f, dirpath=d, save.dups=TRUE)
    #expect_output(object, "Saving .+/example_xs.rda")

    ## testing standard save()
    old_ss <- ss
    save(ss, file=outfile)
    expect_true(file.exists(outfile))
    rm(ss)
    load(outfile)
    expect_equal(as.character(ss), as.character(old_ss))
    expect_true(all(ss == old_ss))
    unlink(outfile)
}

test_that("All XStringSets can be read/written to a FASTA", {
    B_ALPHABET <- strsplit(rawToChar(as.raw(32:126)), "")[[1L]]
    check_rw_strings_fasta(readDNAStringSet, DNAStringSet, DNA_ALPHABET)
    check_rw_strings_fasta(readRNAStringSet, RNAStringSet, RNA_ALPHABET)
    check_rw_strings_fasta(readAAStringSet, AAStringSet, AA_ALPHABET)
    check_rw_strings_fasta(readBStringSet, BStringSet, B_ALPHABET)
})

test_that("All XStringSets can be serialized and subseqeuntly read", {
    B_ALPHABET <- strsplit(rawToChar(as.raw(32:126)), "")[[1L]]
    check_rw_strings_serial(DNAStringSet, DNA_ALPHABET)
    check_rw_strings_serial(RNAStringSet, RNA_ALPHABET)
    check_rw_strings_serial(AAStringSet, AA_ALPHABET)
    check_rw_strings_serial(BStringSet, B_ALPHABET)
})

## Previous RUnit testing
test_that("XStringSet read/write is correct", {
    ## setup initial values
    tmpdir <- tempdir()
    tmpfilefa <- file.path(tmpdir,"whee.fa")
    tmpfilefq <- file.path(tmpdir,"whee.fq")
    if (file.exists(tmpfilefa)) unlink(tmpfilefa)
    if (file.exists(tmpfilefq)) unlink(tmpfilefq)

    x1 <- DNAStringSet(c("TTGA", "CTCN"))
    q1 <- PhredQuality(c("*+,-", "6789"))
    qdna1 <- QualityScaledDNAStringSet(x1, q1)
    names_dna <- c("A", "B")
    names(qdna1) <- names_dna

    names(x1) <- names_dna
    writeXStringSet(x1, format="fastq", file=tmpfilefq)
    qdna2 <- readQualityScaledDNAStringSet(tmpfilefq)
    expect_true(!all(quality(qdna2) == quality(qdna1)))
    expect_true(all(all(as(quality(qdna2),"IntegerList") == 26L)))

    writeXStringSet(x1, format="fastq", file=tmpfilefq, qualities = q1)
    qdna2 <- readQualityScaledDNAStringSet(tmpfilefq)
    expect_true(all(quality(qdna2) == quality(qdna1)))

    writeXStringSet(qdna1, format="fastq", file=tmpfilefq)
    qdna2 <- readQualityScaledDNAStringSet(tmpfilefq)
    expect_true(all(names(qdna2) == names_dna))
    expect_true(all(quality(qdna2) == quality(qdna1)))

    qdna2 <- readDNAStringSet(tmpfilefq, format="fastq")
    expect_true(is.null(mcols(qdna2)))

    qdna2 <- readDNAStringSet(tmpfilefq, format="fastq", with.qualities = TRUE)
    expect_true(all(mcols(qdna2)$qualities == quality(qdna1)))

    unlink(tmpfilefa)
    unlink(tmpfilefq)
})

test_that("Writing XStringSets larger than IOBUF_SIZE functions correctly", {
    ## IOBUF_SIZE is 200,003
    tf <- tempfile()
    set.seed(456L)
    x1 <- paste(sample(DNA_BASES, 400020, replace=TRUE), collapse='')
    x2 <- paste(sample(DNA_BASES, 200127, replace=TRUE), collapse='')
    xss <- DNAStringSet(c(x1, x2))
    names(xss) <- c("seq1", "seq2")
    exp_output_fasta <- paste0('>seq1', x1, '>seq2', x2)
    for(w in c(5, 80, 213, 200003, width(xss)[1], width(xss)[2], -1)){
      writeXStringSet(xss, tf, width = w)
      opt <- readLines(tf)
      opt <- paste(opt, collapse='')
      expect_equal(opt, exp_output_fasta)
      read_seqs <- readDNAStringSet(tf)
      expect_true(all(as.character(read_seqs) == as.character(xss)))
    }
})

test_that("Writing QualityScaledXStringSets larger than IOBUF_SIZE functions correctly", {
    ## IOBUF_SIZE is 200,003
    tf <- tempfile()
    set.seed(456L)
    x1 <- paste(sample(DNA_BASES, 200127, replace=TRUE), collapse='')
    x2 <- paste(sample(DNA_BASES, 400020, replace=TRUE), collapse='')
    q1 <- paste(sample(as.character(1:9), 200127, replace=TRUE), collapse='')
    q2 <- paste(sample(as.character(1:9), 400020, replace=TRUE), collapse='')
    xss <- DNAStringSet(c(x1, x2))
    quals <- PhredQuality(c(q1, q2))
    names(xss) <- c("seq1", "seq2")
    names(quals) <- c("seq1", "seq2")
    exp_output_fastq <- paste0("@seq1", x1, "+seq1", q1, "@seq2", x2, "+seq2", q2)
    qss <- QualityScaledDNAStringSet(xss, quals)
    writeQualityScaledXStringSet(qss, tf)
    opt <- paste(readLines(tf), collapse='')
    expect_equal(opt, exp_output_fastq)
    read_seqs <- readQualityScaledDNAStringSet(tf)
    expect_true(all(as.character(read_seqs) == as.character(xss)))
    expect_true(all(as.character(quality(read_seqs)) == as.character(quals)))
})
