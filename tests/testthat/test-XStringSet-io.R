## TODO: Add in I/O tests for FASTQ files
# (will revisit after doing QualityScaledXStringSet testing)

check_rw_strings_fasta <- function(XStringSetRFUN, XStringSetMAKEFUN, X_ALPH){
	# setup
	set.seed(123L)
	NUM_STRINGS <- 10L
	seq_lens <- sample(100L, NUM_STRINGS, replace=TRUE)
	f <- tempfile()
	if(file.exists(f)) unlink(f)

	# make some test strings
	strings <- character(NUM_STRINGS)
	for(i in seq_len(NUM_STRINGS))
		strings[i] <- paste(sample(X_ALPH, seq_lens[i], replace=TRUE), collapse='')

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
	expect_equal(fasta.seqlengths(f, seqtype=seqtype(ss), nrec=2, use.names=FALSE), seq_lens[seq_len(2)])
	expect_equal(fasta.seqlengths(f, seqtype=seqtype(ss), skip=2, nrec=3, use.names=FALSE), seq_lens[3:5])

	# check indexing
	indices <- fasta.index(f, seqtype=seqtype(ss))
	indices_sub <- fasta.index(f, skip=2L, nrec=3L, seqtype=seqtype(ss))

	expect_s3_class(indices, "data.frame")
	expect_equal(colnames(indices), c("recno", "fileno", "offset", "desc", "seqlength", "filepath"))
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

check_rw_strings_serial <- function(XStringSetMAKEFUN, X_ALPH){
	# setup
	set.seed(345L)
	NUM_STRINGS <- 10L
	seq_lens <- sample(100L, NUM_STRINGS, replace=TRUE)
	d <- tempdir()
	f <- "example_xs"
	outfile <- file.path(d, paste0(f, '.rda'))
	if(file.exists(outfile)) unlink(outfile)

	# make some test strings
	strings <- character(NUM_STRINGS)
	for(i in seq_len(NUM_STRINGS))
		strings[i] <- paste(sample(X_ALPH, seq_lens[i], replace=TRUE), collapse='')

	ss <- XStringSetMAKEFUN(strings)

	expect_output(saveXStringSet(ss, f, dirpath=d, save.dups=FALSE), "Saving .+/example_xs.rda")
	expect_true(file.exists(outfile))

	if("example_xs" %in% ls()) rm("example_xs")
	load(outfile)
	expect_equal(as.character(example_xs), as.character(ss))
	expect_true(all(example_xs == ss))
	unlink(outfile)

	# unsure how this is supposed to work
	# expect_output(saveXStringSet(XStringSetMAKEFUN(rep(strings, 2)), f, dirpath=d, save.dups=TRUE), "Saving .+/example_xs.rda")

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
	B_ALPHABET <- strsplit(rawToChar(as.raw(32:126)), '')[[1]]
	check_rw_strings_fasta(readDNAStringSet, DNAStringSet, DNA_ALPHABET)
	check_rw_strings_fasta(readRNAStringSet, RNAStringSet, RNA_ALPHABET)
	check_rw_strings_fasta(readAAStringSet, AAStringSet, AA_ALPHABET)
	check_rw_strings_fasta(readBStringSet, BStringSet, B_ALPHABET)
})

test_that("All XStringSets can be serialized and subseqeuntly read", {
	B_ALPHABET <- strsplit(rawToChar(as.raw(32:126)), '')[[1]]
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
  if(file.exists(tmpfilefa)) unlink(tmpfilefa)
  if(file.exists(tmpfilefq)) unlink(tmpfilefq)


  x1 <- DNAStringSet(c("TTGA", "CTCN"))
  q1 <- PhredQuality(c("*+,-", "6789"))
  qdna1 <- QualityScaledDNAStringSet(x1, q1)
  names_dna <- c("A", "B")
  names(qdna1) <- names_dna

  names(x1) <- names_dna
  writeXStringSet(x1, format="fastq", file=tmpfilefq)
  expect_warning(qdna2 <- readQualityScaledDNAStringSet(tmpfilefq),
  	"metadata columns on input DNAStringSet object were dropped")

  expect_true(!all(quality(qdna2) == quality(qdna1)))
  expect_true(all(all(as(quality(qdna2),"IntegerList") == 26L)))

  writeXStringSet(x1, format="fastq", file=tmpfilefq, qualities = q1)
  expect_warning(qdna2 <- readQualityScaledDNAStringSet(tmpfilefq),
  	"metadata columns on input DNAStringSet object were dropped")
  expect_true(all(quality(qdna2) == quality(qdna1)))

  writeXStringSet(qdna1, format="fastq", file=tmpfilefq)
  expect_warning(qdna2 <- readQualityScaledDNAStringSet(tmpfilefq),
  	"metadata columns on input DNAStringSet object were dropped")
  expect_true(all(names(qdna2) == names_dna))
  expect_true(all(quality(qdna2) == quality(qdna1)))

  qdna2 <- readDNAStringSet(tmpfilefq, format="fastq")
  expect_true(is.null(mcols(qdna2)))

  qdna2 <- readDNAStringSet(tmpfilefq, format="fastq", with.qualities = TRUE)
  expect_true(all(mcols(qdna2)$qualities == quality(qdna1)))

  unlink(tmpfilefa)
  unlink(tmpfilefq)
})
