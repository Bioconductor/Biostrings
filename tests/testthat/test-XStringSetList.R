## XStringSetList-class.R exports the following:
## - DNAStringSetList
## - AAStringSetList
## note that RNAStringSetList and BStringSetList are not supported


test_that("XStringSetList constuctors work properly", {
	seq_lengths <- c(2,3,4)
	seqs <- list(DNA=paste(DNA_ALPHABET, collapse=''),
							RNA=paste(RNA_ALPHABET, collapse=''),
							AA=paste(AA_ALPHABET, collapse=''),
							B=rawToChar(as.raw(32:126)))

	for(s in c("DNA", "RNA", "AA", "B")){
		seq <- seqs[[s]]
		all_seqs <- lapply(seq_along(seq_lengths), \(i){
			get(paste0(s, "StringSet"))(rep(seq, seq_lengths[i]))
		})
		cl_name <- paste0(s, "StringSetList")
		CONSTRUCTOR <- get(cl_name)
		expect_s4_class(CONSTRUCTOR(all_seqs), cl_name)

		seqlist <- CONSTRUCTOR(all_seqs)
		expect_equal(length(seqlist), length(seq_lengths))
		expect_equal(lengths(seqlist), seq_lengths)
		expect_equal(as.list(seqlist), all_seqs)
		expect_equal(seqtype(seqlist), s)

		# empty constructor
		expect_equal(length(CONSTRUCTOR()), 0L)

		names(all_seqs) <- LETTERS[seq_along(seq_lengths)]
		expect_equal(names(CONSTRUCTOR(all_seqs)), names(all_seqs))
		expect_null(names(CONSTRUCTOR(all_seqs, use.names=FALSE)))
	}
})