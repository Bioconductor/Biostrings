## QualityScaledXStringSet.R exports the following:
## - QualityScaledXStringSet class
## - readQualityScaledXStringSet
## - writeQualityScaledXStringSet
## - quality
## - all accessor methods defined for XStringSet objects
##
## here, X is a seqtype in c("DNA", "RNA", "AA", "B")
## the only accessors with different definitions are:
## - windows
## - narrow
## - reverse
## - reverseComplement
## - show


test_that("QualityScaledXStringSet I/O works properly", {
	set.seed(4L)
	n_samp <- 40L
	n_seq <- 3L
	qs <- list(phred=PhredQuality,
							solexa=SolexaQuality,
							illumina=IlluminaQuality)

	alphs <- list(DNA=DNA_BASES,
								RNA=RNA_BASES,
								AA=AA_STANDARD,
								B=LETTERS)

	for(quality_type in c('phred', 'solexa', 'illumina')){
		q <- do.call(c, lapply(seq_len(n_seq), \(x) qs[[quality_type]](sample(40L, n_samp, replace=TRUE))))
		for(stype in c("DNA", "RNA", "AA", "B")){
			seqs <- vapply(seq_len(n_seq), \(i) {
				paste(sample(alphs[[stype]], n_samp, replace=TRUE), collapse='')
			}, character(1L))
			seqs <- as(seqs, paste0(stype, "StringSet"))
			names(seqs) <- letters[seq_len(n_seq)]

			qss <- get(paste0("QualityScaled", stype, "StringSet"))(seqs, q)

			expect_s4_class(qss, paste0("QualityScaled", stype, "StringSet"))
			expect_s4_class(quality(qss), class(qs[[quality_type]](0L)))

			tf <- tempfile()
			if(file.exists(tf)) file.remove(tf)
			writeQualityScaledXStringSet(qss, tf)
			expect_equal(fastq.geometry(tf), c(n_seq, n_samp))

			if(stype == "DNA"){
				## only readQualityScaledDNAStringSet is currently defined
				qss2 <- readQualityScaledDNAStringSet(tf, quality.scoring=quality_type)

				expect_equal(as.character(qss), as.character(qss2))
				expect_equal(as.character(quality(qss)), as.character(quality(qss2)))
			}
			if(file.exists(tf)) file.remove(tf)
		}
	}
})

test_that("QualityScaledXStringSet properly implements reverse, reverseComplement", {
	.revString <- function(s) paste(strsplit(s, '')[[1]][seq(nchar(s), 1L)], collapse='')
	set.seed(10L)
	n_samp <- 40L
	n_seq <- 3L
	qs <- list(phred=PhredQuality,
							solexa=SolexaQuality,
							illumina=IlluminaQuality)

	alphs <- list(DNA=DNA_BASES,
								RNA=RNA_BASES,
								AA=AA_STANDARD,
								B=LETTERS)

	qualities <- lapply(qs, \(q){
			do.call(c, lapply(seq_len(n_seq), \(i){
					q(sample(40L, n_samp, replace=TRUE))
				}))
		})

	for(stype in c("DNA", "RNA", "AA", "B")){
		seqs <- vapply(seq_len(n_seq), \(x) paste(sample(alphs[[stype]], n_samp, replace=TRUE), collapse=''), character(1L))
		seqs <- get(paste0(stype, "StringSet"))(seqs)
		for(q in names(qs)){
			qual <- qualities[[q]]
			qss <- get(paste0("QualityScaled", stype, "StringSet"))(seqs, qual)
			exp_output <- vapply(as.character(qual), .revString, character(1L), USE.NAMES=FALSE)
			expect_equal(as.character(quality(reverse(qss))), exp_output)
			if(stype %in% c("DNA", "RNA")){
				expect_equal(as.character(quality(reverseComplement(qss))), exp_output)
			}
		}
	}
})