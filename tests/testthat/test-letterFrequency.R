## letterFrequency.R exports the following:
##	- alphabetFrequency
##	- hasOnlyBaseLetters
##	- uniqueLetters
##	- letterFrequency
##	- *letterFrequencyInSlidingView
##	- consensusMatrix
##	- consensusString (matrix, DNAStringSet, RNAStringSet)
##  - mkAllStrings
##  - oligonucleotideFrequency
##	- dinucleotideFrequency
##	- trinucleotideFrequency
## 	- *oligonucleotideTransitions
##  - *twoWayAlphabetFrequency
##
## starred functions are still lacking tests

B_ALPHABET <- strsplit(rawToChar(as.raw(32:126)), '')[[1]]

check_uniqueLetters <- function(XStringFUN, X_ALPH){
	set.seed(123L)
	s <- sample(X_ALPH, 100L, replace=TRUE)
	u_v <- unique(s)

	s <- paste(s, collapse='')
	expect_equal(sort(uniqueLetters(XStringFUN(s))), sort(u_v))
}

check_alphfreq <- function(XStringFUN, X_ALPH, isMatrix=FALSE, checkBaseOnly=NULL){
	set.seed(456L)
	tab <- integer(length(X_ALPH))
	names(tab) <- X_ALPH

	if(isMatrix){
		exp_res <- matrix(0L, nrow=5, ncol=length(tab))
		colnames(exp_res) <- names(tab)
		strs <- character(5L)
		for(i in seq_len(5)){
			tab[] <- 0L
			s <- sample(X_ALPH, 100L, replace=TRUE)
			u_v <- table(s)
			tab[names(u_v)] <- u_v
			exp_res[i,] <- tab
			strs[i] <- paste(s, collapse='')
		}
		norm_res <- exp_res / rowSums(exp_res)
		ss <- XStringFUN(strs)
	} else {
		strs <- sample(X_ALPH, 100L, replace=TRUE)
		u_v <- table(strs)
		tab[names(u_v)] <- u_v
		exp_res <- tab
		ss <- XStringFUN(paste(strs, collapse=''))
		norm_res <- exp_res / sum(exp_res)
	}

	# Check alphabetFrequency
	v <- alphabetFrequency(ss)
	expect_equal(v, exp_res)
	expect_equal(sum(v), sum(nchar(strs)))
	expect_equal(alphabetFrequency(ss, as.prob=TRUE), norm_res)

	# Check letterFrequency
	expect_equal(letterFrequency(ss, X_ALPH), exp_res)
	expect_equal(letterFrequency(ss, X_ALPH, as.prob=TRUE), norm_res)

	if(!is.null(checkBaseOnly)){
		# provide the alphabet
		b <- checkBaseOnly
		if(isMatrix){
			p_in <- match(colnames(exp_res), b, nomatch=0)
			exp_res <- cbind(exp_res[,p_in!=0], other=rowSums(exp_res[,p_in==0]))
			norm_res <- exp_res / rowSums(exp_res)
		} else {
			p_in <- match(names(exp_res), b, nomatch=0)
			exp_res <- c(exp_res[p_in!=0], sum(exp_res[p_in==0]))
			names(exp_res) <- c(b, "other")
			norm_res <- exp_res / sum(exp_res)
		}

		v <- alphabetFrequency(ss, baseOnly=TRUE)
		expect_equal(v, exp_res)
		expect_equal(sum(v), sum(nchar(strs)))
		expect_error(alphabetFrequency(ss, baseOnly=NA), "'baseOnly' must be TRUE or FALSE")
		expect_error(alphabetFrequency(ss, baseOnly=1), "'baseOnly' must be TRUE or FALSE")
		expect_equal(alphabetFrequency(ss, baseOnly=TRUE, as.prob=TRUE), norm_res)
		if(isMatrix){
			expect_equal(letterFrequency(ss, b), exp_res[,seq_along(b)])
			expect_equal(letterFrequency(ss, b, as.prob=TRUE), norm_res[,seq_along(b)])
		} else {
			expect_equal(letterFrequency(ss, b), exp_res[seq_along(b)])
			expect_equal(letterFrequency(ss, b, as.prob=TRUE), norm_res[seq_along(b)])
		}
	}
}

test_that("uniqueLetters works properly", {
	check_uniqueLetters(DNAString, DNA_ALPHABET)
	check_uniqueLetters(RNAString, RNA_ALPHABET)
	check_uniqueLetters(AAString, AA_ALPHABET)
	check_uniqueLetters(BString, B_ALPHABET)

	check_uniqueLetters(DNAStringSet, DNA_ALPHABET)
	check_uniqueLetters(RNAStringSet, RNA_ALPHABET)
	check_uniqueLetters(AAStringSet, AA_ALPHABET)
	check_uniqueLetters(BStringSet, B_ALPHABET)
})

test_that("letterFrequency, alphabetFrequency functions for all input types", {
	check_alphfreq(DNAString, DNA_ALPHABET, checkBaseOnly=DNA_BASES)
	check_alphfreq(RNAString, RNA_ALPHABET, checkBaseOnly=RNA_BASES)
	check_alphfreq(AAString, AA_ALPHABET, checkBaseOnly=AA_STANDARD)

	check_alphfreq(DNAStringSet, DNA_ALPHABET, isMatrix=TRUE, checkBaseOnly=DNA_BASES)
	check_alphfreq(RNAStringSet, RNA_ALPHABET, isMatrix=TRUE, checkBaseOnly=RNA_BASES)
	check_alphfreq(AAStringSet, AA_ALPHABET, isMatrix=TRUE, checkBaseOnly=AA_STANDARD)

	## checking for ambiguity tabulation
	d <- DNAString("ATGCATGC.+Y")
	exp_output <- c(2L, 4L, 2L)
	names(exp_output) <- c("A", "G|C", "T")
	expect_equal(letterFrequency(d, letters=c("A", "GC", "T")), exp_output)
	names(exp_output)[2] <- "G*C"
	expect_equal(letterFrequency(d, letters=c("A", "GC", "T"), OR="*"), exp_output)

	r <- RNAString("AUGCAUGC.+Y")
	exp_output <- 8L
	names(exp_output) <- c("A|U|G|C")
	expect_equal(letterFrequency(r, letters="AUGC"), exp_output)
	exp_output <- rep(2L, 4)
	names(exp_output) <- c("A","U","G","C") # intentionally different from RNA_BASES
	expect_equal(letterFrequency(r, letters="AUGC", OR=0), exp_output)

	## BStrings have their own check
	set.seed(999L)
	bstr <- sample(B_ALPHABET, 100, replace=TRUE)
	bstr_counts <- table(match(bstr, B_ALPHABET))
	r <- integer(256L)
	r[as.integer(names(bstr_counts)) + 32L] <- bstr_counts
	bstr <- paste(bstr, collapse='')
	expect_equal(alphabetFrequency(BString(bstr)), r)
	expect_equal(alphabetFrequency(BStringSet(bstr)), matrix(r, nrow=1L, ncol=256L))

	expect_error(alphabetFrequency(BString(), baseOnly=TRUE), "unused argument (baseOnly = TRUE)", fixed=TRUE)

	## checking for Views objects and maskedString objects
	ds <- sample(DNA_ALPHABET, 100L, replace=TRUE)
	exp_out <- matrix(0L, nrow=2, ncol=length(DNA_ALPHABET))
	colnames(exp_out) <- DNA_ALPHABET
	exp_out[1,names(table(ds[seq(1,50)]))] <- table(ds[seq(1,50)])
	exp_out[2,names(table(ds[seq(51,100)]))] <- table(ds[seq(51,100)])
	v <- Views(DNAString(paste(ds, collapse='')), start=c(1,51), end=c(50,100))
	expect_equal(alphabetFrequency(v), exp_out)
	expect_equal(letterFrequency(v, "A"), exp_out[,"A",drop=FALSE])

	mask_d <- DNAString(paste(ds, collapse=''))
	masks(mask_d) <- Mask(length(ds), start=c(5,50), end=c(10,60))
	exp_out <- rep(0L, length(DNA_ALPHABET))
	names(exp_out) <- DNA_ALPHABET
	subdstr <- ds[-c(5:10,50:60)]
	exp_out[names(table(subdstr))] <- table(subdstr)
	expect_equal(alphabetFrequency(mask_d), exp_out)
	expect_equal(letterFrequency(mask_d, "A"), exp_out["A"])
})

test_that("hasOnlyBaseLetters works as expected", {
	## XStringSet
	expect_true(hasOnlyBaseLetters(DNAStringSet(DNA_BASES)))
	expect_true(hasOnlyBaseLetters(RNAStringSet(RNA_BASES)))
	expect_true(hasOnlyBaseLetters(AAStringSet(AA_STANDARD)))

	expect_false(hasOnlyBaseLetters(DNAStringSet(DNA_ALPHABET)))
	expect_false(hasOnlyBaseLetters(RNAStringSet(RNA_ALPHABET)))
	expect_false(hasOnlyBaseLetters(AAStringSet(AA_ALPHABET)))
	expect_error(hasOnlyBaseLetters(BStringSet()), 'unable to find an inherited method')

	## Views
	d1 <- DNAString(paste(DNA_BASES, collapse=''))
	d2 <- DNAString(paste(DNA_ALPHABET, collapse=''))
	expect_true(hasOnlyBaseLetters(Views(d1, start=1, end=length(DNA_BASES))))
	expect_false(hasOnlyBaseLetters(Views(d2, start=1, end=length(DNA_ALPHABET))))

	## MaskedXString
	md1 <- d1
	md2 <- d2
	masks(md1) <- Mask(length(d1), start=2, end=3)
	masks(md2) <- Mask(length(d2), start=2, end=5)
	expect_true(hasOnlyBaseLetters(md1))
	expect_false(hasOnlyBaseLetters(md2))

})

test_that("miscellaneous letterFrequency methods work", {
	## mkAllStrings
	alph <- c("A", "B", "C")
	exp_output <- expand.grid(alph, alph, alph)
	exp_output <- apply(exp_output, 1L, paste, collapse='')
	expect_equal(sort(mkAllStrings(alph, 3L)), sort(exp_output))

	## uniqueLetters
	d <- paste(rep(DNA_ALPHABET, each=3), collapse='')
	expect_equal(uniqueLetters(DNAString(d)), DNA_ALPHABET)
	expect_equal(uniqueLetters(DNAStringSet(c(d,d))), DNA_ALPHABET)
	expect_equal(uniqueLetters(Views(DNAString(d), start=c(1,6), end=c(5,nchar(d)))), DNA_ALPHABET)
	md <- DNAString(d)
	masks(md) <- Mask(nchar(d), start=seq(1,nchar(d),3), end=seq(2,nchar(d),3))
	expect_equal(uniqueLetters(md), DNA_ALPHABET)
})

test_that("nucleotideFrequency methods work properly", {
	set.seed(999L)
	l <- 102L
	ds <- sample(DNA_BASES, l, replace=TRUE)
	d <- DNAString(paste(ds, collapse=''))
	dss <- DNAStringSet(c(ds,ds))

	## oligonucleotideFrequency(x, 1L) should be the same as alphabetFrequency(...,baseOnly=TRUE)
	## (ignoring the 'other' arg)
	expect_equal(oligonucleotideFrequency(d, 1L), alphabetFrequency(d, baseOnly=TRUE)[seq_len(4L)])
	expect_equal(oligonucleotideFrequency(dss, 1L), alphabetFrequency(dss, baseOnly=TRUE)[,seq_len(4L)])

	## dinucleotide correctness
	all_doubles <- table(paste0(ds[seq_len(l-1L)], ds[seq_len(l-1L)+1L]))
	exp_res <- rep(0L, length(DNA_BASES)**2)
	names(exp_res) <- mkAllStrings(DNA_BASES, 2L)
	exp_res[names(all_doubles)] <- all_doubles
	expect_equal(oligonucleotideFrequency(d, 2L), exp_res)
	expect_equal(dinucleotideFrequency(d), exp_res)

	## trinucleotide correctness
	all_triples <- table(paste0(ds[seq_len(l-2L)], ds[seq_len(l-2L)+1L], ds[seq_len(l-2L)+2L]))
	exp_res <- rep(0L, length(DNA_BASES)**3)
	names(exp_res) <- mkAllStrings(DNA_BASES, 3L)
	exp_res[names(all_triples)] <- all_triples
	expect_equal(oligonucleotideFrequency(d, 3L), exp_res)
	expect_equal(trinucleotideFrequency(d), exp_res)

	## using step
	all_codons <- codons(d)
	exp_res[] <- 0L
	codon_table <- table(as.character(all_codons))
	exp_res[names(codon_table)] <- codon_table
	expect_equal(oligonucleotideFrequency(d, 3L, step=3L), exp_res)

	## frequency on a Views object
	exp_res <- matrix(0L, nrow=l/3, ncol=length(DNA_BASES)**3)
	colnames(exp_res) <- mkAllStrings(DNA_BASES, 3L)
	exp_res[cbind(seq_len(nrow(exp_res)), match(as.character(all_codons), colnames(exp_res)))] <- 1L
	expect_equal(oligonucleotideFrequency(codons(d), 3L), exp_res)

	## frequency on a MaskedXString object
	mask_d <- d
	masks(mask_d) <- Mask(l, start=c(10,50), end=c(19,57))
	mask_str <- strsplit(as.character(mask_d), '')[[1]]
	all_doubles <- table(paste0(mask_str[seq_len(l-1L)], mask_str[seq_len(l-1L)+1L]))
	all_doubles <- all_doubles[!grepl("#", names(all_doubles))]
	exp_res <- rep(0L, length(DNA_BASES)**2)
	names(exp_res) <- mkAllStrings(DNA_BASES, 2L)
	exp_res[names(all_doubles)] <- all_doubles
	expect_equal(oligonucleotideFrequency(mask_d, 2L), exp_res)

	## disallowed types
	expect_error(oligonucleotideFrequency(AAString(), "must contain sequences of type DNA or RNA"))
	expect_error(oligonucleotideFrequency(BStringSet(), "must contain sequences of type DNA or RNA"))
})