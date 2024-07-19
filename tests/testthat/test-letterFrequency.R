## letterFrequency.R exports the following:
##	- alphabetFrequency
##	- hasOnlyBaseLetters
##	- uniqueLetters
##	- letterFrequency
##	- letterFrequencyInSlidingView
##	- consensusMatrix
##	- consensusString (matrix, DNAStringSet, RNAStringSet)

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
})

test_that("hasOnlyBaseLetters works as expected", {
	expect_true(hasOnlyBaseLetters(DNAStringSet(DNA_BASES)))
	expect_true(hasOnlyBaseLetters(RNAStringSet(RNA_BASES)))
	expect_true(hasOnlyBaseLetters(AAStringSet(AA_STANDARD)))

	expect_false(hasOnlyBaseLetters(DNAStringSet(DNA_ALPHABET)))
	expect_false(hasOnlyBaseLetters(RNAStringSet(RNA_ALPHABET)))
	expect_false(hasOnlyBaseLetters(AAStringSet(AA_ALPHABET)))

	expect_error(hasOnlyBaseLetters(BStringSet()), 'unable to find an inherited method')
})