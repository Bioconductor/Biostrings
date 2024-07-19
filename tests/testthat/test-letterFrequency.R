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
		strs <- character(1L)
		for(i in seq_len(5)){
			tab[] <- 0L
			s <- sample(X_ALPH, 100L, replace=TRUE)
			u_v <- table(s)
			tab[names(u_v)] <- u_v
			exp_res[i,] <- tab
			strs[i] <- paste(s, collapse='')
		}
		ss <- XStringFUN(strs)
	} else {
		s <- sample(X_ALPH, 100L, replace=TRUE)
		u_v <- table(s)
		tab[names(u_v)] <- u_v
		exp_res <- tab
		ss <- XStringFUN(paste(s, collapse=''))
	}

	expect_equal(alphabetFrequency(ss), exp_res)

	if(!is.null(checkBaseOnly)){
		# provide the alphabet
		b <- checkBaseOnly
		if(isMatrix){
			p_in <- match(colnames(exp_res), b, nomatch=0)
			exp_res <- cbind(exp_res[,p_in!=0], other=rowSums(exp_res[,p_in==0]))
		} else {
			p_in <- match(names(exp_res), b, nomatch=0)
			exp_res <- c(exp_res[p_in!=0], sum(exp_res[p_in==0]))
			names(exp_res) <- c(b, "other")
		}

		expect_equal(alphabetFrequency(ss, baseOnly=TRUE), exp_res)
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

test_that("alphabetFrequency functions for all input types", {
	check_alphfreq(DNAString, DNA_ALPHABET, checkBaseOnly=DNA_BASES)
	check_alphfreq(RNAString, RNA_ALPHABET, checkBaseOnly=RNA_BASES)
	check_alphfreq(AAString, AA_ALPHABET, checkBaseOnly=AA_STANDARD)

	check_alphfreq(DNAStringSet, DNA_ALPHABET, isMatrix=TRUE, checkBaseOnly=DNA_BASES)
	check_alphfreq(RNAStringSet, RNA_ALPHABET, isMatrix=TRUE, checkBaseOnly=RNA_BASES)
	check_alphfreq(AAStringSet, AA_ALPHABET, isMatrix=TRUE, checkBaseOnly=AA_STANDARD)

	## BStrings have their own check
	set.seed(999L)
	bstr <- sample(B_ALPHABET, 100, replace=TRUE)
	bstr_counts <- table(match(bstr, B_ALPHABET))
	r <- integer(256L)
	r[as.integer(names(bstr_counts)) + 32L] <- bstr_counts
	bstr <- paste(bstr, collapse='')
	expect_equal(alphabetFrequency(BString(bstr)), r)
	expect_equal(alphabetFrequency(BStringSet(bstr)), matrix(r, nrow=1L, ncol=256L))
})