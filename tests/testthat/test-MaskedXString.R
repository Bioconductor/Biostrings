## Testing for MaskedXString-class.R, maskMotif.R, injectHardMask.R
##
## MaskedXString-class.R exports:
## - MaskedXString class (B, AA, DNA, RNA)
## - unmasked
## - masks, masks<-
## - length
## - maskedratio
## - maskedwidth
## - nchar
## - seqtype
## - *Views
## - as.character
## - subseq
## - *toString
## - *show
## - *coercion between masked types
## - *coercion between unmasked types (B <-> MaskedB)
## - *coercion to MaskCollection, NormalIRanges, XStringViews

dna_c <- sample(DNA_BASES, 25L, replace=TRUE)
rna_c <- sample(RNA_BASES, 25L, replace=TRUE)
aaa_c <- sample(AA_STANDARD, 25L, replace=TRUE)
bbb_c <- sample(LETTERS, 25L, replace=TRUE)

dna <- DNAString(paste(dna_c, collapse=''))
rna <- RNAString(paste(rna_c, collapse=''))
aaa <- AAString(paste(aaa_c, collapse=''))
bbb <- BString(paste(bbb_c, collapse=''))

mask_starts <- c(1,5,10)
mask_ends <- c(3,7,12)
m1 <- Mask(25L, start=mask_starts, end=mask_ends)
m2 <- Mask(25L, start=c(15,20), end=c(18,22))
total_mask_width <- sum(mask_ends - mask_starts + 1L)

.check_masked_ascharacter <- function(maskedv, origv, maskstarts, maskends){
	c1 <- strsplit(as.character(maskedv), '')[[1]]
	c2 <- strsplit(as.character(origv), '')[[1]]
	if(sum(c1 == "#") != total_mask_width) return(FALSE)
	if(!all(c1[c1!='#'] == c2[c1!='#'])) return(FALSE)
	TRUE
}

test_that("masking coerces between XString and MaskedXString", {
	l <- list(dna, rna, aaa, bbb)
	for(i in seq_along(l)){
		tmp <- l[[i]]
		masks(tmp) <- m1
		## should convert to MaskedXString
		expect_s4_class(tmp, paste0("Masked", seqtype(l[[i]]), "String"))

		## seqtype should stay the same
		expect_equal(seqtype(tmp), seqtype(l[[i]]))

		## length should stay the same
		expect_equal(length(tmp), length(l[[i]]))

		## maskedwidth should report 9
		expect_equal(maskedwidth(tmp), total_mask_width)
		expect_equal(maskedratio(tmp), total_mask_width/nchar(l[[i]]))
		expect_equal(nchar(tmp), nchar(l[[i]]) - total_mask_width)

		## make sure that as.character preserves masking
		expect_true(.check_masked_ascharacter(tmp, l[[i]]))
		expect_equal(as.character(subseq(tmp, 4, 10)), substr(as.character(tmp), 4, 10))
		expect_equal(toString(tmp), as.character(tmp))

		## should convert back to XString
		expect_s4_class(as(tmp, paste0(seqtype(tmp), "String")), class(l[[i]]))

		## equality should work here because they're from the same pool (?)
		masks(tmp) <- m2
		expect_equal(masks(tmp), m2)
		expect_equal(unmasked(tmp), l[[i]])
	}
})