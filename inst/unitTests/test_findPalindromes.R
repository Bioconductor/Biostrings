###

test_findComplementedPalindromes <- function()
{
    seq1 <- DNAString("AATTT")
    current <- findComplementedPalindromes(seq1, min.armlength=2)
    checkIdentical(IRanges(1, 4:5), ranges(current))
    current <- findComplementedPalindromes(seq1, min.armlength=2,
                                                 max.looplength=0)
    checkIdentical(IRanges(1, 4), ranges(current))

    seq2 <- DNAString("TTTAA")
    current <- findComplementedPalindromes(seq2, min.armlength=2)
    checkIdentical(IRanges(1:2, 5), ranges(current))
    current <- findComplementedPalindromes(seq2, min.armlength=2,
                                                 max.looplength=0)
    checkIdentical(IRanges(2, 5), ranges(current))

    x <- DNAString("AAATTT")
    current <- findComplementedPalindromes(x, min.armlength=2)
    checkIdentical(IRanges(c(1, 1, 2), c(5, 6, 6)), ranges(current))
    current <- findComplementedPalindromes(x, min.armlength=3)
    checkIdentical(IRanges(1, 6), ranges(current))
    current <- findComplementedPalindromes(x, min.armlength=4)
    checkIdentical(IRanges(), ranges(current))

    ## With nested palindromes.

    x <- DNAString("TTGAATTGAA")
    current <- findComplementedPalindromes(x, min.armlength=2, max.looplength=0)
    checkIdentical(IRanges(4, 7), ranges(current))
    for (i in 1:5) {
        current <- findComplementedPalindromes(x, min.armlength=2,
                                                  max.looplength=i)
        checkIdentical(IRanges(c(1, 4, 6), c(5, 7, 10)), ranges(current))
    }
    current <- findComplementedPalindromes(x, min.armlength=2, max.looplength=6)
    checkIdentical(IRanges(c(1, 4, 1, 6), c(5, 7, 10, 10)), ranges(current))

    x <- DNAString("AAGAATTGTT")
    for (i in 0:2) {
        current <- findComplementedPalindromes(x, min.armlength=2,
                                                  max.looplength=i)
        checkIdentical(IRanges(4, 7), ranges(current))
    }
    for (i in 3:5) {
        current <- findComplementedPalindromes(x, min.armlength=2,
                                                  max.looplength=i)
        checkIdentical(IRanges(c(1, 4, 4), c(7, 7, 10)), ranges(current))
    }
    current <- findComplementedPalindromes(x, min.armlength=2, max.looplength=6)
    checkIdentical(IRanges(c(1, 4, 1, 4), c(7, 7, 10, 10)), ranges(current))

    # With IUPAC ambiguity codes. Ambiguity codes are not allowed in the arms
    # of a palindrome!

    x <- DNAString("NNNNNN")
    current <- findComplementedPalindromes(x, min.armlength=2)
    checkIdentical(IRanges(), ranges(current))
    
    x <- DNAString("ACCGNAAATTTNCGGT")
    current <- findComplementedPalindromes(x, min.armlength=2)
    checkIdentical(IRanges(c(6, 6, 7), c(10, 11, 11)), ranges(current))
    current <- findComplementedPalindromes(x, min.armlength=4)
    checkIdentical(IRanges(), ranges(current))
    current <- findComplementedPalindromes(x, min.armlength=3, max.looplength=8)
    checkIdentical(IRanges(c(6, 1), c(11, 16)), ranges(current))
}

