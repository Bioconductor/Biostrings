
test_that("element-wise comparison of XStringSet objects works as expected", {
    dna <- DNAStringSet(DNA_ALPHABET)
    rna <- RNAStringSet(RNA_ALPHABET)
    aaa <- AAStringSet(AA_ALPHABET)
    bbb <- BStringSet(LETTERS)

    expect_equal(match(dna, dna), seq_along(dna))
    expect_equal(aaa[seq_len(26)] < bbb, AA_ALPHABET[seq_len(26L)] < LETTERS)

    expect_equal(match(sort(aaa), bbb, nomatch=0), c(rep(0L, 4L), seq_len(26L)))
    expect_true(all(dna == as.character(dna)))

    expect_error2(aaa == dna, "is not supported")
    expect_error2(aaa == rna, "is not supported")
    expect_true(all(dna == BStringSet(DNA_ALPHABET)))
    expect_true(all(rna == BStringSet(RNA_ALPHABET)))
    expect_true(all(aaa == BStringSet(AA_ALPHABET)))
    expect_equal(dna == NULL, logical(0L))
})

