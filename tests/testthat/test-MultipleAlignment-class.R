## MultipleAlignment.R exports the following:
## - XMultipleAlignment constructors (DNA, RNA, AA)
## - readXMultipleAlignment functions (DNA, RNA, AA)
## - write.phylip

test_that("write.phylip only functions for MultipleAlignment objects", {
    tf <- tempfile()
    dnastr <- "ATGC"
    expect_null(write.phylip(DNAMultipleAlignment(dnastr), tf))
    expect_equal(unname(as.character(readDNAMultipleAlignment(tf, format="phylip"))),
        dnastr)

    expect_error(write.phylip(DNAStringSet(dnastr), tf),
                    "must be a MultipleAlignment object or derivative")
})
