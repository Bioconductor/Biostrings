
test_that("compatible_seqtypes()", {
    compatible_seqtypes <- Biostrings:::compatible_seqtypes
    SUPPORTED_SEQTYPES <- Biostrings:::.SUPPORTED_SEQTYPES
    target <- matrix(TRUE,
                     nrow=length(SUPPORTED_SEQTYPES),
                     ncol=length(SUPPORTED_SEQTYPES),
                     dimnames=list(SUPPORTED_SEQTYPES, SUPPORTED_SEQTYPES))
    target[c("DNA", "RNA"), "AA"] <- target["AA", c("DNA", "RNA")] <- FALSE
    current <- sapply(SUPPORTED_SEQTYPES,
        function(seqtype1) sapply(SUPPORTED_SEQTYPES,
            function(seqtype2) compatible_seqtypes(seqtype1, seqtype2))
    )
    expect_identical(current, target)
})

