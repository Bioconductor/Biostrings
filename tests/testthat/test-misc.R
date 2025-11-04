
test_that("N50() works as expected", {
    set.seed(6L)
    n_seq <- 20L
    possible_len <- seq(100,10000)
    lens <- sample(possible_len, n_seq)
    many_seqs <- DNAStringSet(
        vapply(seq_len(n_seq), function(i)
               paste(sample(DNA_BASES, lens[i], replace=TRUE), collapse=""),
        character(1L)))

    ## N50 is calculated by:
    ##  1. add sizes of contigs until half the total size is reached
    all_widths <- sort(width(many_seqs), decreasing=TRUE)
    total_width <- sum(all_widths)
    cs <- cumsum(all_widths)
    pos_first <- which.max(cs > total_width/2)

    ##  2. take the size of the contig that was added last
    N50_exp <- all_widths[pos_first]

    expect_equal(N50(width(many_seqs)), N50_exp)
})

