test_that("MIndexList creation and properties work", {
    patterns <- DNAStringSet(c("AAC", "TTA", "CGT"))
    pdict <- PDict(patterns)
    subjects <- DNAStringSet(c("AACTA", "TTAC", "CGTTAA", "NNN"))
    names(subjects) <- paste0("seq", 1:4)

    mindex_list <- vmatchPDict(pdict, subjects)
    expect_s4_class(mindex_list, "MIndexList")
    expect_equal(length(mindex_list), 4L)
    expect_equal(names(mindex_list), paste0("seq", 1:4))
})

test_that("MIndexList subsetting works", {
    patterns <- DNAStringSet(c("AAC", "TTA", "CGT"))
    pdict <- PDict(patterns)
    subjects <- DNAStringSet(c("AACTA", "TTAC", "CGTTAA", "NNN"))
    names(subjects) <- paste0("seq", 1:4)
    mindex_list <- vmatchPDict(pdict, subjects)

    # Double bracket subsetting returns MIndex
    m1 <- mindex_list[[1]]
    expect_s4_class(m1, "MIndex")

    # Single bracket subsetting returns MIndexList
    sub_list <- mindex_list[2:3]
    expect_s4_class(sub_list, "MIndexList")
    expect_equal(length(sub_list), 2L)
    expect_equal(names(sub_list), c("seq2", "seq3"))
})

test_that("MIndexList coercion works", {
    patterns <- DNAStringSet(c("AAC", "TTA", "CGT"))
    pdict <- PDict(patterns)
    subjects <- DNAStringSet(c("AACTA", "TTAC", "CGTTAA", "NNN"))
    names(subjects) <- paste0("seq", 1:4)
    mindex_list <- vmatchPDict(pdict, subjects)

    # Coercion to list
    lst <- as.list(mindex_list)
    expect_type(lst, "list")
    expect_equal(length(lst), 4L)
    expect_s4_class(lst[[1]], "MIndex")

    # Coercion from list to MIndexList
    mindex_list2 <- as(lst, "MIndexList")
    expect_s4_class(mindex_list2, "MIndexList")
    expect_equal(length(mindex_list2), 4L)
    expect_equal(names(mindex_list2), paste0("seq", 1:4))

    # Coercion from S4 List to MIndexList
    s4_list <- S4Vectors::List(lst[[1]], lst[[2]])
    mindex_list3 <- as(s4_list, "MIndexList")
    expect_s4_class(mindex_list3, "MIndexList")
    expect_equal(length(mindex_list3), 2L)

    # Invalid coercion triggers error
    expect_error(as(list(1, 2, 3), "MIndexList"),
                 "Only lists of MIndex objects can be coerced to MIndexList")
})

test_that("MIndexList show and showAsCell work", {
    patterns <- DNAStringSet(c("AAC", "TTA", "CGT"))
    pdict <- PDict(patterns)
    subjects <- DNAStringSet(c("AACTA", "TTAC", "CGTTAA", "NNN"))
    names(subjects) <- paste0("seq", 1:4)
    mindex_list <- vmatchPDict(pdict, subjects)

    # showAsCell
    cell_repr <- showAsCell(mindex_list)
    expect_equal(cell_repr, c("3 patterns, 1 matches",
                              "3 patterns, 1 matches",
                              "3 patterns, 2 matches",
                              "3 patterns, 0 matches"))

    # show
    expect_output(show(mindex_list), "MIndexList of length 4")
})
