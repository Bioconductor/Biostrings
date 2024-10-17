## WARNING: Do *NOT* use checkIdentical() on XString, XStringSet, or
## XStringSetList objects. It is *NOT* reliable. Use checkTrue(all.equal())
## instead.

## TODO: Maybe "all.equal" should be made an S4 generic with S4/S3 method
## combos for XVector and XVectorList object?

B_ALPHABET <- strsplit(rawToChar(as.raw(32:126)), '')[[1]]

## Unclass the XStringSet object by setting its "pool" and "ranges" slots
## to NULL first.
.unclass_XStringSet <- function(x)
{
    slot(x, "pool", check=FALSE) <- NULL
    slot(x, "ranges", check=FALSE) <- NULL
    unclass(x)
}

all.equal.XStringSet <- function(target, current, ...)
{
    ## Compare the sequences (and their names if they have).
    target_seqs <- as.character(target)
    current_seqs <- as.character(current)
    ok1 <- all.equal(target_seqs, current_seqs)
    ## Compare the rest.

    ## this was previously commented out, but works with testthat
    target <- .unclass_XStringSet(target)
    current <- .unclass_XStringSet(current)
    ok2 <- all.equal(target, current)
    ##

    ok2 <- identical(class(target), class(current))
    if (!ok2)
        ok2 <- "class mismatch"
    ok3 <- all.equal(target@metadata, current@metadata)
    ok4 <- all.equal(target@elementMetadata, current@elementMetadata)
    ans <- character(0)
    if (!isTRUE(ok1))
        ans <- c(ans, ok1)
    if (!isTRUE(ok2))
        ans <- c(ans, ok2)
    if (!isTRUE(ok3))
        ans <- c(ans, ok3)
    if (!isTRUE(ok4))
        ans <- c(ans, ok4)
    if (length(ans) == 0L)
        return(TRUE)
    ans
}

all.equal.XStringSetList <- function(target, current, ...)
{
    ok1 <- identical(class(target), class(current))
    if (!ok1)
        ok1 <- "class mismatch"
    target2 <- as(target, "CharacterList")
    current2 <- as(current, "CharacterList")
    ## Temporary workaround until coercion from XStringSet to CharacterList
    ## is fixed to propagate metadata and metadata columns.
    metadata(target2) <- metadata(target)
    mcols(target2) <- mcols(target)
    metadata(current2) <- metadata(current)
    mcols(current2) <- mcols(current)
    ok2 <- all.equal(target2, current2)
    ans <- character(0)
    if (!isTRUE(ok1))
        ans <- c(ans, ok1)
    if (!isTRUE(ok2))
        ans <- c(ans, ok2)
    if (length(ans) == 0L)
        return(TRUE)
    ans
}

.XStringSetList_constructor <- function(XStringSetFUN, XStringSetListFUN, XS_ALPHABET){
    xs1 <- XStringSetFUN(XS_ALPHABET[1:8])
    xs2 <- XStringSetFUN(XS_ALPHABET[9:17])

    lst1 <- XStringSetListFUN(xs1, xs2)
    lst2 <- XStringSetListFUN(as.character(XS_ALPHABET[1:8]),
                              as.character(XS_ALPHABET[9:17]))

    expect_true(all.equal.XStringSetList(lst1, lst2))
    expect_equal(length(XStringSetListFUN()), 0L)
}

.XStringSetList_unlist <- function(XStringSetFUN, XStringSetListFUN, XS_ALPHABET){
    lst <- XStringSetListFUN(XS_ALPHABET, XS_ALPHABET)
    expected <- XStringSetFUN(c(XS_ALPHABET, XS_ALPHABET))
    expect_true(all.equal.XStringSet(unlist(lst), expected))
}

.XStringSetList_append <- function(XStringSetFUN, XStringSetListFUN, XS_ALPHABET){
    xs <- XStringSetFUN(XS_ALPHABET)
    lst <- XStringSetListFUN(XS_ALPHABET, XS_ALPHABET)
    elementMetadata(lst) <- DataFrame(C1=c("list1", "list2"))

    xs2a <- c(lst, lst)
    xs2b <- rep(lst, 2L)
    xs2c <- append(lst, lst)
    expect_true(all.equal.XStringSetList(xs2a, xs2b))
    expect_true(all.equal.XStringSetList(xs2a, xs2c))
}

.XStringSetList_subset <- function(XStringSetFUN, XStringSetListFUN, XS_ALPHABET){
    lst1 <- XStringSetListFUN(XS_ALPHABET[-1], XS_ALPHABET[-2], XS_ALPHABET[-3], XS_ALPHABET[-4])
    lst2 <- XStringSetListFUN(XS_ALPHABET[-2], XS_ALPHABET[-1], XS_ALPHABET[-4])

    expect_true(all.equal.XStringSetList(lst1[1:2], lst2[2:1]))
    expect_true(all.equal.XStringSetList(lst1[-c(2:3)], lst2[-1]))
}

.XStringSetList_from_list <- function(XStringSetFUN, XStringSetListFUN, XS_ALPHABET){
    lst1 <- list(XStringSetFUN(XS_ALPHABET), XStringSetFUN(XS_ALPHABET))
    lst2 <- XStringSetListFUN(XS_ALPHABET, XS_ALPHABET)

    expect_true(all.equal.XStringSetList(XStringSetListFUN(lst1), lst2))
}

## DNAStringSet
test_that("XStringSetList constructor works correctly", {
    .XStringSetList_constructor(DNAStringSet, DNAStringSetList, DNA_ALPHABET)
    .XStringSetList_constructor(RNAStringSet, RNAStringSetList, RNA_ALPHABET)
    .XStringSetList_constructor(AAStringSet, AAStringSetList, AA_ALPHABET)
    .XStringSetList_constructor(BStringSet, BStringSetList, B_ALPHABET)
})

test_that("XStringSetList unlist works correctly", {
    .XStringSetList_unlist(DNAStringSet, DNAStringSetList, DNA_ALPHABET)
    .XStringSetList_unlist(RNAStringSet, RNAStringSetList, RNA_ALPHABET)
    .XStringSetList_unlist(AAStringSet, AAStringSetList, AA_ALPHABET)
    .XStringSetList_unlist(BStringSet, BStringSetList, B_ALPHABET)
})

test_that("XStringSetList append works correctly", {
    .XStringSetList_append(DNAStringSet, DNAStringSetList, DNA_ALPHABET)
    .XStringSetList_append(RNAStringSet, RNAStringSetList, RNA_ALPHABET)
    .XStringSetList_append(AAStringSet, AAStringSetList, AA_ALPHABET)
    .XStringSetList_append(BStringSet, BStringSetList, B_ALPHABET)
})

test_that("XStringSetList subset works correctly", {
    .XStringSetList_subset(DNAStringSet, DNAStringSetList, DNA_ALPHABET)
    .XStringSetList_subset(RNAStringSet, RNAStringSetList, RNA_ALPHABET)
    .XStringSetList_subset(AAStringSet, AAStringSetList, AA_ALPHABET)
    .XStringSetList_subset(BStringSet, BStringSetList, B_ALPHABET)
})

test_that("XStringSetList showAsCell works correctly", {
    expect_equal(showAsCell(DNAStringSetList()), character(0L))
    expect_equal(showAsCell(DNAStringSetList(DNA_ALPHABET)), "A,C,G,...")
})

test_that("XStringSetList seqtype is correct", {
    expect_equal(seqtype(DNAStringSetList()), "DNA")
    expect_equal(seqtype(RNAStringSetList()), "RNA")
    expect_equal(seqtype(AAStringSetList()), "AA")
    expect_equal(seqtype(BStringSetList()), "B")

    # downgrade sequences
    x <- BStringSetList(DNA_ALPHABET)
    seqtype(x) <- "DNA"
    all.equal.XStringSetList(x, DNAStringSetList(DNA_ALPHABET))
    seqtype(x) <- "RNA"
    all.equal.XStringSetList(x, RNAStringSetList(RNA_ALPHABET))

    x <- BStringSetList(AA_ALPHABET)
    seqtype(x) <- "AA"
    all.equal.XStringSetList(x, AAStringSetList(AA_ALPHABET))

    # invalid conversions
    x <- DNAStringSetList(DNA_ALPHABET)
    expect_error({seqtype(x) <- "AA"}, 'incompatible sequence types "DNA" and "AA"')
    seqtype(x) <- "RNA"
    expect_error({seqtype(x) <- "AA"}, 'incompatible sequence types "RNA" and "AA"')
})

test_that("XStringSetLists show correctly", {
    ## TODO: add some values in case printing only breaks on output
    expect_output(show(DNAStringSetList(DNA_ALPHABET)), "DNAStringSetList of length 1\\n\\[\\[1\\]\\]")
    expect_output(show(RNAStringSetList(RNA_ALPHABET)), "RNAStringSetList of length 1\\n\\[\\[1\\]\\]")
    expect_output(show(AAStringSetList(AA_ALPHABET)), "AAStringSetList of length 1\\n\\[\\[1\\]\\]")
    expect_output(show(BStringSetList(B_ALPHABET)), "BStringSetList of length 1\\n\\[\\[1\\]\\]")
})

test_that("XStringSetLists coerce from list correctly", {
    .XStringSetList_from_list(DNAStringSet, DNAStringSetList, DNA_ALPHABET)
    .XStringSetList_from_list(RNAStringSet, RNAStringSetList, RNA_ALPHABET)
    .XStringSetList_from_list(AAStringSet, AAStringSetList, AA_ALPHABET)
    .XStringSetList_from_list(BStringSet, BStringSetList, B_ALPHABET)
})

test_that("XStringSetList uses nchar correctly", {
    dna <- DNAStringSetList(DNA_ALPHABET)
    expect_s4_class(nchar(dna), "IntegerList")
    expect_equal(unlist(nchar(dna)), rep(1L, length(DNA_ALPHABET)))
})