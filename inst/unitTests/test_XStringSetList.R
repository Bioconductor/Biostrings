## TODO : add tests for append, subset

## WARNING: Do *NOT* use checkIdentical() on XString, XStringSet, or
## XStringSetList objects. It is *NOT* reliable. Use checkTrue(all.equal())
## instead.  

## TODO: Maybe "all.equal" should be made an S4 generic with S4/S3 method
## combos for XVector and XVectorList object?

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
    # .unclass_XStringSet() works interactively but fails when run in the
    # context of the unit tests.
    #target <- .unclass_XStringSet(target)
    #current <- .unclass_XStringSet(current)
    #ok2 <- all.equal(target, current)
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

.XStringSetList_constructor <- 
    function(XStringSetFUN, XStringSetListFUN, XS_ALPHABET)
{
    xs1 <- XStringSetFUN(XS_ALPHABET[1:8])
    xs2 <- XStringSetFUN(XS_ALPHABET[9:17])

    lst1 <- XStringSetListFUN(xs1, xs2) 
    lst2 <- XStringSetListFUN(as.character(XS_ALPHABET[1:8]), 
                              as.character(XS_ALPHABET[9:17])) 
    checkTrue(all.equal.XStringSetList(lst1, lst2))

    checkTrue(length(XStringSetListFUN()) == 0)
}

.XStringSetList_unlist <-
    function(XStringSetFUN, XStringSetListFUN, XS_ALPHABET)
{
    lst <- XStringSetListFUN(XS_ALPHABET, XS_ALPHABET)
    expected <- XStringSetFUN(c(XS_ALPHABET, XS_ALPHABET))
    checkTrue(all.equal.XStringSet(unlist(lst), expected))
}

.XStringSetList_append <-
    function(XStringSetFUN, XStringSetListFUN, XS_ALPHABET)
{
    xs <- XStringSetFUN(XS_ALPHABET)
    lst <- XStringSetListFUN(XS_ALPHABET, XS_ALPHABET)
    elementMetadata(lst) <- DataFrame(C1=c("list1", "list2"))

    xs2a <- c(lst, lst)
    xs2b <- rep(lst, 2L)
    xs2c <- append(lst, lst) 
    checkTrue(all.equal.XStringSetList(xs2a, xs2b))
    checkTrue(all.equal.XStringSetList(xs2a, xs2c))
}

## DNAStringSet
test_DNAStringSetList_constructor <- function()
    .XStringSetList_constructor(DNAStringSet, DNAStringSetList, DNA_ALPHABET)

test_DNAStringSetList_unlist <- function()
    .XStringSetList_unlist(DNAStringSet, DNAStringSetList, DNA_ALPHABET)

test_DNAStringSetList_append <- function()
    .XStringSetList_append(DNAStringSet, DNAStringSetList, DNA_ALPHABET)

## AAStringSet
test_AAStringSetList_constructor <- function()
    .XStringSetList_constructor(AAStringSet, AAStringSetList, AA_ALPHABET)

test_AAStringSetList_unlist <- function()
    .XStringSetList_unlist(AAStringSet, AAStringSetList, AA_ALPHABET)

test_AAStringSetList_append <- function()
    .XStringSetList_append(AAStringSet, AAStringSetList, AA_ALPHABET)


test_DNAStringSetList_showAsCell <- function()
{
    dna <- showAsCell(DNAStringSetList())
    checkTrue(is(dna, "character"))
    dna <- showAsCell(DNAStringSetList(DNA_ALPHABET))
    checkTrue(is(dna, "character"))
}
