## TODO : add tests for append, subset

## WARNING: Do *NOT* use checkIdentical() on XString, XStringSet, or
## XStringSetList objects. It is *NOT* reliable. Use checkTrue(all.equal())
## instead.  

## TODO: Maybe "all.equal" should be made an S4 generic with S4/S3 method
## combos for XVector and XVectorList object?

## Unclass the DNAStringSet object and sets its "pool" and "ranges" slots
## to NULL.
.unclass_DNAStringSet <- function(x)
{
    x <- unclass(x)
    slot(x, "pool", check=FALSE) <- NULL
    slot(x, "ranges", check=FALSE) <- NULL
    x
}

all.equal.DNAStringSet <- function(target, current, ...)
{
    ## Compare the sequences (and their names if they have).
    target_seqs <- as.character(target)
    current_seqs <- as.character(current)
    ok1 <- all.equal(target_seqs, current_seqs)
    ## Compare the rest.
    target <- .unclass_DNAStringSet(target)
    current <- .unclass_DNAStringSet(current)
    ok2 <- all.equal(target, current)
    ans <- character(0)
    if (!isTRUE(ok1))
        ans <- c(ans, ok1)
    if (!isTRUE(ok2))
        ans <- c(ans, ok2)
    if (length(ans) == 0L)
        return(TRUE)
    ans
}

test_DNAStringSetList_constructor <- function()
{
    dna1 <- DNAStringSet(DNA_ALPHABET[1:8])
    dna2 <- DNAStringSet(DNA_ALPHABET[9:17])

    lst1 <- DNAStringSetList(dna1, dna2) 
    lst2 <- DNAStringSetList(as.character(DNA_ALPHABET[1:8]), 
                             as.character(DNA_ALPHABET[9:17])) 
    checkTrue(all.equal(lst1, lst2))

    checkTrue(length(DNAStringSetList()) == 0)
}

test_DNAStringSetList_unlist <- function()
{
    lst <- DNAStringSetList(DNA_ALPHABET, DNA_ALPHABET)
    expected <- DNAStringSet(c(DNA_ALPHABET, DNA_ALPHABET))
    checkTrue(all.equal(unlist(lst), expected))
}

test_DNAStringSetList_append <- function()
{
    dna <- DNAStringSet(DNA_ALPHABET)
    lst <- DNAStringSetList(DNA_ALPHABET, DNA_ALPHABET)
    elementMetadata(lst) <- DataFrame(C1=c("list1", "list2"))

    dna2a <- c(lst, lst)
    dna2b <- rep(lst, 2L)
    dna2c <- append(lst, lst) 
    checkTrue(all.equal(dna2a, dna2b))
    checkTrue(all.equal(dna2a, dna2c))
}

