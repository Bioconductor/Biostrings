## TODO : add tests for append, subset

test_DNAStringSetList_constructor <- function()
{
    dna1 <- DNAStringSet(DNA_ALPHABET[1:8])
    dna2 <- DNAStringSet(DNA_ALPHABET[9:17])

    lst1 <- DNAStringSetList(dna1, dna2) 
    lst2 <- DNAStringSetList(as.character(DNA_ALPHABET[1:8]), 
                             as.character(DNA_ALPHABET[9:17])) 
    checkIdentical(lst1, lst2)

    checkTrue(length(DNAStringSetList()) == 0)
}

test_DNAStringSetList_unlist <- function()
{
    lst <- DNAStringSetList(DNA_ALPHABET, DNA_ALPHABET)
    expected <- c(DNA_ALPHABET, DNA_ALPHABET)
    checkIdentical(as.character(unlist(lst)), expected)
}

test_DNAStringSetList_append <- function()
{
    dna <- DNAStringSet(DNA_ALPHABET)
    lst <- DNAStringSetList(DNA_ALPHABET, DNA_ALPHABET)
    elementMetadata(lst) <- DataFrame(C1=c("list1", "list2"))

    dna2a <- c(lst, lst)
    dna2b <- rep(lst, 2L)
    dna2c <- append(lst, lst) 
    checkIdentical(dna2a, dna2b)
    checkIdentical(dna2a, dna2c)
}

