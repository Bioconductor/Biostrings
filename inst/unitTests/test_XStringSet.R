### -------------------------------------------------------------------------
### Helper functions
###

### In R 2.14 (and maybe before that), 2 external pointers are always
### considered identical so the identical() function cannot be used to
### compare the "pool" slots of 2 XVectorList objects. The workaround we
### use below is to extract the adresses in each pool as a character vector,
### and then to compare the 2 character vectors.

### '.eltAddresses(x)' collects the addresses of the elements in 'x' (in
### practice 'x' will be a list of external pointers or environments).
.eltAddresses <- function(x) sapply(x, XVector:::address)

### 'x' and 'y' must be XVectorList vectors.
.haveIdenticalPools <- function(x, y)
    identical(.eltAddresses(x@pool@xp_list), .eltAddresses(y@pool@xp_list))

### 'x' must be an XVectorList vector.
.poolEltLengths <- function(x)
{
    pool_len <- length(x@pool)
    if (pool_len == 0L)
        return(integer(0))
    sapply(seq_len(pool_len), function(i) length(x@pool[[i]]))
}


### -------------------------------------------------------------------------

test_width_character <- function()
{
    x <- safeExplode(rawToChar(as.raw(1:255)))
    checkIdentical(width(x), rep.int(1L, 255))
}

test_DNAStringSet_constructor <- function()
{
    dna <- DNAStringSet(DNA_ALPHABET)

    ## Checking internal representation.
    checkIdentical(.poolEltLengths(dna), length(DNA_ALPHABET))
}

test_DNAStringSet_width <- function()
{
    dna <- DNAStringSet(DNA_ALPHABET)
    checkIdentical(width(dna), width(DNA_ALPHABET))
}

test_DNAStringSet_subsetting <- function()
{
    dna <- DNAStringSet(DNA_ALPHABET)
    elementMetadata(dna) <- DataFrame(C1=dna)

    dna0 <- dna[FALSE]
    checkIdentical(length(dna0), 0L)
    ## Checking internal representation.
    checkIdentical(.poolEltLengths(dna0), integer(0))
    checkIdentical(.haveIdenticalPools(elementMetadata(dna0)$C1, dna0),
                   TRUE)
    checkIdentical(elementMetadata(dna0)$C1@ranges, dna0@ranges)

    idx <- rep.int((8:6)*2L, 100L)
    dna300 <- dna[idx]
    checkIdentical(length(dna300), length(idx))
    ## Checking internal representation.
    checkIdentical(.haveIdenticalPools(dna300, dna), TRUE)
    checkIdentical(.haveIdenticalPools(elementMetadata(dna300)$C1, dna300),
                   TRUE)
    checkIdentical(elementMetadata(dna300)$C1@ranges, dna300@ranges)
}

test_DNAStringSet_combining <- function()
{
    dna <- DNAStringSet(DNA_ALPHABET)
    elementMetadata(dna) <- DataFrame(C1=dna)

    dna2a <- c(dna, dna)
    dna2b <- rep(dna, 2L)
    checkIdentical(dna2a, dna2b)
    ## Checking internal representation.
    checkIdentical(.haveIdenticalPools(dna2a, dna), TRUE)
    checkIdentical(.haveIdenticalPools(dna2a, dna2b), TRUE)
    checkIdentical(.haveIdenticalPools(elementMetadata(dna2a)$C1, dna2a),
                   TRUE)
    checkIdentical(elementMetadata(dna2a)$C1@ranges, dna2a@ranges)
}

test_DNAStringSet_unlist <- function()
{
    dna <- DNAStringSet(DNA_ALPHABET)
    checkIdentical(as.character(unlist(dna)), paste(DNA_ALPHABET, collapse=""))
}

test_DNAStringSet_compaction <- function()
{
    dna <- DNAStringSet(DNA_ALPHABET)
    elementMetadata(dna) <- DataFrame(C1=dna)

    idx <- rep.int((8:6)*2L, 100L)
    dna300 <- dna[idx]
    compact_dna300 <- compact(dna300)
    checkIdentical(as.character(compact_dna300), as.character(dna300))
    ## Checking internal representation.
    checkIdentical(.poolEltLengths(compact_dna300), 3L)
    checkIdentical(.poolEltLengths(elementMetadata(compact_dna300)$C1),
                   .poolEltLengths(compact_dna300))
    checkIdentical(elementMetadata(compact_dna300)$C1@ranges,
                   compact_dna300@ranges)
}

test_DNAStringSet_showAsCell <- function()
{
    dna <- showAsCell(DNAStringSet())
    checkTrue(is(dna, "character"))
    dna <- showAsCell(DNAStringSet(DNA_ALPHABET))
    checkTrue(is(dna, "character"))
}
