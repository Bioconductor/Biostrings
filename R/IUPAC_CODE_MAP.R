### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The IUPAC extended genetic alphabet.
###

IUPAC_CODE_MAP <- c(
    A="A",
    C="C",
    G="G",
    T="T",
    M="AC",
    R="AG",
    W="AT",
    S="CG",
    Y="CT",
    K="GT",
    V="ACG",
    H="ACT",
    D="AGT",
    B="CGT",
    N="ACGT"
)

### DNA and RNA alphabets
DNA_ALPHABET <- c(names(IUPAC_CODE_MAP), "-")
RNA_ALPHABET <- DNA_ALPHABET
RNA_ALPHABET[RNA_ALPHABET == "T"] <- "U"

