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

lettersToAmbiguity <- function(x)
{
    if (!is.character(x) || any(is.na(x)) || any(nchar(x) == 0))
        stop("'x' must be a vector of non-empty character strings")
    if (length(x) == 0)
        return(character(0))
    IUPAC_code_revmap <- names(IUPAC_CODE_MAP)
    names(IUPAC_code_revmap) <- IUPAC_CODE_MAP
    toAmbiguity <- function(xx)
    {
        all_bases <- unique(unlist(strsplit(IUPAC_CODE_MAP[xx], "", fixed=TRUE)))
        if (any(is.na(all_bases)))
            stop("some strings in 'x' contain non IUPAC letters")
        ans <- IUPAC_code_revmap[paste(sort(all_bases), collapse="")]
        names(ans) <- NULL
        ans
    }
    sapply(strsplit(toupper(x), "", fixed=TRUE), toAmbiguity)
}

