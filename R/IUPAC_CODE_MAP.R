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

mergeIUPACLetters <- function(x)
{
    if (!is.character(x) || any(is.na(x)) || any(nchar(x) == 0))
        stop("'x' must be a vector of non-empty character strings")
    x <- CharacterList(strsplit(toupper(x), "", fixed=TRUE))
    yy <- unname(IUPAC_CODE_MAP[unlist(x, use.names=FALSE)])
    if (any(is.na(yy)))
        stop("some strings in 'x' contain non IUPAC letters")
    yy <- CharacterList(strsplit(yy, "", fixed=TRUE))
    y <- unstrsplit(sort(unique(IRanges:::regroupBySupergroup(yy, x))))
    names(IUPAC_CODE_MAP)[match(y, IUPAC_CODE_MAP)]
}

