### =========================================================================
### add_colors()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###

### Placeholder, initialized in .onLoad()
DNA_AND_RNA_COLORED_LETTERS <- NULL
AA_COLORED_LETTERS <- NULL

### Return a named character vector where all the names are single letters.
### Colors for A, C, G, and T were inspired by
###   https://en.wikipedia.org/wiki/Nucleotide#Structure
### Called in .onLoad() to initialize DNA_AND_RNA_COLORED_LETTERS.
make_DNA_AND_RNA_COLORED_LETTERS <- function()
{
    ## Not sure why but the built-in white() style in the crayon package
    ## produces some kind of light grey text color. So we define a style
    ## that produces a text color that is 100% white.
    whiter <- make_style(rgb(1, 1, 1))
    dark_grey_bg <- make_style(rgb(0.5,0.5,0.5), bg=TRUE)

    ## All the IUPAC ambiguity letters minus N.
    dark_grey_bg_letters <- c("M", "R", "W", "S", "Y", "K", "V", "H", "D", "B")

    c(
        A=make_style(rgb(1, 0.5, 0.5), bg=TRUE)(black("A")),
        C=make_style(rgb(0.5, 1, 0.5), bg=TRUE)(black("C")),
        G=make_style(rgb(0.5, 1, 1), bg=TRUE)(black("G")),
        T=make_style(rgb(1, 0.8, 0.5), bg=TRUE)(black("T")),
        U=make_style(rgb(1, 1, 0.5), bg=TRUE)(black("U")),
        setNames(sprintf(dark_grey_bg(whiter("%s")), dark_grey_bg_letters),
                 dark_grey_bg_letters),
        N=make_style("grey", bg=TRUE)(whiter("N"))
    )
}

### 'x' must be a character vector.
.add_dna_and_rna_colors <- function(x)
{
    if (!isTRUE(getOption("Biostrings.coloring", default=FALSE)))
        return(x)
    ans <- vapply(x,
        function(xi) {
            xi <- safeExplode(xi)
            m <- match(xi, names(DNA_AND_RNA_COLORED_LETTERS))
            match_idx <- which(!is.na(m))
            xi[match_idx] <- DNA_AND_RNA_COLORED_LETTERS[m[match_idx]]
            paste0(xi, collapse="")
        },
        character(1),
        USE.NAMES=FALSE
    )
    x_names <- names(x)
    if (!is.null(x_names))
        names(ans) <- x_names
    ans
}

### Return a named character vector where all the names are single letters.
### Colors amino acids by similarity
### Colors groupins by
###   https://www.jalview.org/help/html/colourSchemes/zappo.html
### Called in .onLoad() to initialize AA_COLORED_LETTERS.
make_AA_COLORED_LETTERS <- function(){
    whiter <- make_style(rgb(1, 1, 1))
    dark_grey_bg <- make_style(rgb(0.5,0.5,0.5), bg=TRUE)

    ## All the IUPAC ambiguity letters minus X.
    dark_grey_bg_letters <- c("U","O","B","J","Z")

    cp <- c("#fbf8cc", "#ffcfd2", "#cfbaf0",
            "#a3c4f3", "#8eecf5", "#b9fbc0", "#f1c0e8")

    c(
        # Cysteine
        C=make_style(cp[1], bg=TRUE)(black("C")),

        # Aliphatic/hydrophobic
        A=make_style(cp[2], bg=TRUE)(black("A")),
        V=make_style(cp[2], bg=TRUE)(black("V")),
        M=make_style(cp[2], bg=TRUE)(black("M")),
        L=make_style(cp[2], bg=TRUE)(black("L")),
        I=make_style(cp[2], bg=TRUE)(black("I")),

        # Conformationally Special
        P=make_style(cp[3], bg=TRUE)(black("P")),
        G=make_style(cp[3], bg=TRUE)(black("G")),

        # Positive
        K=make_style(cp[4], bg=TRUE)(black("K")),
        R=make_style(cp[4], bg=TRUE)(black("R")),
        H=make_style(cp[4], bg=TRUE)(black("H")),

        # Hydrophilic
        N=make_style(cp[5], bg=TRUE)(black("N")),
        T=make_style(cp[5], bg=TRUE)(black("T")),
        Q=make_style(cp[5], bg=TRUE)(black("Q")),
        S=make_style(cp[5], bg=TRUE)(black("S")),

        # Aromatic
        F=make_style(cp[6], bg=TRUE)(black("F")),
        Y=make_style(cp[6], bg=TRUE)(black("Y")),
        W=make_style(cp[6], bg=TRUE)(black("W")),

        # Negative
        E=make_style(cp[7], bg=TRUE)(black("E")),
        D=make_style(cp[7], bg=TRUE)(black("D")),

        # Ambiguity
        setNames(sprintf(dark_grey_bg(whiter("%s")), dark_grey_bg_letters),
                 dark_grey_bg_letters),

        # Any code
        X=make_style("grey", bg=TRUE)(whiter("X"))
    )
}

### 'x' must be a character vector.
.add_aa_colors <- function(x)
{
    if (!isTRUE(getOption("Biostrings.coloring", default=FALSE)))
        return(x)
    ans <- vapply(x,
        function(xi) {
            xi <- safeExplode(xi)
            m <- match(xi, names(AA_COLORED_LETTERS))
            match_idx <- which(!is.na(m))
            xi[match_idx] <- AA_COLORED_LETTERS[m[match_idx]]
            paste0(xi, collapse="")
        },
        character(1),
        USE.NAMES=FALSE
    )
    x_names <- names(x)
    if (!is.null(x_names))
        names(ans) <- x_names
    ans
}

add_colors <- function(x) UseMethod("add_colors")
add_colors.default <- identity
add_colors.DNA <- add_colors.RNA <- .add_dna_and_rna_colors
add_colors.AA <- .add_aa_colors
