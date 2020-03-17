### =========================================================================
### add_colors()
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###

### Placeholder, initialized in .onLoad()
DNA_AND_RNA_COLORED_LETTERS <- NULL

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

add_colors <- function(x) UseMethod("add_colors")
add_colors.default <- identity
add_colors.DNA <- add_colors.RNA <- .add_dna_and_rna_colors

