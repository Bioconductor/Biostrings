## Tests for coloring.R
## - make_DNA_AND_RNA_COLORED_LETTERS
## - make_AA_COLORED_LETTERS
## - update_X_palette (X is one of DNA, RNA, AA, B)

test_that("coloring works for DNA, RNA, and AA", {
    ## not a super important test
    make_DNA_AND_RNA_COLORED_LETTERS <-
        Biostrings:::make_DNA_AND_RNA_COLORED_LETTERS
    make_AA_COLORED_LETTERS <- Biostrings:::make_AA_COLORED_LETTERS

    dna_rna_expected <- c(DNA_BASES, "U", DNA_ALPHABET[-c(1:4,16:18)])
    expect_true(!any(duplicated(make_DNA_AND_RNA_COLORED_LETTERS())))
    expect_equal(sort(names(make_DNA_AND_RNA_COLORED_LETTERS())),
                 sort(dna_rna_expected))

    aa_expected <- AA_ALPHABET[-c(27:30)]
    expect_true(!any(duplicated(make_AA_COLORED_LETTERS())))
    expect_equal(sort(names(make_AA_COLORED_LETTERS())), sort(aa_expected))
})

test_that("users can update color palettes", {
    colored_letter <- \(letter, fg, bg){
        crayon::make_style(bg, bg=TRUE)(crayon::make_style(fg)(letter))
    }

    dnapalette <- get("DNA_AND_RNA_COLORED_LETTERS", envir=Biostrings:::.pkgenv)
    aapalette <- get("AA_COLORED_LETTERS", envir=Biostrings:::.pkgenv)
    bpalette <- get("BSTRING_COLORED_LETTERS", envir=Biostrings:::.pkgenv)

    origdna_palette <- Biostrings:::make_DNA_AND_RNA_COLORED_LETTERS()
    origaa_palette <- Biostrings:::make_AA_COLORED_LETTERS()
    origb_palette <- character(0L)

    ## check initialization
    expect_identical(dnapalette, origdna_palette)
    expect_identical(aapalette, origaa_palette)
    expect_identical(bpalette, origb_palette)

    ## check DNA update
    DNA_palette <- list(
      A=list(fg="blue",bg="black"),
      T=list(fg="red",bg='black'),
      G=list(fg='green',bg='black'),
      C=list(fg='yellow',bg='black')
    )
    update_DNA_palette(DNA_palette)

    dnapalette <- get("DNA_AND_RNA_COLORED_LETTERS", envir=Biostrings:::.pkgenv)
    expect_identical(dnapalette[c("A","T","G","C")],
                    c(A=colored_letter("A", "blue", "black"),
                      T=colored_letter("T", "red", "black"),
                      G=colored_letter("G", "green", "black"),
                      C=colored_letter("C", "yellow", "black")))
    update_DNA_palette()
    dnapalette <- get("DNA_AND_RNA_COLORED_LETTERS", envir=Biostrings:::.pkgenv)
    expect_identical(dnapalette, origdna_palette)

    ## Check AA update
    AA_palette <- list(
      A=list(fg="white", bg="purple"),
      B=list(fg=rgb(1,1,1), bg='orange')
    )
    update_AA_palette(AA_palette)
    aapalette <- get("AA_COLORED_LETTERS", envir=Biostrings:::.pkgenv)
    expect_identical(aapalette[c("A","B")],
                    c(A=colored_letter("A","white","purple"),
                      B=colored_letter("B", rgb(1,1,1), "orange")))
    update_AA_palette()
    aapalette <- get("AA_COLORED_LETTERS", envir=Biostrings:::.pkgenv)
    expect_identical(aapalette, origaa_palette)

    B_palette <- list(
      A=list(bg='green'),
      B=list(bg="red"),
      C=list(bg='blue'),
      D=list(fg="orange"),
      E=list(fg="yellow")
    )
    update_B_palette(B_palette)
    bpalette <- get("BSTRING_COLORED_LETTERS", envir=Biostrings:::.pkgenv)
    expect_identical(bpalette[c("A","B","C","D","E")],
                    c(A=colored_letter("A", rgb(1,1,1), "green"),
                      B=colored_letter("B", rgb(1,1,1), "red"),
                      C=colored_letter("C", rgb(1,1,1), "blue"),
                      D=crayon::make_style("orange")("D"),
                      E=crayon::make_style("yellow")("E")))

    multibyte_char_palette <- list()
    multibyte_char_palette[[rawToChar(as.raw(239L))]] <- list(fg="red")
    expect_no_condition(update_B_palette(multibyte_char_palette))

    update_B_palette()
    bpalette <- get("BSTRING_COLORED_LETTERS", envir=Biostrings:::.pkgenv)
    expect_identical(bpalette, origb_palette)

    ## sad path testing
    expect_error(update_DNA_palette(list(E=list(fg="yellow"))),
                  "Invalid codes specified.")
    expect_error(update_AA_palette(list(test=list(fg="yellow"))),
                  "Invalid codes specified.")
    expect_error(update_B_palette(list(test=list(fg="yellow"))),
                  "Invalid codes specified.")
    expect_error(update_DNA_palette(10), "should be NULL or a named list")
    expect_error(update_AA_palette(10), "should be NULL or a named list")
    expect_error(update_B_palette(10), "should be NULL or a named list")
})
