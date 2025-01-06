###

.pkgenv <- new.env(parent=emptyenv())

.onLoad <- function(libname, pkgname)
{
    .Call2("init_DNAlkups",
          DNA_STRING_CODEC@enc_lkup, DNA_STRING_CODEC@dec_lkup,
          PACKAGE=pkgname)
    .Call2("init_RNAlkups",
          RNA_STRING_CODEC@enc_lkup, RNA_STRING_CODEC@dec_lkup,
          PACKAGE=pkgname)
    .Call2("init_AAlkups",
          AA_STRING_CODEC@enc_lkup, AA_STRING_CODEC@dec_lkup,
          PACKAGE=pkgname)
    
    assign("DNA_AND_RNA_COLORED_LETTERS", make_DNA_AND_RNA_COLORED_LETTERS(), envir=.pkgenv)
    assign("AA_COLORED_LETTERS", make_AA_COLORED_LETTERS(), envir=.pkgenv)
    assign("BSTRING_COLORED_LETTERS", character(0L), envir=.pkgenv)

    option_name <- "Biostrings.coloring"
    if (!(option_name %in% names(.Options)))
        options(setNames(list(TRUE), option_name))
}

.onUnload <- function(libpath)
{
    library.dynam.unload("Biostrings", libpath)
}

.test <- function() BiocGenerics:::testPackage("Biostrings")
