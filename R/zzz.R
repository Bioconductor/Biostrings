###

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
    DNA_AND_RNA_COLORED_LETTERS <<- make_DNA_AND_RNA_COLORED_LETTERS()
    AA_COLORED_LETTERS <<- make_AA_COLORED_LETTERS()
    option_name <- "Biostrings.coloring"
    if (!(option_name %in% names(.Options)))
        options(setNames(list(TRUE), option_name))

    option_name <- "Biostrings.showRaw"
    if (!(option_name %in% names(.Options)))
        options(setNames(list(FALSE), option_name))

    ## BString lookup for raw strings
    ## 256 char lookup table for 0:255 (note off by one)
    ## characters 0-31 and 127-255 are not displayable
    ## so positions 1-32 and 128-256 should be replaced
    encoding_details <- l10n_info()
    bstring_lookup <- rawToChar(as.raw(0:255), multiple=TRUE)
    invalid_chars <- c(1:32,128:256)
    if(encoding_details$`UTF-8`){
      # braille is nice if supported
      # allows for char comparisons after as.character() comparisons
      bstring_lookup[invalid_chars] <-
        as.character(parse(text=paste0("'\\U28", as.raw(95:255), "'")))
    } else if (encoding_details$MBCS){
      # use multibyte question mark if supported
      compact_unknown <- rawToChar(as.raw(c(0xef, 0xbf, 0xbd)))
      bstring_lookup[invalid_chars] <- compact_unknown
    }
    BSTRING_RAW_LOOKUP <<- bstring_lookup
}

.onUnload <- function(libpath)
{
    library.dynam.unload("Biostrings", libpath)
}

.test <- function() BiocGenerics:::testPackage("Biostrings")
