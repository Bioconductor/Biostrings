###

.onLoad <- function(libname, pkgname)
{
    .Call("init_DNAlkups",
          DNA_STRING_CODEC@enc_lkup, DNA_STRING_CODEC@dec_lkup,
          PACKAGE=pkgname)
    .Call("init_RNAlkups",
          RNA_STRING_CODEC@enc_lkup, RNA_STRING_CODEC@dec_lkup,
          PACKAGE=pkgname)
}

.onUnload <- function(libpath)
{
    library.dynam.unload("Biostrings", libpath)
}

