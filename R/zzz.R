###

.onLoad <- function(libname, pkgname)
{
    .Call("init_DNAlkups",
          DNA_STRING_CODEC@enc_lkup, DNA_STRING_CODEC@dec_lkup,
          PACKAGE="Biostrings")
    .Call("init_RNAlkups",
          RNA_STRING_CODEC@enc_lkup, RNA_STRING_CODEC@dec_lkup,
          PACKAGE="Biostrings")
    if (interactive() && .Platform$OS.type == "windows" &&
             .Platform$GUI == "Rgui") {
        addVigs2WinMenu("Biostrings")
    }
}

.onUnload <- function(libpath)
{
    library.dynam.unload("Biostrings", libpath)
}

