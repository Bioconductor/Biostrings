.onLoad <- function(libname, pkgname)
{
    require("methods")
}

.onUnload <- function(libpath)
{
    library.dynam.unload("Biostrings", libpath)
}

