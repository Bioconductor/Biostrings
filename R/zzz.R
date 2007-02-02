.onLoad <- function(libname, pkgname)
{
    require("methods")
    #initCharBuffer0()
}

.onUnload <- function(libpath)
{
    library.dynam.unload("Biostrings", libpath)
}

