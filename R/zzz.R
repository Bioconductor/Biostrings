### A tricky way to create a "CharBuffer" object at load time
### (to enable the trick, uncomment the 4 lines below + the call
### to initCharBuffer0() in .onLoad() and add bsGlobals to the NAMESPACE)

#bsGlobals <- new.env(parent=emptyenv())
#initCharBuffer0 <- function() {
#    bsGlobals$CharBuffer0 <- new("CharBuffer", 40)
#}

### Below is another way to create a "CharBuffer" object before
### the package is loaded. We have to anticipate the load
### of the shared object (which normally occurs at loading
### time, and before the .onLoad() hook is called) because
### the "CharBuffer" class needs it!
library.dynam("Biostrings", package="Biostrings")

DNA_STRING_CODEC <- BStringCodec.DNA()
DNA_ALPHABET <- alphabet(DNA_STRING_CODEC)

RNA_STRING_CODEC <- BStringCodec.RNA()
RNA_ALPHABET <- alphabet(RNA_STRING_CODEC)

.onLoad <- function(libname, pkgname)
{
    require("methods")
    #initCharBuffer0()
}

.onUnload <- function(libpath)
{
    library.dynam.unload("Biostrings", libpath)
}

