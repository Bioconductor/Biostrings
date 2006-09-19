# A tricky way to create a "ByteBuffer" object at load time
# (to enable the trick, uncomment the 4 lines below + the call
# to initByteBuffer0() in .onLoad() and add bsGlobals to the NAMESPACE)

#bsGlobals <- new.env(parent=emptyenv())
#initByteBuffer0 <- function() {
#    bsGlobals$ByteBuffer0 <- new("ByteBuffer", 40)
#}

# Below is another way to create a "ByteBuffer" object before
# the package is loaded. We have to anticipate the load
# of the shared object (which normally occurs at loading
# time, and before the .onLoad() hook is called) because
# the "ByteBuffer" class needs it!
library.dynam("Biostrings", package="Biostrings")

DNA_STRING_CODEC <- DNAcodec()
RNA_STRING_CODEC <- RNAcodec()

.onLoad <- function(libname, pkgname)
{
    require("methods")
    #initByteBuffer0()
}

.onUnload <- function(libpath)
{
    library.dynam.unload("Biostrings", libpath)
}

