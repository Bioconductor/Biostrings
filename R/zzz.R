# A tricky way to create a "bbuf" object at load time
# (to enable the trick, uncomment the 4 lines below + the call
# to initBbuf0() in .onLoad() and add bsGlobals to the NAMESPACE)

#bsGlobals <- new.env(parent=emptyenv())
#initBbuf0 <- function() {
#    bsGlobals$bbuf0 <- new("bbuf", 40)
#}

# Below is another way to create a "bbuf" object before
# the package is loaded. We have to anticipate the load
# of the shared object (which normally occurs at loading
# time, and before the .onLoad() hook is called) because
# the "bbuf" class needs it!
library.dynam("Biostrings", package="Biostrings")

DNA_STRING_CODEC <- DNAcodec()
RNA_STRING_CODEC <- RNAcodec()

.onLoad <- function(libname, pkgname) {
    require("methods")
    #initBbuf0()
}

.onUnload <- function(libpath) {
    library.dynam.unload("Biostrings", libpath)
}

