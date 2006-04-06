# ===========================================================================
# The BStringGenome class
# ---------------------------------------------------------------------------

setClass(
    "BStringGenome",
    representation(
        specie="character",
        UCSCnames="character",
        data_env="environment",  # Env where we define the data objects
        cache_env="environment"  # Private data store
    )
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Constructor-like functions and generics

# IMPORTANT: This function will NOT work on Windows if the "SaveImage"
# feature is on in the DESCRIPTION file of the client package
# (i.e. the package that creates BStringGenome objects at load-time
# by having something like
#   cel <- new("BStringGenome", UCSCnames, "CelegansGenome")
# in its zzz.R file).
# You NEED to have "SaveImage: no" in the DESCRIPTION file.
# The problem with "SaveImage: yes" is that "R CMD INSTALL" will execute
# this function BEFORE installing the "inst files", but this function needs
# the "inst files" to be installed BEFORE it can be called!
assignDataToNames <- function(x, package, subdir)
{
    names <- x@UCSCnames
    data_env <- x@data_env
    cache_env <- x@cache_env

    addCachedItem <- function(name, file)
    {
        ## Add a lazy-loading active binding thingie named 'name'
        ## containing the single object serialized in file 'file'.
        force(name)
        force(file)
        getter <- function()
        {
            if (!is.null(cache_env[[name]]))
              return(cache_env[[name]])
            found <- load(file, envir=cache_env)
            if (name != found)
              cache_env[[name]] <- get(found, envir=cache_env)
              #stop("bad data file: name of objects must match:\n",
              #     sQuote(name), " != ", sQuote(found))
            cache_env[[name]]
        }
        makeActiveBinding(name, getter, data_env)
    }

    for (name in names) {
        file <- system.file(subdir, paste(name, ".rda", sep=""), package=package)
        addCachedItem(name, file)
    }
}

setMethod("initialize", "BStringGenome",
    function(.Object, specie, UCSCnames, package, subdir)
    {
        .Object@specie <- specie
        .Object@UCSCnames <- UCSCnames
        .Object@data_env <- new.env(parent=emptyenv())
        .Object@cache_env <- new.env(parent=emptyenv())
        assignDataToNames(.Object, package, subdir)
        .Object
    }
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Standard generic methods

setMethod("show", "BStringGenome",
    function(object)
    {
        cat(object@specie, "chromosomes:\n")
        print(ls(object@data_env))
        cat("(use the '$' operator to access a given chromosome)\n") 
    }
)

setMethod("$", "BStringGenome",
    function(x, name) get(name, envir=x@data_env)
)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Other functions and generics

setGeneric("unload", function(x, what) standardGeneric("unload"))

setMethod("unload", "BStringGenome",
    function(x, what)
    {
        if (missing(what))
            what <- ls(x@cache_env) ## everything
        remove(list=what, envir=x@cache_env)
    }
)
