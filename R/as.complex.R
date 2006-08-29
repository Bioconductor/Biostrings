# 'base_values' must be a named vector containing the complex values
# associated to each base e.g.
#   base_values=c("A"=1+0i, "G"=0+1i, "T"=-1+0i, "C"=0-1i)
setMethod("as.complex", "DNAString",
    function(x, base_values)
    {
        z <- complex(nchar(x))
        for (i in 1:nchar(x)) {
            z[i] <- base_values[as.character(x[i])]
        }
        z
    }
)

