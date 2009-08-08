\name{xscat}

\alias{xscat}


\title{Concatenate sequences contained in XString, XStringSet
  and/or XStringViews objects}

\description{
  This function mimics the semantic of \code{paste(..., sep="")}
  but accepts \link{XString}, \link{XStringSet} or \link{XStringViews}
  arguments and returns an \link{XString} or \link{XStringSet} object.
}

\usage{
  xscat(...)
}

\arguments{
  \item{\dots}{
    One or more character vectors (with no NAs), \link{XString},
    \link{XStringSet} or \link{XStringViews} objects.
  }
}

\value{
  An \link{XString} object if all the arguments are either \link{XString}
  objects or character strings.
  An \link{XStringSet} object otherwise.
}

\author{H. Pages}

\seealso{
  \link{XString-class},
  \link{XStringSet-class},
  \link{XStringViews-class},
  \code{\link[base]{paste}}
}

\examples{
  ## Return a BString object:
  xscat(BString("abc"), BString("EF"))
  xscat(BString("abc"), "EF")
  xscat("abc", "EF")

  ## Return a BStringSet object:
  xscat(BStringSet("abc"), "EF")

  ## Return a DNAStringSet object:
  xscat(c("t", "a"), DNAString("N"))

  ## Arguments are recycled to the length of the longest argument:
  xscat("x", LETTERS, c("3", "44", "555")) 

  ## Concatenating big XStringSet objects:
  library(drosophila2probe)
  probes <- DNAStringSet(drosophila2probe)
  mm <- complement(narrow(probes, start=13, end=13))
  left <- narrow(probes, end=12)
  right <- narrow(probes, start=14)
  xscat(left, mm, right)

  ## Collapsing an XStringSet (or XStringViews) object with a small
  ## number of elements:
  probes1000 <- as.list(probes[1:1000])
  y1 <- do.call(xscat, probes1000)
  y2 <- do.call(c, probes1000)  # slightly faster than the above
  y1 == y2  # TRUE
  ## Note that this method won't be efficient when the number of
  ## elements to collapse is big (> 10000) so we need to provide a
  ## collapse() (or xscollapse()) function in Biostrings that will
  ## be efficient at doing this. Please complain on the Bioconductor
  ## mailing list (http://bioconductor.org/docs/mailList.html) if you
  ## need this.
}

\keyword{methods}
\keyword{manip}
