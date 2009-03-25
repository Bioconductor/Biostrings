\name{xscat}

\alias{xscat}


\title{Concatenate sequences contained in XString and/or XStringSet objects}

\description{
  This function mimics the semantic of \code{paste(..., sep="")}
  but for \link{XString} or \link{XStringSet} objects.
}

\usage{
  xscat(...)
}

\arguments{
  \item{\dots}{
    One or more \link{XString} and/or \link{XStringSet} objects.
  }
}

\seealso{
  \link{XString-class},
  \link{XStringSet-class},
  \code{\link[base]{paste}}
}

\examples{
  ## Arguments are recycled to the length of the longest argument:
  xscat("a", LETTERS, c("1", "12", "123")) 

  library(drosophila2probe)
  probes <- DNAStringSet(drosophila2probe$sequence)
  mm <- complement(narrow(probes, start=13, end=13))
  left <- narrow(probes, end=12)
  right <- narrow(probes, start=14)
  xscat(left, mm, right)

  ## Collapsing an XStringSet (or XStringViews) object with a small
  ## number of elements:
  probes1000 <- as.list(probes[1:1000])
  do.call(xscat, probes1000)
  ## Note that this method won't be efficient when the number of
  ## elements to collapse is big (> 10000) so we need to provide a
  ## collapse() (or xscollapse()) function in Biostrings that will
  ## be efficient at doing this. Please complain on the Bioconductor
  ## mailing list (http://bioconductor.org/docs/mailList.html) if you
  ## need this.
}

\keyword{methods}
\keyword{manip}
